/* lb-maple.C
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_lb_maple_C
#define __LINBOX_lb_maple_C

#include <fstream>
#include <iostream>
#include <lb-driver.h>
#include <lb-maple-utilities.h>
#include <linbox/util/timer.h>

extern "C"{
#include <maplec.h>
}


extern "C" {
 
	/***********************************
	 * Initializer of LinBox interface *
	 ***********************************/
	
	ALGEB lbStart (MKernelVector kv, ALGEB *argv){		
		//MaplePrintf(kv, "LinBox driver initialization...\n");	
		try {LinBoxInit();} 
		catch(lb_runtime_error &t){lbRaiseError(kv, t);}	
		lb_kv = kv;
		return ToMapleNULL(kv);
	}
	
	ALGEB lbStop (MKernelVector kv, ALGEB *argv){
		//MaplePrintf(kv, "LinBox driver termination...");	
		std::cout<<"terminating LinBox...";
		//try{LinBoxEnd();}
		//catch(lb_runtime_error &t){lbRaiseError(kv, t);}	
		std::cout<<"done\n";
		return ToMapleNULL(kv);
	}


	/**************************************
	 * Information from the LinBox driver *
	 **************************************/
	ALGEB lbDataInfo (MKernelVector kv, ALGEB *argv){
		std::ostringstream out;
		LinBoxDataInfo(out);
		size_t l=out.str().length();
		char* msg = new char[l];
		strncpy(msg, out.str().c_str(), l); 
		MaplePrintf(kv, msg);
		delete  msg;
		return ToMapleNULL(kv);
	}	


	/******************************************
	 ******************************************
	 *** Function to create LinBox's object ***
	 ******************************************
	 ******************************************/
	
	/**********************************
	 * Interface to create an element *
	 **********************************/
	ALGEB lbCreateElement (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
			return ToMapleNULL(kv);
		}		
		const DomainKey *key = &MapleToDomainKey(kv, argv[1]);	       
		try { 
			LB_GMP_SET();	
			const EltKey *k = &createElement(*key);
			LB_GMP_RESTORE();	
			return ElementKeyToMaple(kv, *k); 
		} 
		catch ( lb_runtime_error &t ) 
			{ lbRaiseError(kv, t); return ToMapleNULL(kv); }				
	}

	/********************************
	 * Interface to create a domain *
	 ********************************/
	ALGEB lbCreateDomain (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);		
		const DomainKey *key;
		try {
			if (argc < 1){
				LB_GMP_SET();	
				key = &createDomain(0);
				LB_GMP_RESTORE();	
				return DomainKeyToMaple(kv, *key);
			}
			//LinBox::integer *p = GMPMapleToLinBox(kv, argv[1]);
			LinBox::integer p;
			GMPMapleToLinBox(p, kv, argv[1]);

			if (argc == 1){
				LB_GMP_SET();	
				key = &createDomain(p);
				LB_GMP_RESTORE();	
				return DomainKeyToMaple(kv, *key);
			}			
			if (argc == 2){
				LB_GMP_SET();	
				key = &createDomain(p, MapleToString(kv, argv[2]));
				LB_GMP_RESTORE();	
				return DomainKeyToMaple(kv, *key); 
			}
			if (argc > 2){
				MapleRaiseError(kv, "wrong number of argument");
				return ToMapleNULL(kv);
			}
		}
		catch ( lb_runtime_error &t ) 
			{ lbRaiseError(kv, t); return ToMapleNULL(kv);}	
		return ToMapleNULL(kv);
	}

	/**********************************
	 * Interface to create a blackbox *
	 **********************************/
	ALGEB lbCreateBlackboxFromMatrix(MKernelVector kv, ALGEB A, const LinBox::integer &p){
		
	
		
		RTableSettings setting;
		RTableGetSettings(kv,&setting,A);
		size_t m,n;
		m = RTableUpperBound(kv, A, 1);	
		n = RTableUpperBound(kv, A, 2);

	
		std::stringstream *buffer= new std::stringstream();//std::string(buffer_data, m*n));

		//LinBox::Timer chrono;
		//chrono.start();
		if (setting.storage == RTABLE_RECT)
			DenseMatrixToBuffer(kv, A, *buffer, m, n, setting);
		else
			if (setting.storage == RTABLE_SPARSE)
				SparseMatrixToBuffer(kv, A, *buffer, m, n, setting);	
			else
				MapleRaiseError(kv, "Matrix storage must be either dense or sparse");
		//chrono.stop();

		//std::ofstream FILE("MAPLE_FILE.TXT");
		//FILE<<buffer->str()<<"\n";
		//FILE.close();

		//std::cout<<"buffering in <- : "<<chrono;
		//chrono.clear();
		//chrono.start();
		LB_GMP_SET();	
		const DomainKey   *Dkey = &createDomain(p);
		const BlackboxKey *Bkey = &createBlackbox(*Dkey, *buffer); 
		deleteDomain (*Dkey);
		LB_GMP_RESTORE();
		//chrono.stop();
		//std::cout<<"buffering out -> : "<<chrono;
		delete buffer;
		
	
		return BlackboxKeyToMaple(kv, *Bkey);
	}
	


	ALGEB lbCreateBlackbox (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		const BlackboxKey *key;
		if ((argc < 1) || (argc > 4)){
			MapleRaiseError(kv, "wrong number of arguments");
			return ToMapleNULL(kv);
		}
		try {		
			if (argc == 1){
				if (IsMapleRTable(kv, argv[1]))
					return lbCreateBlackboxFromMatrix(kv, argv[1], LinBox::integer(0));
				
			}			
			if (argc == 2){
				if (IsMapleInteger(kv, argv[1]) && IsMapleRTable(kv, argv[2])){
					//LinBox::integer *p = GMPMapleToLinBox(kv, argv[1]);
					LinBox::integer p;
					GMPMapleToLinBox(p, kv, argv[1]);
					LB_GMP_SET();
					ALGEB ret = lbCreateBlackboxFromMatrix(kv, argv[2], p); 
					LB_GMP_RESTORE();
					return ret; 
				}
				    
				if (IsMapleDomainKey(kv, argv[1]) && IsMapleString(kv, argv[2])) {
					const DomainKey *k = &MapleToDomainKey(kv, argv[1]);
					std::ifstream in(MapleToString(kv, argv[2]));
					LB_GMP_SET();
					key = &createBlackbox(*k, in);
					LB_GMP_RESTORE();
					return BlackboxKeyToMaple(kv, *key);
			
				}				
			}						
			if (argc == 3){
				if (IsMapleDomainKey(kv, argv[1])){
					const DomainKey *k = &MapleToDomainKey(kv, argv[1]);
					if (IsMapleInteger(kv, argv[2])){
						LB_GMP_SET();
						key = &createBlackbox(*k, MapleToInteger32(kv, argv[2]), MapleToInteger32(kv, argv[3]));	
						LB_GMP_RESTORE();
					}
					else {
						std::ifstream in( MapleToString(kv, argv[2]));
						LB_GMP_SET();
						key = &createBlackbox(*k, in, MapleToString(kv, argv[3]));
						LB_GMP_RESTORE();
					}
					return BlackboxKeyToMaple(kv, *key);
				}
			}			
			if (argc == 4){
				if (IsMapleDomainKey(kv, argv[1])) {
					const DomainKey *k = &MapleToDomainKey(kv, argv[1]);
					LB_GMP_SET();
					key = &createBlackbox(*k, MapleToInteger32(kv, argv[2]), MapleToInteger32(kv, argv[3]), MapleToString(kv, argv[4]));
					LB_GMP_RESTORE();
					return BlackboxKeyToMaple(kv, *key);
				}
			}
			MapleRaiseError(kv, "wrong types of arguments");
		}
		catch ( lb_runtime_error &t ) 
			{ lbRaiseError(kv, t); }
		return ToMapleNULL(kv);
	}
	
	/********************************
	 * Interface to create a vector *
	 ********************************/

	ALGEB lbCreateVectorFromVector(MKernelVector kv, ALGEB V, const LinBox::integer &p){
		
		LB_GMP_SET();
		const DomainKey *Dkey = &createDomain(p);
		LB_GMP_RESTORE();
		std::stringstream buffer;
		RTableSettings setting; 
		RTableData tmp;	
		RTableGetSettings(kv,&setting,V);
		size_t n;
		n = RTableUpperBound(kv, V, 1);		
		buffer<<n<<"\n";
		M_INT index[1];
		if (setting.data_type == RTABLE_INTEGER8)
			for (size_t i=1;i<n+1; ++i){index[0]=(M_INT)i;
				tmp=RTableSelect(kv, V, index);
				buffer<<tmp.int8<<"\n";
			}		
		if (setting.data_type == RTABLE_INTEGER16)
			for (size_t i=1;i<n+1; ++i){index[0]=(M_INT)i;
				tmp=RTableSelect(kv, V, index);
				buffer<<tmp.int16<<"\n";
			}
		if (setting.data_type == RTABLE_INTEGER32)
			for (size_t i=1;i<n+1; ++i){index[0]=(M_INT)i;
				tmp=RTableSelect(kv, V, index);
				buffer<<tmp.int32<<"\n";
			}
		if (setting.data_type == RTABLE_INTEGER64)
			for (size_t i=1;i<n+1; ++i){index[0]=(M_INT)i;
				tmp=RTableSelect(kv, V, index);
				buffer<<tmp.int64<<"\n";
			}
		if (setting.data_type == RTABLE_FLOAT32)
			for (size_t i=1;i<n+1; ++i){index[0]=(M_INT)i;
				tmp=RTableSelect(kv, V, index);
				buffer<<tmp.float32<<"\n";
			}
		if (setting.data_type == RTABLE_FLOAT64)
			for (size_t i=1;i<n+1; ++i){index[0]=(M_INT)i;
				tmp=RTableSelect(kv, V, index);
				buffer<<tmp.float64<<"\n";
			}
		if (setting.data_type == RTABLE_DAG) {
			LinBox::integer ibuf;
			for (size_t i=1;i<n+1; ++i){index[0]=(M_INT)i;
				tmp=RTableSelect(kv, V, index);
				GMPMapleToLinBox(ibuf, kv,tmp.dag);
				buffer<<ibuf<<"\n";
			}
		}
		    if ((setting.data_type == RTABLE_COMPLEX)|| (setting.data_type == RTABLE_CXDAG))
			MapleRaiseError(kv, "data type format in the matrix is not yet recognized by LinBox ");

		LB_GMP_SET();
		const VectorKey *Vkey = &createVector(*Dkey, buffer);	
		deleteDomain (*Dkey);	
		LB_GMP_RESTORE();
		
		return VectorKeyToMaple(kv, *Vkey);		       
	}
	
	ALGEB lbCreateVector (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		const VectorKey *key;
		if ((argc < 1) || (argc > 4)){
			MapleRaiseError(kv, "wrong number of arguments");
			return ToMapleNULL(kv);
		} 
		try {  
			if (argc == 1){
				if (IsMapleRTable(kv, argv[1]))
					return lbCreateVectorFromVector(kv, argv[1], LinBox::integer(0));
			}					
			if (argc == 2){
				if (IsMapleInteger(kv, argv[1]) && IsMapleRTable(kv, argv[2])){
					//LinBox::integer *p = GMPMapleToLinBox(kv, argv[1]);
					LinBox::integer p;
					GMPMapleToLinBox(p, kv, argv[1]);
					LB_GMP_SET();
					ALGEB ret = lbCreateVectorFromVector(kv, argv[2], p);
					LB_GMP_RESTORE();
					return ret;
				}
				
				if ( IsMapleDomainKey(kv, argv[1]) && IsMapleInteger(kv, argv[2])) {
					const DomainKey *k = &MapleToDomainKey(kv, argv[1]);
					LB_GMP_SET();
					key = &createVector(*k, MapleToInteger32(kv, argv[2]));	
					LB_GMP_RESTORE();
					return VectorKeyToMaple(kv, *key);
				}
			}						
			if (argc == 3){
				if ( IsMapleDomainKey(kv, argv[1]) && IsMapleInteger(kv, argv[2]) && IsMapleString(kv, argv[3])){
					const DomainKey *k = &MapleToDomainKey(kv, argv[1]);
					LB_GMP_SET();
					key = &createVector(*k, MapleToInteger32(kv, argv[2]), MapleToString(kv, argv[3]));	
					LB_GMP_RESTORE();
					return VectorKeyToMaple(kv, *key);
				}
			}
			MapleRaiseError(kv, "wrong types of arguments");
		}
		catch ( lb_runtime_error &t ) 
			{ lbRaiseError(kv, t); }
		return ToMapleNULL(kv);
	}



	/***************************
	 ***************************
	 *** Domain's Functions ****
	 ***************************
	 ***************************/	

	/******************************
	 * Interface to copy a Domain *
	 *****************************/
	ALGEB lbCopyDomain (MKernelVector kv, ALGEB *argv){
		try {
			const DomainKey *key = &MapleToDomainKey(kv, argv[1]);
			return DomainKeyToMaple(kv, copyDomain(*key));
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}
		return ToMapleNULL(kv);
	}
	
	/*****************************************************
	 * Interface to change globally the prime field type *
	 *****************************************************/
	ALGEB lbSetPrimeField (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
			return ToMapleNULL(kv);
		}
		if (!IsMapleString(kv, argv[1])){
			MapleRaiseError(kv, "String expected for 1st argument");
			return ToMapleNULL(kv);
		}
		setPrimeField(MapleToString(kv, argv[1]));
		return ToMapleNULL(kv);		    
	}
	
	/********************************************************
	 * Interface to change globally the rational field type *
	 ********************************************************/
	ALGEB lbSetRationalField (MKernelVector kv, ALGEB *argv){
	M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
			return ToMapleNULL(kv);
		}
		if (!IsMapleString(kv, argv[1])){
			MapleRaiseError(kv, "String expected for 1s argument");
			return ToMapleNULL(kv);
		}
		setRationalField(MapleToString(kv, argv[1]));
		return ToMapleNULL(kv);		    
	}

	/******************************************************
	 * Interface to change globally the integer ring type *
	 ******************************************************/
	ALGEB lbSetIntegerRing (MKernelVector kv, ALGEB *argv){
	M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
			return ToMapleNULL(kv);
		}
		if (!IsMapleString(kv, argv[1])){
			MapleRaiseError(kv, "String expected for 1s argument");
			return ToMapleNULL(kv);
		}
		setIntegerRing(MapleToString(kv, argv[1]));
		return ToMapleNULL(kv);		    
	}
	

	
	/*****************************
	 *****************************
	 *** Blackbox's Functions ****
	 *****************************
	 *****************************/	
		
	/********************************
	 * Interface to copy a blackbox *
	 ********************************/
	ALGEB lbCopyBlackbox (MKernelVector kv, ALGEB *argv){
		try {
			const BlackboxKey *key = &MapleToBlackboxKey(kv, argv[1]);
			return BlackboxKeyToMaple(kv, copyBlackbox(*key));
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}

	/************************************************
	 * Interface to get the dimension of a blackbox *
	 ************************************************/
	ALGEB lbGetBlackboxDimension (MKernelVector kv, ALGEB *argv){
		// not yet  : pb with return value ?
		return ToMapleNULL(kv);
	}

	/*******************************************
	 * Interface to write a blackbox in a file *
	 *******************************************/
	ALGEB lbWriteBlackbox (MKernelVector kv, ALGEB *argv){
		try {		
			std::ofstream os(MapleToString(kv, argv[2]));
			const BlackboxKey *key = &MapleToBlackboxKey(kv, argv[1]);
			writeBlackbox(*key, os);
			return ToMapleNULL(kv);			
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}


	/**************************************************
	 * Interface to fill a blackbox with random value *
	 **************************************************/
	ALGEB lbSetBlackboxAtRandom (MKernelVector kv, ALGEB *argv){
		try {			
			const BlackboxKey *key = &MapleToBlackboxKey(kv, argv[1]);
			setBlackboxAtRandom(*key);
			return argv[1];
		}	
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}
	

	/******************************************************
	 * Interface to rebind a blackbox over another domain *
	 ******************************************************/
	ALGEB lbRebindBlackbox (MKernelVector kv, ALGEB *argv){		
		try {		
			const DomainKey   *Dkey = &MapleToDomainKey(kv, argv[1]);;
			const BlackboxKey *Bkey = &MapleToBlackboxKey(kv, argv[2]);			
			rebindBlackbox(*Bkey, *Dkey);
			return ToMapleNULL(kv);
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}


	/**************************************************
	 * Interface to change globally the blackbox type *
	 **************************************************/
	ALGEB lbSetBlackbox (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of parameter in lbSetBlackbox");
			return ToMapleNULL(kv);
		}
		if (!IsMapleString(kv, argv[1])){
			MapleRaiseError(kv, "lbSetBlackbox expects a string as parameter");
			return ToMapleNULL(kv);
		}
		setBlackbox(MapleToString(kv, argv[1]));
		return ToMapleNULL(kv);		    
	}



	/***************************
	 ***************************
	 *** Vector's Functions ****
	 ***************************
	 ***************************/
	
	/******************************
	 * Interface to copy a vector *
	 ******************************/
	ALGEB lbCopyVector (MKernelVector kv, ALGEB *argv){
		try {
			const VectorKey *key = &MapleToVectorKey(kv, argv[1]);
			return VectorKeyToMaple(kv, copyVector(*key));
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}

	/**********************************************
	 * Interface to get the dimension of a vector *
	 **********************************************/
	ALGEB lbGetVectorDimension (MKernelVector kv, ALGEB *argv){
		try {
			const VectorKey *key = &MapleToVectorKey(kv, argv[1]);
			size_t n = getVectorDimension(*key);
			return ToMapleInteger(kv, n);
		}
		catch (lb_runtime_error &t)
			{lbRaiseError(kv, t);}
		return ToMapleNULL(kv);
	}

	/*****************************************
	 * Interface to write a vector in a file *
	 *****************************************/
	ALGEB lbWriteVector (MKernelVector kv, ALGEB *argv){
		try {		
			std::ofstream os(MapleToString(kv, argv[2]));
			const VectorKey *key = &MapleToVectorKey(kv, argv[1]);
			writeVector(*key, os);
			return ToMapleNULL(kv);
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}

	/************************************************
	 * Interface to fill a vector with random value *
	 ************************************************/
	ALGEB lbSetVectorAtRandom (MKernelVector kv, ALGEB *argv){
		try {			
			const VectorKey *key = &MapleToVectorKey(kv, argv[1]);
			setVectorAtRandom(*key);
			return argv[1];
		}	
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);	
	}

	/******************************************************
	 * Interface to rebind a vector over another domain *
	 ******************************************************/
	ALGEB lbRebindVector (MKernelVector kv, ALGEB *argv){
		try {		
			const DomainKey   *Dkey = &MapleToDomainKey(kv, argv[1]);;
			const VectorKey *Bkey = &MapleToVectorKey(kv, argv[2]);			
			rebindVector(*Bkey, *Dkey);
			return argv[2];
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}

	/************************************************
	 * Interface to change globally the vector type *
	 ************************************************/
	ALGEB lbSetVector (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
			return ToMapleNULL(kv);
		}
		if (!IsMapleString(kv, argv[1])){
			MapleRaiseError(kv, "String expected as 1st argument");
			return ToMapleNULL(kv);
		}
		setVector(MapleToString(kv, argv[1]));
		return ToMapleNULL(kv);		    
	}



	/*******************************
	 *******************************
	 *** Polynomial's Functions ****
	 *******************************
	 *******************************/

	/***********************************
	 * Interface to write a polynomial *
	 ***********************************/
	ALGEB lbWritePolynomial (MKernelVector kv, ALGEB *argv){
		try {	
			std::ofstream os(MapleToString(kv, argv[2]));
			const VectorKey *key = &MapleToVectorKey(kv, argv[1]);
			writePolynomial(*key, os);
			return ToMapleNULL(kv);
		}
       		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}


	/*********************************
	 *********************************
	 **** list of LinBox solutions ***
	 *********************************
	 *********************************/

      	/******************************************************
	 * Interface to compute the determinant of a blackbox *
	 ******************************************************/
	ALGEB lbDeterminant (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		try{				
			if (argc == 1){
				LB_GMP_SET();
				const BlackboxKey *key = &MapleToBlackboxKey(kv, argv[1]);
				const EltKey        *k = &lb_determinant(*key);
				LB_GMP_RESTORE();
				return ElementKeyToMaple(kv, *k);
			}
			if (argc == 2){	
				LB_GMP_SET();
				const EltKey      *k   = &MapleToElementKey(kv, argv[1]);
				const BlackboxKey *key = &MapleToBlackboxKey(kv, argv[2]);
				lb_determinant(*k, *key);
				LB_GMP_RESTORE();
				return ToMapleNULL(kv);
			}
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}

	/***********************************************
	 * Interface to compute the rank of a blackbox *
	 ***********************************************/
	ALGEB lbRank (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		try{			
			if (argc == 1){
				LB_GMP_SET();
				const BlackboxKey *key = &MapleToBlackboxKey(kv, argv[1]);
				size_t r = lb_rank(*key);
				LB_GMP_RESTORE();
				return ToMapleInteger(kv, r);
			}
			if (argc == 2){	
				LB_GMP_SET();
				size_t *r = (size_t*) MapleToPointer(kv, argv[1]);
				const BlackboxKey *key = &MapleToBlackboxKey(kv, argv[2]);
				lb_determinant(*r, *key);
				LB_GMP_RESTORE();
				return ToMapleNULL(kv);
			}
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}



	/*************************************************************
	 * Interface to compute the minimal polynomial of a blackbox *
	 *************************************************************/
	ALGEB lbMinpoly (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		try{		
			if (argc == 1){
				LB_GMP_SET();
				const BlackboxKey    *key = &MapleToBlackboxKey(kv, argv[1]);
				const PolynomialKey  *k   = &lb_minpoly(*key);
				LB_GMP_RESTORE();
				return PolynomialKeyToMaple(kv, *k);
			}
			if (argc == 2){	
				LB_GMP_SET();	
				const PolynomialKey  *k   = &MapleToPolynomialKey(kv, argv[1]);
				const BlackboxKey    *key = &MapleToBlackboxKey(kv, argv[2]);
				lb_minpoly(*k, *key);
				LB_GMP_RESTORE();
				return ToMapleNULL(kv);
			}
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}


	/********************************************************************
	 * Interface to compute the characteristic polynomial of a blackbox *
	 ********************************************************************/
	ALGEB lbCharpoly (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		try{		
			if (argc == 1){
				LB_GMP_SET();
				const BlackboxKey    *key = &MapleToBlackboxKey(kv, argv[1]);
				const PolynomialKey  *k   = &lb_charpoly(*key);
				LB_GMP_RESTORE();
				return PolynomialKeyToMaple(kv, *k);
			}
			if (argc == 2){	
				LB_GMP_SET();			
				const PolynomialKey  *k   = &MapleToPolynomialKey(kv, argv[1]);
				const BlackboxKey    *key = &MapleToBlackboxKey(kv, argv[2]);
				lb_charpoly(*k, *key);	
				LB_GMP_RESTORE();
				return ToMapleNULL(kv);
			}
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}


	/**************************************
	 * Interface to solve a linear system *
	 **************************************/
	ALGEB lbSolve (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		try{		
			if (argc == 2){	
				LB_GMP_SET();
				const BlackboxKey *Bkey = &MapleToBlackboxKey(kv, argv[1]);
				const VectorKey   *Vkey = &MapleToVectorKey(kv, argv[2]);
				const VectorKey   *Rkey = &lb_solve(*Bkey, *Vkey);
				LB_GMP_RESTORE();
				return VectorKeyToMaple(kv, *Rkey);
			}
			if (argc == 3){	
				LB_GMP_SET();
				const VectorKey   *Rkey = &MapleToVectorKey(kv, argv[1]);
				const BlackboxKey *Bkey = &MapleToBlackboxKey(kv, argv[2]);
				const VectorKey   *Vkey = &MapleToVectorKey(kv, argv[3]);				
				lb_solve(*Rkey, *Bkey, *Vkey);	
				LB_GMP_RESTORE();
				return ToMapleNULL(kv);
			}
		}
		catch (lb_runtime_error &t)
			{ lbRaiseError(kv, t);}	
		return ToMapleNULL(kv);
	}

	

	/*******************************************
	 *******************************************
	 **** Conversion from LinBox to Maple   ****
	 *******************************************
	 *******************************************/

	/***********************************************************
	 * API to convert a domain element to its Maple equivalent *
	 ***********************************************************/
	ALGEB lbConvertElement (MKernelVector kv, ALGEB *argv){
		try {
			const EltKey *k = &MapleToElementKey(kv, argv[1]);
			SerialElement s;
			SerializeElement(s, *k);

			if (strcmp(s.type, "integer")==0)
				return LinBoxToGMPMaple(kv, s.list.front());
			else 
				if (strcmp(s.type,"rational")==0){
					ALGEB n,d,f;
					n = LinBoxToGMPMaple(kv, s.list.front());
					d = LinBoxToGMPMaple(kv, s.list.back());
					f = EvalMapleStatement(kv,"Fraction:");
					return EvalMapleProc(kv, f, 2, n, d);
				}
				else
					MapleRaiseError(kv, "LinBox internal error (serializing element problem)");
		}
		catch (lb_runtime_error &t)
			{lbRaiseError(kv, t);}
		return ToMapleNULL(kv);
	}


	/*****************************************************
	 * API to convert a blackbox to its Maple equivalent *
	 *****************************************************/
	ALGEB lbConvertBlackbox (MKernelVector kv, ALGEB *argv){
		return ToMapleNULL(kv);
	}

	/***************************************************
	 * API to convert a Vector to its Maple equivalent *
	 ***************************************************/
	ALGEB lbConvertVector (MKernelVector kv, ALGEB *argv){
		try {
			const VectorKey *k = &MapleToVectorKey(kv, argv[1]);
			SerialVector s;
			SerializeVector(s, *k);
						
			if (strcmp(s.type, "integer")==0){
				size_t n = s.list.size();
				RTableSettings setting;
				M_INT bounds[2];bounds[0]=1;bounds[1]=n;
				RTableGetDefaults(kv, &setting);
				setting.num_dimensions=1;
				setting.subtype=RTABLE_COLUMN;
				setting.data_type=RTABLE_DAG;
				
				ALGEB vector = RTableCreate(kv, &setting, NULL, bounds);
				M_INT index[1];
				RTableData tmp;
				for (size_t i=1; i<n+1; ++i){
					index[0]=i;
					tmp.dag = LinBoxToGMPMaple(kv, s.list[i-1]);
					RTableAssign(kv, vector, index, tmp);
				}
				return vector;
			}
			if (strcmp(s.type, "rational")==0){
				size_t n = (s.list.size());
				if (n & 0x1) MapleRaiseError(kv, "LinBox internal error (serializing vector problem)");
				n = n>>1;
				RTableSettings setting;
				M_INT bounds[2];bounds[0]=1;bounds[1]=n;
				RTableGetDefaults(kv, &setting);
				setting.num_dimensions=1;
				setting.subtype=RTABLE_COLUMN;
				setting.data_type=RTABLE_DAG;
			
				ALGEB vector = RTableCreate(kv, &setting, NULL, bounds);
				ALGEB f = EvalMapleStatement(kv,"Fraction:");				
				M_INT index[1];
				RTableData tmp;
				for (size_t i=1; i<n+1; ++i){
					index[0]=i;
					tmp.dag =  EvalMapleProc(kv, f, 2, LinBoxToGMPMaple(kv, s.list[(i<<1)-2]) , LinBoxToGMPMaple(kv, s.list[(i<<1)-1]));
					RTableAssign(kv, vector, index, tmp);						    
				}
				return vector;
			}
			
			else 			
				MapleRaiseError(kv, "LinBox internal error (serializing vector problem)");
		}
		catch (lb_runtime_error &t)
			{lbRaiseError(kv, t);}
		return ToMapleNULL(kv);
	}

	/*******************************************************
	 * API to convert a Polynomial to its Maple equivalent *
	 *******************************************************/
	ALGEB lbConvertPolynomial (MKernelVector kv, ALGEB *argv){
		try {
			const PolynomialKey *k = &MapleToPolynomialKey(kv, argv[1]);
			SerialPolynomial s;
			SerializePolynomial(s, *k);
						
			if (strcmp(s.type, "integer")==0){
				size_t n = s.list.size();
				ALGEB f, listcoeff;
				// ALGEB x = EvalMapleProc(kv,EvalMapleStatement(kv,"proc(s) return s; end proc;"),1,argv[2]);
				listcoeff= MapleListAlloc(kv,n);
				for (size_t i=0; i<n; ++i)
					MapleListAssign(kv,listcoeff,i+1,LinBoxToGMPMaple(kv, s.list[i]));
				f = EvalMapleStatement(kv,"proc(l, name) local i, p; p:=0; for  i from 1 to nops(l) do  p:=p+ l[i]*name^(i-1); end do; return p; end proc;");
				return EvalMapleProc(kv, f, 2, listcoeff, argv[2]);
			}
			else 			
				MapleRaiseError(kv, "LinBox internal error (serializing polynomial problem)");
		}
		catch (lb_runtime_error &t)
			{lbRaiseError(kv, t);}
		return ToMapleNULL(kv);
	}


	
	/*******************************************
	 *******************************************
	 **** Higher level API on LinBox object ****
	 *******************************************
	 *******************************************/

	/***************************
	 * Copy of a LinBox object *
	 ***************************/
	ALGEB lbCopy (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
		}
		if (IsMapleDomainKey(kv, argv[1]))
			return lbCopyDomain(kv, argv);		
		if (IsMapleBlackboxKey(kv, argv[1]))
			return lbCopyBlackbox(kv, argv);
		if (IsMapleVectorKey(kv, argv[1]))
			return lbCopyVector(kv, argv);
	
		MapleRaiseError(kv, "LinBox object (lbDomain, lbBlackbox, lbVector) expected for 1st argument");
		return ToMapleNULL(kv);
	}

	/******************************
	 * writing of a LinBox object *
	 ******************************/
	ALGEB lbWrite (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 2){
			MapleRaiseError(kv, "wrong number of arguments");
		}
		if (!IsMapleString(kv, argv[2]))
			MapleRaiseError(kv, "Filename expected for 2nd argument");
		
		if (IsMapleBlackboxKey(kv, argv[1]))
			return lbWriteBlackbox(kv, argv);
		if (IsMapleVectorKey(kv, argv[1]))
			return lbWriteVector(kv, argv);
		if (IsMaplePolynomialKey(kv, argv[1]))
			return lbWritePolynomial(kv ,argv);

		MapleRaiseError(kv, "LinBox object (lbBlackbox, lbVector, lbPolynomial) expected for 1st argument");
		return ToMapleNULL(kv);
	}

	/****************************************
	 * get the dimension of a LinBox object *
	 ****************************************/
	ALGEB lbDimension (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
		}
		if (IsMapleBlackboxKey(kv, argv[1]))
			return lbGetBlackboxDimension(kv, argv);
		if (IsMapleVectorKey(kv, argv[1]))
			return lbGetVectorDimension(kv, argv);
		
		MapleRaiseError(kv, "LinBox object (lbBlackbox, lbVector) expected for 1st argument");
		return ToMapleNULL(kv);
	}

	/**********************************************
	 * rebind a linbox object over another domain *
	 **********************************************/
	ALGEB lbRebind (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 2){
			MapleRaiseError(kv, "wrong number of arguments");
		}
		if (!IsMapleDomainKey(kv, argv[1]))
			MapleRaiseError(kv, "LinBox object (lbDomain) expected for 1st argument");
		
		if (IsMapleBlackboxKey(kv, argv[2]))
			return lbRebindBlackbox(kv, argv);
		if (IsMapleVectorKey(kv, argv[2]))
			return lbRebindVector(kv, argv);

		MapleRaiseError(kv, "LinBox object (lbBlackbox, lbVector) expected for 1st argument");
		return ToMapleNULL(kv);
	}

	/**********************************
	 * fill randomly a linbox object  *
	 **********************************/
	ALGEB lbRandom (MKernelVector kv, ALGEB *argv){
		M_INT argc = MapleNumArgs(kv, (ALGEB) argv);
		if (argc != 1){
			MapleRaiseError(kv, "wrong number of arguments");
		}
		if (IsMapleBlackboxKey(kv, argv[1]))
			return lbSetBlackboxAtRandom(kv, argv);
		if (IsMapleVectorKey(kv, argv[1]))
			return lbSetVectorAtRandom(kv, argv);

		MapleRaiseError(kv, "LinBox object (lbBlackbox, lbVector) expected for 1st argument");
		return ToMapleNULL(kv);
	}


} // end of extern "C"

#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
