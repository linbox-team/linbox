/* lb-maple-utilities.h
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


#ifndef __LINBOX_lb_maple_utilities_H
#define __LINBOX_lb_maple_utilities_H

#include <linbox-config.h>
#include <lb-driver.h>
extern "C" {
#include <gmp.h>
#include <maplec.h>
#include "mpltable.h"
}

//#define __LB_PRINT_GC
#define __LINBOX_GC


 
extern "C" {

#ifdef __LINBOX_MAPLE_GMP_ACCESS
	extern void *  (*__gmp_allocate_func)(size_t);
	extern void *  (*__gmp_reallocate_func)(void*, size_t,size_t);
	extern void    (*__gmp_free_func)(void*, size_t);
	static void * (*LB_GMP_ALLOC)(size_t)                = __gmp_allocate_func;
	static void * (*LB_GMP_REALLOC)(void*,size_t,size_t) = __gmp_reallocate_func;
	static void   (*LB_GMP_FREE)(void*,size_t)           = __gmp_free_func;
#endif

	void LB_GMP_SET(){
#ifdef __LINBOX_MAPLE_GMP_ACCESS
	mp_set_memory_functions(NULL,NULL,NULL);
#endif
	}

	void LB_GMP_RESTORE() {
#ifdef __LINBOX_MAPLE_GMP_ACCESS
		mp_set_memory_functions(LB_GMP_ALLOC,LB_GMP_REALLOC, LB_GMP_FREE);
#endif
	}



	/*********************************
	 * Raising LinBox error in Maple *
	 *********************************/
	static void lbRaiseError(MKernelVector kv, lb_runtime_error &t){
		std::ostringstream out;
		out<<t;
		size_t l=out.str().length();
		char* msg = new char[l];
		strncpy(msg, out.str().c_str(), l); 
		MapleRaiseError(kv, msg);
		delete msg;
	}


	/******************************
	 * global LinBox Maple kernel *
	 ******************************/
	static MKernelVector lb_kv;


	/*********************************
	 * Domain Key Handling           *
	 * use maple garbage collection  *
	 *********************************/	

	/*************************************************
	 * functions for garbage collection and printing *
	 *************************************************/

	static void M_DECL MarkDomainKey      (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage marked a LinBox Domain\n");
#endif
	}
	 
	static void M_DECL DisposeDomainKey   (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage collected a LinBox Domain\n");
#endif
#ifdef __LINBOX_GC	
		const DomainKey* k = (const DomainKey*) MapleToPointer(lb_kv, key);	
		try { 
			LB_GMP_SET();
			deleteDomain(*k);
			LB_GMP_RESTORE();
		}
		catch (lb_runtime_error &t)
			{;}
#endif
	}

	static ALGEB M_DECL PrintDomainKey     (ALGEB key){
		std::ostringstream out;
		const DomainKey *k =(const DomainKey*) MapleToPointer(lb_kv, key);
		writeDomainInfo(*k, out);
		char msg[out.str().length()];
		strcpy(msg, out.str().c_str());
		return  ToMapleName(lb_kv, msg, FALSE);
	}



	/*********************************
	 * function over Maple DomainKey *
	 *********************************/

	bool IsMapleDomainKey(MKernelVector kv, ALGEB k){
		return (IsMaplePointer(kv, k) && (MaplePointerType(kv, k) == (M_INT)&DisposeDomainKey));
	}


	ALGEB DomainKeyToMaple (MKernelVector kv, const DomainKey& key){
		ALGEB val;
		val = ToMaplePointer(kv, (void*) (&key), (M_INT)&DisposeDomainKey);
		MaplePointerSetMarkFunction     (kv, val, MarkDomainKey);
		MaplePointerSetDisposeFunction  (kv, val, DisposeDomainKey);
		MaplePointerSetPrintFunction    (kv, val, PrintDomainKey);
		return val;
	}

	const DomainKey& MapleToDomainKey (MKernelVector kv, ALGEB k){
		if (!IsMapleDomainKey(kv,k))
			MapleRaiseError(kv, "wrong type argument, must be a domain key");
		const DomainKey * val = (const DomainKey*) MapleToPointer(kv, k);
		return *val;
	}




	/*********************************
	 * Blackbox Key Handling         *
	 * use maple garbage collection  *
	 *********************************/

	/*************************************************
	 * functions for garbage collection and printing *
	 *************************************************/

	static void M_DECL MarkBlackboxKey      (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage marked a LinBox Blackbox\n");
#endif		
	}

	static void M_DECL DisposeBlackboxKey   (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage collected a LinBox Blackbox\n");
#endif
#ifdef __LINBOX_GC	
		//printf("collect blackbox");
		const BlackboxKey* k = (const BlackboxKey*) MapleToPointer(lb_kv, key);
		//printf("  %d\n",*k);
		try { 
			LB_GMP_SET();
			deleteBlackbox(*k);
			LB_GMP_RESTORE();
		}
		catch (lb_runtime_error &t)
			{;}
#endif
	}

	static ALGEB M_DECL PrintBlackboxKey     (ALGEB key){
		std::ostringstream out;
		const BlackboxKey *k =(const BlackboxKey*) MapleToPointer(lb_kv, key);
		if (0) {
			try {
				writeBlackbox(*k, out);
			}
			catch (lb_runtime_error &t)
				{lbRaiseError(lb_kv, t);}
		}
		else 
			writeBlackboxInfo(*k, out);
		char msg[out.str().length()];
		strcpy(msg, out.str().c_str());
		return  ToMapleName(lb_kv, msg, FALSE);
	}


	/***************************************************************
	 * conversion between DomainKey and corresponding Maple object *
	 ***************************************************************/

	bool IsMapleBlackboxKey(MKernelVector kv, ALGEB k){
		return (IsMaplePointer(kv, k) && (MaplePointerType(kv, k) == (M_INT)&DisposeBlackboxKey));
	}


	ALGEB BlackboxKeyToMaple (MKernelVector kv, const BlackboxKey& key){
		ALGEB val;
		val = ToMaplePointer(kv, (void*) (&key), (M_INT)&DisposeBlackboxKey);
		MaplePointerSetMarkFunction     (kv, val, MarkBlackboxKey);
		MaplePointerSetDisposeFunction  (kv, val, DisposeBlackboxKey);
		MaplePointerSetPrintFunction    (kv, val, PrintBlackboxKey);
		return val;
	}

	const BlackboxKey& MapleToBlackboxKey (MKernelVector kv, ALGEB k){
		if (!IsMapleBlackboxKey(kv,k))
			MapleRaiseError(kv, "wrong type argument, must be a blackbox key");
		const BlackboxKey * val = (const BlackboxKey*) MapleToPointer(kv, k);
		return *val;
	}




	/*********************************
	 * Vector Key Handling         *
	 * use maple garbage collection  *
	 *********************************/

	/*************************************************
	 * functions for garbage collection and printing *
	 *************************************************/

	static void M_DECL MarkVectorKey      (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage marked a LinBox Vector\n");
#endif
	}

	static void M_DECL DisposeVectorKey   (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage collected a LinBox Vector\n");
#endif
#ifdef __LINBOX_GC	
		const VectorKey* k = (const VectorKey*) MapleToPointer(lb_kv, key);
		try {
			LB_GMP_SET();
			deleteVector(*k);
			LB_GMP_RESTORE();
		}
		catch (lb_runtime_error &t)
			{;}
#endif
	}

	static ALGEB M_DECL PrintVectorKey     (ALGEB key){
		std::ostringstream out;
		const VectorKey *k =(const VectorKey*) MapleToPointer(lb_kv, key);
		if (false){
			try{
				writeVector(*k, out);
			}
			catch(lb_runtime_error &t)
				{lbRaiseError(lb_kv, t);}
		}
		else
			writeVectorInfo(*k, out);
		char msg[out.str().length()];
		strcpy(msg, out.str().c_str());
		return  ToMapleName(lb_kv, msg, FALSE);
	}


	/***************************************************************
	 * conversion between DomainKey and corresponding Maple object *
	 ***************************************************************/

	bool IsMapleVectorKey(MKernelVector kv, ALGEB k){
		return (IsMaplePointer(kv, k) && (MaplePointerType(kv, k) == (M_INT)&DisposeVectorKey));
	}

	ALGEB VectorKeyToMaple (MKernelVector kv, const VectorKey& key){
		ALGEB val;
		val = ToMaplePointer(kv, (void*) (&key), (M_INT)&DisposeVectorKey);
		MaplePointerSetMarkFunction     (kv, val, MarkVectorKey);
		MaplePointerSetDisposeFunction  (kv, val, DisposeVectorKey);
		MaplePointerSetPrintFunction    (kv, val, PrintVectorKey);
		return val;
	}

	const VectorKey& MapleToVectorKey (MKernelVector kv, ALGEB k){
		if (!IsMapleVectorKey(kv,k))
			MapleRaiseError(kv, "wrong type argument, must be a vector key");
		const VectorKey * val = (const VectorKey*) MapleToPointer(kv, k);
		return *val;
	}



	/*********************************
	 * Poylynomial Key Handling      *
	 * use maple garbage collection  *
	 *********************************/
	/*************************************************
	 * functions for garbage collection and printing *
	 *************************************************/

	static void M_DECL MarkPolynomialKey      (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage marked a LinBox Polynomial\n");
#endif
	}

	static void M_DECL DisposePolynomialKey   (ALGEB key){
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage collected a LinBox Polynomial\n");
#endif
#ifdef __LINBOX_GC		
		const VectorKey* k = (const VectorKey*) MapleToPointer(lb_kv, key);
		try {
			LB_GMP_SET();
			deleteVector(*k);
			LB_GMP_RESTORE();
		}
		catch (lb_runtime_error &t)
			{;}
#endif
	}

	static ALGEB M_DECL PrintPolynomialKey   (ALGEB key){
		std::ostringstream out;
		const VectorKey *k =(const VectorKey*) MapleToPointer(lb_kv, key);
		try {
			writePolynomial(*k, out);
		}
		catch(lb_runtime_error &t)
			{lbRaiseError(lb_kv ,t);}
		char msg[out.str().length()];
		strcpy(msg, out.str().c_str());
		return  ToMapleName(lb_kv, msg, FALSE);
	}


	/*******************************************************************
	 * conversion between PolynomialKey and corresponding Maple object *
	 *******************************************************************/
		
	bool IsMaplePolynomialKey(MKernelVector kv, ALGEB k){
		return (IsMaplePointer(kv, k) && (MaplePointerType(kv, k) == (M_INT)&DisposePolynomialKey));
	}

	ALGEB PolynomialKeyToMaple (MKernelVector kv, const VectorKey& key){
		ALGEB val;
		val = ToMaplePointer(kv, (void*) (&key), (M_INT)&DisposePolynomialKey);
		MaplePointerSetMarkFunction     (kv, val, MarkPolynomialKey);
		MaplePointerSetDisposeFunction  (kv, val, DisposePolynomialKey);
		MaplePointerSetPrintFunction    (kv, val, PrintPolynomialKey);
		return val;
	}

	const VectorKey& MapleToPolynomialKey (MKernelVector kv, ALGEB k){
		if (!IsMaplePolynomialKey(kv,k))
			MapleRaiseError(kv, "wrong type argument, must be a polynomial key");
		const VectorKey * val = (const VectorKey*) MapleToPointer(kv, k);
		return *val;
	}


	/*********************************
	 * Element Key Handling          *
	 * use maple garbage collection  *
	 *********************************/

	/*************************************************
	 * functions for garbage collection and printing *
	 *************************************************/

	static void M_DECL MarkElementKey      (ALGEB key){
#ifdef __LB_PRINT_GC	
		MaplePrintf(lb_kv, "Maple Garbage marked a LinBox Element\n");
#endif
	}

	static void M_DECL DisposeElementKey   (ALGEB key){	
#ifdef __LB_PRINT_GC
		MaplePrintf(lb_kv, "Maple Garbage collected a LinBox Element\n");
#endif
#ifdef __LINBOX_GC	
		const EltKey* k = (const EltKey*) MapleToPointer(lb_kv, key);
		try {
			LB_GMP_SET();
			deleteElement(*k);		
			LB_GMP_RESTORE();
		}
		catch (lb_runtime_error &t)
			{;}
#endif
	}

	static ALGEB M_DECL PrintElementKey   (ALGEB key){
		std::ostringstream out;
		const EltKey *k =(const EltKey*) MapleToPointer(lb_kv, key);
		try{
			writeElement(*k, out);
		}
		catch (lb_runtime_error &t)
			{lbRaiseError(lb_kv, t);}
		char msg[out.str().length()];
		strcpy(msg, out.str().c_str());
		return  ToMapleName(lb_kv, msg, FALSE);
	}


	/****************************************************************
	 * conversion between EkementKey and corresponding Maple object *
	 ****************************************************************/

	bool IsMapleElementKey(MKernelVector kv, ALGEB k){
		return (IsMaplePointer(kv, k) && (MaplePointerType(kv, k) == (M_INT)&DisposeElementKey));
	}

	ALGEB ElementKeyToMaple (MKernelVector kv, const EltKey& key){
		ALGEB val;
		val = ToMaplePointer(kv, (void*) (&key), (M_INT)&DisposeElementKey);
		MaplePointerSetMarkFunction     (kv, val, MarkElementKey);
		MaplePointerSetDisposeFunction  (kv, val, DisposeElementKey);
		MaplePointerSetPrintFunction    (kv, val, PrintElementKey);
		return val;
	}

	const EltKey& MapleToElementKey (MKernelVector kv, ALGEB k){
		if (!IsMapleElementKey(kv,k))
			MapleRaiseError(kv, "wrong type argument, must be an Element key");
		const EltKey * val = (const EltKey*) MapleToPointer(kv, k);
		return *val;
	}


	/*********************************
	 * Miscleannous functionnalities *
	 *********************************/
	  
	/*************************************************************
	 * Conversion between Maple GMP integers and LinBox integers *
	 *************************************************************/

	void GMPMapleToLinBox(LinBox::integer& x, MKernelVector kv, ALGEB p){
#ifdef  __LINBOX_MAPLE_GMP_ACCESS
		mpz_ptr ptr = MapleToGMPInteger(kv, p);	
		//LB_GMP_SET();		
		mpz_set(LinBox::SpyInteger::get_mpz(x), ptr);
		//LB_GMP_RESTORE();
	
#else
 		// convert integer to string in order to convert to gmp integer
		ALGEB f  = EvalMapleStatement(kv,"proc(n) return convert(n,string);end proc;");
		ALGEB pp = EvalMapleProc(kv, f, 1, p);
		char *ptr = MapleToString(kv, pp);		
		mpz_set_str(LinBox::SpyInteger::get_mpz(x), ptr, 10);
#endif
	}



	ALGEB LinBoxToGMPMaple(MKernelVector kv, const LinBox::integer &p){	
#ifdef __LINBOX_MAPLE_GMP_ACCESS
		return GMPIntegerToMaple(kv, LinBox::SpyInteger::get_mpz(p));
#else	
		std::string tp(p);
		char *ptr = new char[tp.size()+1];
		std::strcpy(ptr, tp.c_str());
		//mpz_get_str(ptr,10,p.get_mpz());
		ALGEB f = EvalMapleStatement(kv,"proc(n) return convert(n, decimal, 10);end proc;");
		ALGEB x = ToMapleString(kv, ptr);
		ALGEB y = EvalMapleProc(kv, f, 1, x);
		delete[] ptr;
		return y;
#endif
	}

	/*************************************************************
	 * Function to pass a Matrix or Vector through buffer        *
	 *************************************************************/
	inline size_t idx_fortran (size_t i, size_t j, size_t stride) {  return j*stride+i;}
	inline size_t idx_c      (size_t i, size_t j, size_t stride) {  return i*stride+j;}

	
	void DenseMatrixToBuffer (MKernelVector kv, ALGEB A, std::ostream& buffer, size_t m, size_t n, RTableSettings &setting) {
		 
		//RTableData tmp;	
		
		buffer<<m<<" "<<n;
		M_INT index[2];
		if (setting.data_type == RTABLE_INTEGER8){buffer<<"\n";
			INTEGER8 *data = (INTEGER8*) RTableDataBlock(kv,A);
			if (setting.order == RTABLE_FORTRAN)			
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_fortran(i,j,m)]<<"\n";
			else
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_c(i,j,n)]<<"\n";											
		}
		if (setting.data_type == RTABLE_INTEGER16){buffer<<"\n";
			INTEGER16 *data = (INTEGER16*) RTableDataBlock(kv,A);
			if (setting.order == RTABLE_FORTRAN)			
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_fortran(i,j,m)]<<"\n";
			else
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_c(i,j,n)]<<"\n";
		}
		if (setting.data_type == RTABLE_INTEGER32){buffer<<"\n";
			INTEGER32 *data = (INTEGER32*) RTableDataBlock(kv,A);		
			if (setting.order == RTABLE_FORTRAN)			
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_fortran(i,j,m)]<<"\n";
			else
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j){
						buffer<<data[idx_c(i,j,n)]<<"\n";
					}
		}
		if (setting.data_type == RTABLE_INTEGER64){buffer<<"\n";
			INTEGER64 *data = (INTEGER64*) RTableDataBlock(kv,A);
			if (setting.order == RTABLE_FORTRAN)			
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_fortran(i,j,m)]<<"\n";
			else
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_c(i,j,n)]<<"\n";
		}
		if (setting.data_type == RTABLE_FLOAT32){buffer<<"\n";
			FLOAT32 *data = (FLOAT32*) RTableDataBlock(kv,A);
			if (setting.order == RTABLE_FORTRAN)			
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_fortran(i,j,m)]<<"\n";
			else
				for (size_t i=1;i<m+1; ++i)
					for (size_t j=1;j<n+1;++j)
						buffer<<data[idx_c(i,j,n)]<<"\n";
		}
		if (setting.data_type == RTABLE_FLOAT64){buffer<<"\n";
			FLOAT64 *data = (FLOAT64*) RTableDataBlock(kv,A);
			if (setting.order == RTABLE_FORTRAN)			
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_fortran(i,j,m)]<<"\n";
			else
				for (size_t i=0;i<m; ++i)
					for (size_t j=0;j<n;++j)
						buffer<<data[idx_c(i,j,n)]<<"\n";
		}
		if (setting.data_type == RTABLE_DAG){	buffer<<" \n";		
			RTableData tmp;	
			LinBox::integer ibuf;
			for (size_t i=1;i<m+1; ++i){index[0]=(M_INT)i;
				for (size_t j=1;j<n+1;++j){index[1]=(M_INT)j;
				//buffer<<i<<" "<<j<<" ";
					tmp=RTableSelect(kv, A, index);				
					GMPMapleToLinBox(ibuf, kv,tmp.dag);
					buffer<<ibuf<<"\n";
				}
			}
			buffer<<"0 0 0 \n";
		}

		if ((setting.data_type == RTABLE_COMPLEX)|| (setting.data_type == RTABLE_CXDAG))
			MapleRaiseError(kv, "data type format in the matrix is not yet recognized by LinBox ");
		
	}

	void SparseMatrixToBuffer (MKernelVector kv, ALGEB A, std::ostream &buffer, size_t m, size_t n, RTableSettings &setting) {
		
		buffer<<m<<" "<<n<<" M \n";
		
		// special case for DAG data type
		if (setting.data_type == RTABLE_DAG){
			LinBox::integer ibuf;	
			M_INT index[2];
			for (size_t i=1;i<m+1; ++i){index[0]=(M_INT)i;
				for (size_t j=1;j<n+1;++j){
					index[1]=(M_INT)j;
					RTableData tmp=RTableSelect(kv, A, index);				
					GMPMapleToLinBox(ibuf, kv,tmp.dag);
					if (ibuf != 0)
						buffer<<i<<" "<<j<<" "<<ibuf<<"\n";
				}
			}		
		}
		else {
			M_INT numelem;
			NAG_INT *row, *col;
			row = RTableSparseIndexRow(kv, A, 1);
			col = RTableSparseIndexRow(kv, A, 2);
			numelem = RTableNumElements(kv,A);

			if (setting.data_type == RTABLE_INTEGER8){
				INTEGER8 *data = (INTEGER8*) RTableDataBlock(kv,A);
				for (M_INT i=0;i<numelem;++i){
					buffer<<row[i]<<" "<<col[i]<<" "<<data[i]<<"\n";
				}
			}
			if (setting.data_type == RTABLE_INTEGER16){
				INTEGER16 *data = (INTEGER16*) RTableDataBlock(kv,A);
				for (M_INT i=0;i<numelem;++i){
					buffer<<row[i]<<" "<<col[i]<<" "<<data[i]<<"\n";
				}
			}
			if (setting.data_type == RTABLE_INTEGER32){
				INTEGER32 *data = (INTEGER32*) RTableDataBlock(kv,A);
				for (M_INT i=0;i<numelem;++i){
					buffer<<row[i]<<" "<<col[i]<<" "<<data[i]<<"\n";
				}
			}
			if (setting.data_type == RTABLE_INTEGER64){
				INTEGER64 *data = (INTEGER64*) RTableDataBlock(kv,A);
				for (M_INT i=0;i<numelem;++i){
					buffer<<row[i]<<" "<<col[i]<<" "<<data[i]<<"\n";
				}
			}
			if (setting.data_type == RTABLE_FLOAT32){
				FLOAT32 *data = (FLOAT32*) RTableDataBlock(kv,A);
				for (M_INT i=0;i<numelem;++i){
					buffer<<row[i]<<" "<<col[i]<<" "<<data[i]<<"\n";
				}
			}
			if (setting.data_type == RTABLE_FLOAT64){
				FLOAT64 *data = (FLOAT64*) RTableDataBlock(kv,A);
				for (M_INT i=0;i<numelem;++i){
					buffer<<row[i]<<" "<<col[i]<<" "<<data[i]<<"\n";
				}
			}	
			if ((setting.data_type == RTABLE_COMPLEX)|| (setting.data_type == RTABLE_CXDAG))
				MapleRaiseError(kv, "data type format in the matrix is not yet recognized by LinBox ");
		}
		buffer<<" 0 0 0 \n";
	}
	

} // end of extern "C"



#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
