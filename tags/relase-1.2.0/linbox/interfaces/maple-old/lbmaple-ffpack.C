/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) LinBox
 *
 *
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




//#define DEBUG 3


#include "linbox/field/modular.h"
#include "linbox/field/givaro-zpz.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "maplec.h"

using namespace LinBox;

typedef Modular<double> Field;
//typedef GivaroZpz<Std32> Field;
typedef Field::Element Element;


extern "C"
{


	ALGEB fgemm(MKernelVector kv, ALGEB* argv )
	{

		/* expect arguments in following order:
		   1- p (int prime)
		   2- m (int)
		   3- n (int)
		   4- k (int)
		   5- alpha (int)
		   6- A (matrix)
		   7- B (matrix)
		   8- beta (int)
		   9- C (matrix)
		*/

		//MaplePrintf(kv,"Modular Matrix Mutiplication using LinBox library\n");

		int p,m,n,k,alpha,beta;
		p     = MapleToInteger32(kv,(ALGEB)argv[1]);
		m     = MapleToInteger32(kv,(ALGEB)argv[2]);
		n     = MapleToInteger32(kv,(ALGEB)argv[3]);
		k     = MapleToInteger32(kv,(ALGEB)argv[4]);
		alpha = MapleToInteger32(kv,(ALGEB)argv[5]);
		beta  = MapleToInteger32(kv,(ALGEB)argv[8]);

		Field F(p);
		Element a,b;
		F.init(a,alpha);
		F.init(b,beta);
		Element *A, *B, *C;

		ALGEB Matrix;
		RTableSettings settings;

		// Check that argument 6,7,9 are float[8] rtable and return a pointer on it
		Matrix = (ALGEB)argv[6];
		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				MaplePrintf(kv,"Making a copy \n");
				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}
		A = (Element*)RTableDataBlock(kv,Matrix);
		//Ae=new Element[m*k];
		//for (int i=0;i<m*k;i++)
		//  *(Ae+i)=*(A+i);

		Matrix = (ALGEB)argv[7];
		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				MaplePrintf(kv,"Making a copy \n");
				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}
		B = (Element*)RTableDataBlock(kv,Matrix);
		//Be=new Element[k*n];
		//for (int i=0;i<k*n;i++)
		//  *(Be+i)=*(B+i);

		Matrix = (ALGEB)argv[9];
		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				MaplePrintf(kv,"Making a copy \n");
				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}
		C = (Element*)RTableDataBlock(kv,Matrix);

		//Ce=new Element[m*n];
		FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,m,n,k,a,A,k,B,n,b,C,n);

		Matrix = (ALGEB)argv[9];
		RTableGetSettings(kv,&settings,Matrix);
		switch (settings.data_type) {

		case RTABLE_FLOAT64:
			break;

		case RTABLE_INTEGER32:
			MaplePrintf(kv,"Converting the result to int[4]\n");
			int* Ce;
			Ce= (int *)RTableDataBlock(kv,Matrix);
			for (int i=0;i<m*n;i++)
				*(Ce+i)= (int) *(C+i);
			break;

		case RTABLE_DAG:
			MaplePrintf(kv,"Converting the result to maple int\n");
			ALGEB* Cee;
			Cee= (ALGEB *)RTableDataBlock(kv,Matrix);
			for (int i=0;i<m*n;i++)
				*(Cee+i)= ToMapleInteger(kv,(long)*(C+i));
			break;
		}



		return Matrix;

	}
}


extern "C" {

	ALGEB lsp(MKernelVector kv, ALGEB* argv )
	{

		/* expect arguments in following order:
		   1- p (int prime)
		   2- m (int)
		   3- n (int)
		   4- A (matrix)
		   5- P (vector)
		   6- k (method: k=1 -> LSP : k=2 -> LQUP)
		   7- Q (vector) only needed when k=2

		*/

		int p=MapleToInteger32(kv,argv[1]);
		int m,n,k;
		m=MapleToInteger32(kv,(ALGEB)argv[2]);
		n=MapleToInteger32(kv,(ALGEB)argv[3]);
		k=MapleToInteger32(kv,(ALGEB)argv[6]);

		Field F(p);


		ALGEB Matrix,OutMatrix,Perm1,Perm2;
		RTableSettings settings,oldsettings;

		Matrix = (ALGEB)argv[4];
		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				RTableGetSettings(kv,&oldsettings,Matrix);
				MaplePrintf(kv,"Making a copy \n");
				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}

		Element *TMP1= (Element*) RTableDataBlock(kv,Matrix);

		Perm1= (ALGEB)argv[5];
		RTableGetSettings(kv,&settings,Perm1);


		size_t *TMP2= (size_t*) RTableDataBlock(kv,Perm1);

		if (k == 2) {
			Perm2= (ALGEB) argv[7];
			size_t *TMP3= (size_t*) RTableDataBlock(kv,Perm2);
			FFPACK::LUdivine( F, FFLAS::FflasNonUnit, m, n, TMP1, n, TMP2, FFPACK::FfpackLQUP,TMP3);
		}
		else{
			FFPACK::LUdivine( F, FFLAS::FflasNonUnit, m, n, TMP1, n, TMP2, FFPACK::FfpackLSP);
		}



		OutMatrix = (ALGEB)argv[4];

		switch (settings.data_type) {

		case RTABLE_FLOAT64:
			OutMatrix=Matrix;
			break;

		case RTABLE_INTEGER32:
			MaplePrintf(kv,"Converting the result to int[4]\n");
			OutMatrix = RTableCopy(kv,&oldsettings,Matrix);
			break;

		case RTABLE_DAG:
			MaplePrintf(kv,"Converting the result to maple int\n");
			OutMatrix = RTableCopy(kv,&oldsettings,Matrix);
			break;
		}


		return OutMatrix;
	}

}

extern "C" {

	ALGEB rank(MKernelVector kv, ALGEB* argv )
	{
		/* expect arguments in following order:
		   1- p (int prime)
		   2- m (int)
		   3- n (int)
		   4- A (matrix)
		*/

		int p=MapleToInteger32(kv,argv[1]);
		int m,n;
		m=MapleToInteger32(kv,(ALGEB)argv[2]);
		n=MapleToInteger32(kv,(ALGEB)argv[3]);

		Field F(p);
		Element *A;

		RTableSettings settings;
		ALGEB Matrix = (ALGEB)argv[4];

		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				MaplePrintf(kv,"Making a copy \n");



				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}
		A= (Element*) RTableDataBlock(kv,Matrix);

		int rank=FFPACK::Rank(F,m,n,A,n);

		return ToMapleInteger(kv,rank);
	}
}



extern "C" {

	ALGEB determinant(MKernelVector kv, ALGEB* argv )
	{
		/* expect arguments in following order:
		   1- p (int prime)
		   2- m (int)
		   3- n (int)
		   4- A (matrix)
		*/

		int p=MapleToInteger32(kv,argv[1]);
		int m,n;
		m=MapleToInteger32(kv,(ALGEB)argv[2]);
		n=MapleToInteger32(kv,(ALGEB)argv[3]);

		Field F(p);
		Element *A;

		RTableSettings settings;
		ALGEB Matrix = (ALGEB)argv[4];

		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				MaplePrintf(kv,"Making a copy \n");
				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}
		A= (Element*) RTableDataBlock(kv,Matrix);

		Element d=FFPACK::Det(F,m,n,A,n);

		return ToMapleInteger(kv,long(d));
	}
}


extern "C" {

	ALGEB inverse(MKernelVector kv, ALGEB* argv )
	{
		/* expect arguments in following order:
		   1- p (int prime)
		   2- m (int)
		   3- A (matrix)
		   4- X (matrix)
		*/

		int p=MapleToInteger32(kv,argv[1]);
		int m;
		m=MapleToInteger32(kv,(ALGEB)argv[2]);

		Field F(p);
		Element *A,*X;

		RTableSettings settings;
		ALGEB Matrix;

		Matrix= (ALGEB)argv[3];
		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				MaplePrintf(kv,"Making a copy \n");
				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}
		A= (Element*) RTableDataBlock(kv,Matrix);


		Matrix= (ALGEB)argv[4];
		RTableGetSettings(kv,&settings,Matrix);
		if( settings.data_type != RTABLE_FLOAT64
		    || settings.storage != RTABLE_RECT
		    || !IsMapleNULL(kv,settings.index_functions) )
			{
				MaplePrintf(kv,"Making a copy \n");
				/* not the format we wanted -- make a copy */
				settings.data_type = RTABLE_FLOAT64;
				settings.storage = RTABLE_RECT;
				settings.index_functions = ToMapleNULL(kv);
				settings.foreign = FALSE;
				Matrix = RTableCopy(kv,&settings,Matrix);
			}
		X= (Element*) RTableDataBlock(kv,Matrix);

		FFPACK::Invert(F,m,A,m,X,m);

		Matrix = (ALGEB)argv[4];
		RTableGetSettings(kv,&settings,Matrix);
		switch (settings.data_type) {

		case RTABLE_FLOAT64:
			break;

		case RTABLE_INTEGER32:
			MaplePrintf(kv,"Converting the result to int[4]\n");
			int* Ae;
			Ae= (int *)RTableDataBlock(kv,Matrix);
			for (int i=0;i<m*m;i++)
				*(Ae+i)= (int) *(X+i);
			break;

		case RTABLE_DAG:
			MaplePrintf(kv,"Converting the result to maple int\n");
			ALGEB* Cee;
			Cee= (ALGEB *)RTableDataBlock(kv,Matrix);
			for (int i=0;i<m*m;i++)
				*(Cee+i)= ToMapleInteger(kv,(long)*(X+i));
			break;
		}

		return Matrix;
	}
}


