/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/FFLAS/fflas.inl
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
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
 *
 */

#ifndef __FFLAS_INL
#define __FFLAS_INL

extern "C" {
#include "cblas.h"
}
#include <math.h>

#include "linbox/FFLAS/FFFMMBLAS.h"
#include "linbox/integer.h"

namespace LinBox {



// be careful , at the end of this function T is modified (it should not be modified) !!!
template <class Field>
inline void Field_trsm_unit (const Field& F,
			     int m, int n,
			     typename Field::Element * B, int ldb,
			     typename Field::Element * T, int ldt,
			     typename Field::Element * A, int lda) {
	
	LinBox::integer p;
	LinBox::integer max_double("9007199254740991");	F.characteristic(p);
	//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
	LinBox::integer norm= integer(n)*pow(integer(2),n-1)*pow(p,n);

	if (norm > max_double )
		{
			int n2= n >>1;
			int n1= n-n2;
						
			Field_trsm_unit (F,m,n1,B,ldb,T,ldt,A,lda);
			Field_dgemm (F,m,n2,n1,-1,A,lda,T+n1,ldt,1,B+n1,ldb);
			Field_trsm_unit (F,m,n2,B+n1,ldb,T+n1*ldt+n1,ldt,A+n1,lda);		       
		}
	else {		
		
		double Td[n*n];
		double Bd[m*n];
	
		MatGFq2MatDouble_Triangular (F,n,n,Td,n,T,ldt);
		MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);
		
		cblas_dtrsm(CblasRowMajor,CblasRight,CblasUpper, CblasNoTrans,CblasUnit,
			    m,n,1.0,Td,n,Bd,n);
		MatDouble2MatGFq (F,m,n,A,lda,Bd,n);
				
		
		/*
		  for (int i=0;i<m;i++)
		  F.assign(*(A+i*lda),*(B+i*ldb));
		*/
	}
		

}
		
// be careful , at the end of this function T is modified (it should not be modified) !!!
template <class Field>
inline void Field_trsm (const Field& F,
			int m, int n,
			typename Field::Element * B, int ldb,
			typename Field::Element * T, int ldt,
			typename Field::Element * A, int lda) {
	
	LinBox::integer p;
	LinBox::integer max_double("9007199254740991");	F.characteristic(p);
	//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
	LinBox::integer norm= integer(n)*pow(integer(2),n-1)*pow(p,n);

	if (norm > max_double )
		{
			int n2= n >>1;
			int n1= n-n2;
						
			Field_trsm (F,m,n1,B,ldb,T,ldt,A,lda);
			Field_dgemm (F,m,n2,n1,-1,A,lda,T+n1,ldt,1,B+n1,ldb);
			Field_trsm (F,m,n2,B+n1,ldb,T+n1*ldt+n1,ldt,A+n1,lda);		       
		}
	else {		
		
		double Td[n*n];
		double Bd[m*n];

		typedef typename Field::Element Element;	
		Element tmp[n];	 

		// normalisation of T
		for (int i=0;i<n;i++) {
			F.inv(tmp[i],*(T+i*ldt+i));
			for (int j=i;j<n;j++)
				F.mulin(*(T+j+ldt*i),tmp[i]);
		}
				
		MatGFq2MatDouble_Triangular (F,n,n,Td,n,T,ldt);
		MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);
		
		cblas_dtrsm(CblasRowMajor,CblasRight,CblasUpper, CblasNoTrans,CblasUnit,
			    m,n,1.0,Td,n,Bd,n);
		MatDouble2MatGFq (F,m,n,A,lda,Bd,n);
				
		// denormalisation of the result A
		for (int i=0;i<m;i++) 
			for (int j=0;j<n;j++)
				F.mulin(*(A+i*lda+j),tmp[j]);

		/*
		  for (int i=0;i<m;i++)
		  F.assign(*(A+i*lda),*(B+i*ldb));
		*/
	}

}


template <class Field>
inline typename Field::Element*  Field_dgemm (const Field& F,
				       int m, int n, int k,
				       int alpha,
				       typename Field::Element * B,
				       int ldb,
				       typename Field::Element * C,
				       int ldc,
				       int beta,
				       typename Field::Element * A,
				       int lda,
				       int nbe=0) 
{ return FFFMMBLAS() (F,m,n,k,alpha,B,ldb,C,ldc,beta,A,lda,nbe);}

} // end of namespace LinBox

#endif
