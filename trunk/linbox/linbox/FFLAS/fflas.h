/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/FFLAS/fflas.h
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

#ifndef __FFLAS_H
#define __FFLAS_H


namespace LinBox {


#ifdef BLAS_AVAILABLE

	/* Wrapper of the function dgemm of blas .
	 * work done by J.G Dumas and C.Pernet with FFLAS routine
	 * A:= beta * A + alpha * B*C
	 */
	template <class Field>
	typename Field::Element*  Field_dgemm (const Field& F,
					       int m, int n, int k,
					       int alpha,
					       typename Field::Element * B,
					       int ldb,
					       typename Field::Element * C,
					       int ldc,
					       int beta,
					       typename Field::Element * A,
					       int lda,					       
					       int nbe=0);



	// Definition of MACRO for field_trsm function.
	enum Triangular {UPPER,LOWER};
	enum Unitary    {UNIT,NOUNIT};
	enum Side       {LEFT,RIGHT};



	/* Field wrapper of BLAS dtrsm function
	 * F is the field
	 * A and B are matrices of size m*n
	 * the size of T depends on the side  A = T^(-1).B   or   A = B.T^(-1) 
	 * ldb is the stride of B
	 * ldt is the stride of T
	 * lda is the stride of A
	 */

	template <class Field>
	void Field_trsm (const Field& F,
			 int m, int n,
			 typename Field::Element * B, int ldb,
			 typename Field::Element * T, int ldt,
			 typename Field::Element * A, int lda,
			 Triangular  tr,
			 Unitary     un,
			 Side        si);



	// in place function
	template <class Field>
	void Field_trsm (const Field& F,
			 int m, int n,
			 typename Field::Element * B, int ldb,
			 typename Field::Element * T, int ldt,
			 Triangular  tr,
			 Unitary     un,
			 Side        si) {
		
		typename Field::Element tmp[m*n];
		Field_trsm (F,m,n,B,ldb,T,ldt,tmp,n,tr,un,si);
		for (int i=0;i<m;i++) 
			for (int j=0;j<n;j++)
				F.assign(*(B+j+i*ldb),*(tmp+j+i*n));
		
	}

#endif




} // end of namespace LinBox.

#include "linbox/FFLAS/fflas.inl"
#endif
