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

	/* Field wrapper of BLAS dtrsm function
	 * A := B*T^-1.
	 * F is the field
	 * A,B are matrices of size m*n
	 * T is an upper triangular matrix of size n*n with full rank.
	 * ldb is the stride of B
	 * ldt is the stride of T
	 * lda is the stride of A
	 */

	template <class Field>
	void Field_trsm (const Field& F,
			 int m, int n,
			 typename Field::Element * B, int ldb,
			 typename Field::Element * T, int ldt,
			 typename Field::Element * A, int lda);
	
	/* Same function as above but diag(T) =Id.
	 * A := B*T^-1  
	 */

	template <class Field>
	void Field_trsm_unit (const Field& F,
			      int m, int n,
			      typename Field::Element * B, int ldb,
			      typename Field::Element * T, int ldt,
			      typename Field::Element * A, int lda);


	// in place Field_trsm
	// B := B*T^-1 with T upper triangular matrix and det T != 0 .
	template <class Field>
	void Field_trsm (const Field& F,
			 int m, int n,
			 typename Field::Element * B, int ldb,
			 typename Field::Element * T, int ldt)   { Field_trsm<Field> (F,m,n,B,ldb,T,ldt,B,ldb);}

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

} // end of namespace LinBox.

#include "linbox/FFLAS/fflas.inl"
#endif
