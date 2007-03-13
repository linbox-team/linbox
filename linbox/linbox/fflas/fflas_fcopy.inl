/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* fflas/fflas_fcopy.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

template<class Field>
inline void
FFLAS::fcopy (const Field& F, const size_t N, 
	      typename Field::Element * X, const size_t incX,
	      const typename Field::Element * Y, const size_t incY ){
	
	typename Field::Element * Xi = X;
	const typename Field::Element * Yi=Y;
	for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
		F.assign(*Xi,*Yi);
}

// template<>
// inline void
// FFLAS::fcopy (const Modular<double>& F, const size_t N, 
// 	      double * X, const size_t incX,
// 	      const double * Y, const size_t incY ){
	
// 	cblas_dcopy(N,Y,incY,X,incX);
// }
