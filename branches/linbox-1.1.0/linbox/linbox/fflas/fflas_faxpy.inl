/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflas/fflas_faxpy.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */
namespace LinBox{
template<class Field>
inline void
FFLAS::faxpy( const Field& F, const size_t N, 
	      const typename Field::Element a,
		      const typename Field::Element * X, const size_t incX,
	      typename Field::Element * Y, const size_t incY ){
	
	const typename Field::Element * Xi = X;
	typename Field::Element * Yi=Y;
	for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
		F.axpyin( *Yi, a, *Xi );
}

template<>
inline void
FFLAS::faxpy( const DoubleDomain& D, const size_t N, 
	      const DoubleDomain::Element a,
	      const DoubleDomain::Element * x, const size_t incx,
	      DoubleDomain::Element * y, const size_t incy ){

	cblas_daxpy( N, a, x, incx, y, incy);
}

}//namespace LinBox
