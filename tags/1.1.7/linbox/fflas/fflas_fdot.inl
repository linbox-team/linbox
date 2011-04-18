
/* fflas_fdot.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */
// Default implementation
// Specializations should be written
// to increase efficiency
template<class Field>
inline typename Field::Element
FFLAS::fdot( const Field& F, const size_t N, 
		     const typename Field::Element * x, const size_t incx,
		     const typename Field::Element * y, const size_t incy ){
	
	typename Field::Element d;
	const typename Field::Element* xi = x;
	const typename Field::Element* yi = y;
	F.init( d, 0 );
	for ( ; xi < x+N*incx; xi+=incx, yi+=incy )
		F.axpyin( d, *xi, *yi );
	return d;
}

template<>
inline FFLAS::DoubleDomain::Element
FFLAS::fdot( const DoubleDomain& , const size_t N, 
	     const DoubleDomain::Element * x, const size_t incx,
	     const DoubleDomain::Element * y, const size_t incy ){
	
	return cblas_ddot( N, x, incx, y, incy );
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
