
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

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
