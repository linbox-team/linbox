
/* fflas/fflas_fger.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

template<class Field>
inline void
FFLAS::fger (const Field& F, const size_t M, const size_t N,
	     const typename Field::Element alpha, 
	     const typename Field::Element * x, const size_t incx,
	     const typename Field::Element * y, const size_t incy, 
	     typename Field::Element * A, const size_t lda){
	
	static typename Field::Element one, mone, tmp;
	F.init( one, 1UL );
	F.neg (mone, one);
	const typename Field::Element* xi=x, *yj=y;
	typename Field::Element* Ai=A;
	
	if ( M < N ){
		if ( F.areEqual( alpha, one ) )
			for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
				yj = y;
				for (size_t j = 0; j < N; ++j, yj+=incy )
					F.axpyin( *(Ai+j), *xi, *yj );
			}
		else if ( F.areEqual( alpha, mone ) )
			for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
				F.neg( tmp, *xi );
				yj = y;
				for (size_t j = 0; j < N; ++j, yj+=incy )
					F.axpyin( *(Ai+j), tmp, *yj );
			}
		else
			for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
				F.mul( tmp, alpha, *xi );
				yj = y;
				for (size_t j = 0; j < N; ++j, yj+=incy )
					F.axpyin( *(Ai+j), tmp, *yj );
			}
	} else {
		if ( F.areEqual( alpha, one ) ){
			for ( ; Ai < A+N; ++Ai, yj+=incy ){
				xi = x;
				for (size_t i = 0; i < M; ++i, xi+=incx )
					F.axpyin( *(Ai+i*lda), *xi, *yj );
			}
		}
		else if ( F.areEqual( alpha, mone ) )
			for ( ; Ai < A+N; ++Ai, yj+=incy ){
				F.neg( tmp, *yj );
				xi = x;
				for (size_t i = 0; i < M; ++i, xi+=incx )
					F.axpyin( *(Ai+i*lda), *xi, tmp );
			}
		else
			for ( ; Ai < A+N; ++Ai, yj+=incy ){
				F.mul( tmp, alpha, *yj );
				xi = x;
				for (size_t i = 0; i < M; ++i, xi+=incx )
					F.axpyin( *(Ai+i*lda), *xi, tmp );
			}
	}
			
}

template<>
inline void
FFLAS::fger( const DoubleDomain& , const size_t M, const size_t N,
		     const DoubleDomain::Element alpha, 
		     const DoubleDomain::Element * x, const size_t incx,
		     const DoubleDomain::Element * y, const size_t incy, 
		     DoubleDomain::Element * A, const size_t lda){
	
	cblas_dger( CblasRowMajor, M, N, alpha, x, incx, y, incy, A, lda );
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
