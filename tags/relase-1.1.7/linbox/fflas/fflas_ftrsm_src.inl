/* Copyright (C) 2005 LinBox
 * Written by C. Pernet
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

#define Mjoin(pre, nam) my_join(pre, nam)
#define my_join(pre, nam) pre ## nam

#ifdef __FFLAS__TRANSPOSE
 #define __FFLAS__Acolinc lda
 #define __FFLAS__Arowinc 1
#else
 #define __FFLAS__Acolinc 1
 #define __FFLAS__Arowinc lda
#endif

#ifdef __FFLAS__LEFT
 #define __FFLAS__SIDE Left
 #define __FFLAS__Na M
 #define __FFLAS__Nb N
 #define __FFLAS__Mb nsplit
 #define __FFLAS__Nb2 N
 #define __FFLAS__Mb2 M-nsplit
 #define __FFLAS__Mbrest nrestsplit
 #define __FFLAS__Nbrest N
 #define __FFLAS__Mupdate M-(i+1)*nsplit
 #define __FFLAS__Nupdate N
 #define __FFLAS__Anorminc __FFLAS__Acolinc
 #define __FFLAS__Bnorminc 1
 #define __FFLAS__Bnormnext ldb
 #define __FFLAS__Bdim N
 #ifdef __FFLAS__LOW
  #define __FFLAS__Atriang A + i * nsplit * (lda + 1)			
  #define __FFLAS__Aupdate A + i * nsplit * (lda + 1) + nsplit*__FFLAS__Arowinc
  #define __FFLAS__Arest A + (__FFLAS__Na - nrestsplit) * (lda + 1)
  #define __FFLAS__Anormnext __FFLAS__Arowinc
  #define __FFLAS__Bupdate B + (i+1)*nsplit*ldb
  #define __FFLAS__Brec B + i * nsplit * ldb
  #define __FFLAS__Brest B + (M - nrestsplit) * ldb
  #define __FFLAS__A1 A
  #define __FFLAS__A2 A + nsplit * __FFLAS__Arowinc
  #define __FFLAS__A3 A + nsplit * (lda + 1)
  #define __FFLAS__B1 B
  #define __FFLAS__B2 B + nsplit * ldb
  #define __FFLAS__Normdim i
#else
  #define __FFLAS__Atriang A + (__FFLAS__Na - (i + 1) * nsplit) * (lda + 1)
  #define __FFLAS__Aupdate A + (__FFLAS__Na  - (i + 1) * nsplit) * __FFLAS__Acolinc
  #define __FFLAS__Arest A 
  #define __FFLAS__Anormnext lda + 1
  #define __FFLAS__Bupdate B
  #define __FFLAS__Brec B + (M - (i + 1) * nsplit) * ldb
  #define __FFLAS__Brest B
  #define __FFLAS__A1 A + (__FFLAS__Na - nsplit) * (lda + 1)
  #define __FFLAS__A2 A + (__FFLAS__Na - nsplit) * __FFLAS__Acolinc
  #define __FFLAS__A3 A 
  #define __FFLAS__B1 B + (M - nsplit)*ldb
  #define __FFLAS__B2 B
  #define __FFLAS__Normdim __FFLAS__Na-i-1
 #endif
#else	
 #define __FFLAS__SIDE Right
 #define __FFLAS__Na N
 #define __FFLAS__Nb nsplit
 #define __FFLAS__Mb M
 #define __FFLAS__Mb2 M
 #define __FFLAS__Nb2 N-nsplit
 #define __FFLAS__Mbrest M
 #define __FFLAS__Nbrest nrestsplit
 #define __FFLAS__Mupdate M
 #define __FFLAS__Nupdate N - (i + 1) * nsplit
 #define __FFLAS__Anorminc __FFLAS__Arowinc
 #define __FFLAS__Bnorminc ldb
 #define __FFLAS__Bnormnext 1
 #define __FFLAS__Bdim M
 #ifdef __FFLAS__UP
  #define __FFLAS__Atriang A + i * nsplit * (lda + 1)
  #define __FFLAS__Aupdate A + i * nsplit * (lda + 1) + nsplit * __FFLAS__Acolinc
  #define __FFLAS__Arest A + (__FFLAS__Na - nrestsplit) * (lda + 1)
  #define __FFLAS__Anormnext __FFLAS__Acolinc
  #define __FFLAS__Bupdate B + (i + 1) * nsplit
  #define __FFLAS__Brec B + i * nsplit
  #define __FFLAS__Brest B + (N - nrestsplit)
  #define __FFLAS__A1 A
  #define __FFLAS__A2 A + nsplit * __FFLAS__Acolinc
  #define __FFLAS__A3 A + nsplit * (lda + 1)
  #define __FFLAS__B1 B
  #define __FFLAS__B2 B + nsplit
  #define __FFLAS__Normdim i
#else
  #define __FFLAS__Atriang A + (__FFLAS__Na - (i + 1) * nsplit) * (lda + 1)
  #define __FFLAS__Aupdate A + (__FFLAS__Na - (i + 1) * nsplit) * __FFLAS__Arowinc
  #define __FFLAS__Arest A
  #define __FFLAS__Anormnext lda + 1
  #define __FFLAS__Bupdate B
  #define __FFLAS__Brec B + N - (i + 1) * nsplit
  #define __FFLAS__Brest B
  #define __FFLAS__A1 A + (__FFLAS__Na - nsplit) * (lda + 1)
  #define __FFLAS__A2 A + (__FFLAS__Na - nsplit) * __FFLAS__Arowinc
  #define __FFLAS__A3 A 
  #define __FFLAS__B1 B + N - nsplit
  #define __FFLAS__B2 B
  #define __FFLAS__Normdim __FFLAS__Na - i -1
 #endif
#endif

#ifdef __FFLAS__UP
 #define __FFLAS__UPLO Upper
#else
 #define __FFLAS__UPLO Lower
#endif

#ifdef __FFLAS__UNIT
 #define __FFLAS__DIAG Unit
#else
 #define __FFLAS__DIAG NonUnit
#endif

#ifdef __FFLAS__TRANSPOSE
 #define __FFLAS__TRANS Trans
#else
 #define __FFLAS__TRANS NoTrans
#endif

#ifdef __FFLAS__DOUBLE
 #define __FFLAS__ELEMENT double
 #define __FFLAS__DOMAIN DoubleDomain
 #define __FFLAS__BLAS_PREFIX d
#endif

#ifdef __FFLAS__FLOAT
 #define __FFLAS__ELEMENT float
 #define __FFLAS__DOMAIN FloatDomain
 #define __FFLAS__BLAS_PREFIX s
#endif

#ifdef __FFLAS__GENERIC
 #define __FFLAS__ELEMENT Element
#endif



#ifndef __FFLAS__GENERIC
template <>
class FFLAS::Mjoin(ftrsm, Mjoin(__FFLAS__SIDE, Mjoin(__FFLAS__UPLO, Mjoin(__FFLAS__TRANS, __FFLAS__DIAG))))<__FFLAS__ELEMENT>{
public:

// TRSM with delayed updates: assumes input in Zp and ensures output in Zp.
// The multiple MatMul updates (recursive sequence) are done over Z
template<class Field>
void delayed (const Field& F, const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda,
	      typename Field::Element * B, const size_t ldb,
	      const size_t nblas, size_t nbblocsblas) {
	
	static typename Field::Element Mone;
	static typename Field::Element one;

	static __FFLAS__DOMAIN D;

	F.init(one, 1.0);
	F.neg(Mone,one);

	if ( __FFLAS__Na <= nblas ){
		for (size_t i=0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				F.init( *(B + i*ldb + j), *( B + i*ldb + j));
#ifndef __FFLAS__UNIT
		typename Field::Element inv;
		typename Field::Element *  Ai = A, * Bi = B;
#ifdef __FFLAS__LEFT
#ifdef __FFLAS__UP
		Ai += __FFLAS__Acolinc;
#endif
#endif
#ifdef __FFLAS__RIGHT
#ifdef __FFLAS__LOW
		Ai += __FFLAS__Arowinc; 
#endif
#endif
		for (size_t i = 0; i < __FFLAS__Na; ++i){
			F.inv (inv, *(A + i * (lda+1)));
			fscal (F, __FFLAS__Normdim, inv, Ai, __FFLAS__Anorminc);
			fscal (F, __FFLAS__Bdim, inv, Bi, __FFLAS__Bnorminc);
			Ai += __FFLAS__Anormnext;
			Bi += __FFLAS__Bnormnext;
		}
#endif // __FFLAS__UNIT
		Mjoin(cblas_,Mjoin(__FFLAS__BLAS_PREFIX,trsm))
			(CblasRowMajor,
			 Mjoin (Cblas, __FFLAS__SIDE),
			 Mjoin (Cblas, __FFLAS__UPLO),
			 Mjoin (Cblas, __FFLAS__TRANS),
			 CblasUnit,
			 M, N, 1.0, A, lda, B, ldb );
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				F.init (*(B + i*ldb + j), *(B + i*ldb + j));
				
#ifndef __FFLAS__UNIT
		Ai = A;
#ifdef __FFLAS__LEFT
#ifdef __FFLAS__UP
		Ai += __FFLAS__Acolinc;
#endif
#endif
#ifdef __FFLAS__RIGHT
#ifdef __FFLAS__LOW
		Ai += __FFLAS__Arowinc; 
#endif
#endif
		for (size_t i = 0; i < __FFLAS__Na; ++i){
			fscal( F, __FFLAS__Normdim, *(A + i * (lda+1)) , Ai, __FFLAS__Anorminc);
			Ai += __FFLAS__Anormnext;
		}
#endif // __FFLAS__UNIT
	} else { // __FFLAS__Na <= nblas
		size_t nbblocsup = (nbblocsblas + 1) / 2;
		size_t nsplit = nbblocsup * nblas;

		this->delayed (F, __FFLAS__Mb, __FFLAS__Nb,
			       __FFLAS__A1, lda, __FFLAS__B1, ldb, nblas, nbblocsup);

#ifdef __FFLAS__RIGHT
		fgemm (D, FflasNoTrans, Mjoin (Fflas, __FFLAS__TRANS), __FFLAS__Mb2, __FFLAS__Nb2, nsplit,
		       -1.0, __FFLAS__B1, ldb, __FFLAS__A2, lda, one, __FFLAS__B2, ldb);
#else
		fgemm (D, Mjoin (Fflas, __FFLAS__TRANS), FflasNoTrans, __FFLAS__Mb2, __FFLAS__Nb2, nsplit,
		       -1.0, __FFLAS__A2, lda, __FFLAS__B1, ldb, one, __FFLAS__B2, ldb);
#endif

		this->delayed (F, __FFLAS__Mb2, __FFLAS__Nb2,
			       __FFLAS__A3, lda, __FFLAS__B2, ldb, nblas, nbblocsblas - nbblocsup);
	}
}
template <class Field>
void operator () (const Field& F, const size_t M, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  typename Field::Element * B, const size_t ldb) {
	
	if (!M || !N ) return;
		
	static typename Field::Element one, Mone;
	F.init(one, 1.0);
	F.neg(Mone, one);
	
	static __FFLAS__DOMAIN D;
	size_t nblas = TRSMBound<Field> (F);

	size_t ndel = DotProdBound (F, 0, one,
#ifdef __FFLAS__DOUBLE
				    FflasDouble);
#else
	                            FflasFloat);
#endif
        ndel = (ndel / nblas)*nblas;
	size_t nsplit = ndel;
	size_t nbblocsplit = (__FFLAS__Na-1) / nsplit;
	size_t nrestsplit = ((__FFLAS__Na-1) % nsplit) +1;
	
	for ( size_t  i = 0; i < nbblocsplit; ++i) {
		this->delayed (F, __FFLAS__Mb, __FFLAS__Nb,
			       __FFLAS__Atriang, lda, __FFLAS__Brec, ldb, nblas, nsplit / nblas);
					 
#ifdef __FFLAS__RIGHT
		fgemm (F, FflasNoTrans, Mjoin (Fflas, __FFLAS__TRANS),
		       __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, Mone,
		       __FFLAS__Brec, ldb, __FFLAS__Aupdate, lda, one, __FFLAS__Bupdate, ldb);
#else
		fgemm (F, Mjoin (Fflas, __FFLAS__TRANS),  FflasNoTrans,
		       __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, Mone,
		       __FFLAS__Aupdate, lda, __FFLAS__Brec, ldb, one, __FFLAS__Bupdate, ldb);
#endif
	}
	if (nrestsplit)
		this->delayed (F, __FFLAS__Mbrest, __FFLAS__Nbrest,
			       __FFLAS__Arest, lda, __FFLAS__Brest, ldb, nblas, nrestsplit / nblas);
}


}; //class ftrsm....

#else // __FFLAS__GENERIC

template <class Element>
class FFLAS::Mjoin(ftrsm, Mjoin(__FFLAS__SIDE, Mjoin(__FFLAS__UPLO, Mjoin(__FFLAS__TRANS, __FFLAS__DIAG)))) {
public:

template<class Field>
void operator()	(const Field& F, const size_t M, const size_t N,
		 typename Field::Element * A, const size_t lda,
		 typename Field::Element * B, const size_t ldb) {
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(one, 1.0);
	F.neg(Mone,one);
	if (__FFLAS__Na == 1){
		
#ifndef __FFLAS__UNIT					
		typename Field::Element inv;
		F.inv(inv, *A);
		FFLAS::fscal(F, __FFLAS__Bdim, inv, B, __FFLAS__Bnorminc);
#endif //__FFLAS__UNIT
	} else { // __FFLAS__Na > 1
		size_t nsplit = __FFLAS__Na >> 1;
		this->operator() (F, __FFLAS__Mb, __FFLAS__Nb, __FFLAS__A1, lda, __FFLAS__B1, ldb);

#ifdef __FFLAS__RIGHT
		fgemm (F, FflasNoTrans , Mjoin (Fflas, __FFLAS__TRANS),
		       __FFLAS__Mb2, __FFLAS__Nb2, nsplit, Mone,
		       __FFLAS__B1, ldb, __FFLAS__A2, lda, one, __FFLAS__B2, ldb);
#else
		fgemm (F, Mjoin (Fflas, __FFLAS__TRANS), FflasNoTrans,
		       __FFLAS__Mb2, __FFLAS__Nb2, nsplit, Mone,
		       __FFLAS__A2, lda, __FFLAS__B1, ldb, one, __FFLAS__B2, ldb);
#endif
		this->operator() (F, __FFLAS__Mb2, __FFLAS__Nb2, __FFLAS__A3, lda, __FFLAS__B2, ldb);
	}
}
};

#endif // __FFLAS__GENERIC


#undef __FFLAS__UPLO
#undef __FFLAS__DIAG
#undef __FFLAS__SIDE
#undef __FFLAS__TRANS
#undef __FFLAS__Na
#undef __FFLAS__Mb
#undef __FFLAS__Nb
#undef __FFLAS__Mbrest
#undef __FFLAS__Nbrest 
#undef __FFLAS__Mupdate 
#undef __FFLAS__Nupdate
#undef __FFLAS__Atriang 
#undef __FFLAS__Aupdate
#undef __FFLAS__Arest 
#undef __FFLAS__Bupdate 
#undef __FFLAS__Brec 
#undef __FFLAS__Brest
#undef __FFLAS__Bnorminc
#undef __FFLAS__Bnormnext
#undef __FFLAS__Anormnext
#undef __FFLAS__Anorminc
#undef __FFLAS__ELEMENT
#undef __FFLAS__BLAS_PREFIX
#undef __FFLAS__DOMAIN
#undef __FFLAS__A1
#undef __FFLAS__A2
#undef __FFLAS__A3
#undef __FFLAS__B1
#undef __FFLAS__B2
#undef __FFLAS__Nb2
#undef __FFLAS__Mb2
#undef __FFLAS__Bdim
#undef __FFLAS__Normdim
#undef __FFLAS__Acolinc
#undef __FFLAS__Arowinc
#undef Mjoin
#undef my_join


/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
