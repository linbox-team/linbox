/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2011 LinBox
 * Written by BB <brice.boyer@imag.fr>
 *
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file benchmarks/benchmark-crafixed.C
 * @ingroup benchmarks
 * @brief Benchmarking fixed CRA routines.
 * Here we make benchmarks for CRT (Chinese Remaindering Theorem/Algorithm) in
 * the following case.  Let \f$\mathbf{v}\f$ be a vector of size \f$n\f$ and
 * whose entries have at most \f$l\f$ bits (signed or unsigned).  Suppose that
 * we only know  \f$\mathbf{v} \mod p_i\f$ for many primes \f$p_i\f$. We try
 * and reconstruct  \f$\mathbf{v}\f$ from these residues.
 *
 * We benchmark for one vector or \f$m\f$ repetitions on different vectors.
 *
 * We use the implementations in LinBox, Givaro, IML and NTL (if the latter two
 * are available).
 *
 * @warning this is not a benchmark for one integer to reconstruct or a for BlasMatrix.
 */

#include "benchmarks/benchmark.h"
#include "linbox/util/debug.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-balanced.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/algorithms/rns.h"

#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-early-multip.h"
#include "linbox/integer.h"

#ifdef __LINBOX_HAVE_IML
#include "linbox/util/iml_wrapper.h"
#endif

#define _LB_LOG2 0.69314718055994530941

using Givaro::Timer ;
using LinBox::integer;

// typedef std::vector<LinBox::integer> Ivect ;
// typedef std::vector<double>          Fvect ;

struct ReductVect {
	std::vector<integer>              & _v ;

	ReductVect(std::vector<integer> &v) :
		_v(v)
   	{ } ;

	template<typename Field>
	std::vector<typename Field::Element> &
	operator()(std::vector<typename Field::Element> & r, const Field & F) const
	{
		typedef typename Field::Element          Element;
		typedef std::vector<Element>               Fvect;
		typedef typename Fvect::iterator            Iter;

		//! @todo LinBox hom or magic here ?
		std::vector<integer>::iterator j = _v.begin();

		for (Iter i =  r.begin() ; i !=  r.end() ; ++i,++j) {
			F.init(*i,*j);
		}
		return r ;
	}

};

struct ReductVectIterator ;
namespace LinBox
{
	template<class Element> struct CRATemporaryVectorTrait<ReductVectIterator ,Element> {
		typedef typename std::vector<double>::iterator Type_t;
	};
}

struct ReductVectIterator  {
	std::vector<integer>              & _v ;
	mutable std::vector<double>         _r ;

	ReductVectIterator(std::vector<integer> &v) :
		_v(v)
   	{
		_r.resize(v.size());
	} ;

	template<typename Iterator, typename Field>
	Iterator &
	operator()(Iterator & r, const Field & F) const
	{
		typedef typename Field::Element          Element;
		typedef std::vector<Element>               Fvect;
		typedef typename Fvect::iterator            Iter;

		//! @todo LinBox hom or magic here ?
		std::vector<integer>::iterator j = _v.begin();
		for (std::vector<double>::iterator i = _r.begin() ; i != _r.end() ; ++i,++j) {
			F.init(*i,*j);
		}
		return r= _r.begin() ;
	}


};

template<class Field>
struct ReductPoint {
	integer & _v ;
	typedef typename Field::Element Element;
	ReductPoint(integer &v) :
		_v(v)
   	{} ;
	Element & operator()(Element & r, const Field & F)
	{
		F.init(r,_v);
		return r ;
	}
};



/*! Bench CRA.
 * @param n size of vector to reconstruct (>1)
 * @param m number of vectors to reconstruct
 * @param l size of the integers
 * @tparam Unsigned use >=0 random integers or not.
 * @param Data collects timings
 */
template<bool Unsigned>
int bench_cra(index_t  n, index_t m, index_t l
	      , LinBox::PlotData<index_t> & Data)
{

	Timer tim, chrono ;
	typedef std::vector<LinBox::integer> Ivect ;
	{ /* LinBox CRA */
		tim.clear();
		typedef LinBox::Modular<double> ModularField ;
		typedef LinBox::FullMultipFixedCRA< ModularField > CRAbase ;
		unsigned int PrimeSize = 22;
		double logV = l*_LB_LOG2 ;
		if (!Unsigned) logV += _LB_LOG2 ;
		std::cout << "size to reconstruct : " << logV << std::endl;
		LinBox::RandomPrimeIterator genprime( PrimeSize );
		for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
			// create the vector to reconstruct
			Ivect V(n),R(n) ;
			for (Ivect::iterator it = V.begin() ; it != V.end() ; ++it) {
				integer::random_lessthan<Unsigned>(*it,l) ;
			}
#ifdef _LB_DEBUG
			for (Ivect::iterator it = V.begin() ; it != V.end() ; ++it) {
				if (naturallog(*it) > logV) {
					std::cout << *it << " too big (" << naturallog(*it) << ")" << std::endl;
				}
			}
#endif

			LinBox::ChineseRemainder<  CRAbase >  cra( std::pair<size_t,double>(n, logV) );
			ReductVectIterator iteration(V);
			chrono.clear(); chrono.start();
			Ivect::iterator Rit = R.begin();
			cra(Rit, iteration, genprime);
			chrono.stop();
			tim += chrono ;
			if (!std::equal(R.begin(),R.end(),V.begin())) {
				std::cerr << "*** LinBox CRA failed " << (Unsigned?"positive":"general") << " ***" << std::endl;
				std::cerr << R << std::endl << "expecting " << std::endl << V << std::endl;
			}
			// else
			// std::cerr << "ok" << std::endl;
		}
		std::cout << "LinBox CRA :" << tim << std::endl;
	}

	{ /* LinBox Early CRA */
		tim.clear();
		typedef LinBox::Modular<double> ModularField ;
		typedef LinBox::EarlyMultipCRA< ModularField > CRAbase ;
		unsigned int PrimeSize = 22;
		LinBox::RandomPrimeIterator genprime( PrimeSize );
		for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
			// create the vector to reconstruct
			Ivect V(n),R(n) ;
			for (Ivect::iterator it = V.begin() ; it != V.end() ; ++it) {
				integer::random_lessthan<Unsigned>(*it,l) ;
			}

			LinBox::ChineseRemainder<  CRAbase >  cra (4);
			ReductVect iteration(V);
			chrono.clear(); chrono.start();
			// Ivect::iterator Rit = R.begin();
			cra(R, iteration, genprime);
			chrono.stop();
			tim += chrono ;
			if (!std::equal(R.begin(),R.end(),V.begin())) {
				std::cerr << "*** LinBox early CRA failed " << (Unsigned?"positive":"general") << " ***" << std::endl;
				std::cerr << R << std::endl << "expecting " << std::endl << V << std::endl;
			}
			// else
			// std::cerr << "ok" << std::endl;
		}
		std::cout << "LinBox early CRA :" << tim << std::endl;
	}

	{ /*  do givaro crt */
		// Init RNS
		typedef LinBox::Modular<double> ModularField ;
		tim.clear();
		LinBox::RNS<Unsigned> rns(l) ;
		chrono.clear() ; chrono.start() ;
		rns.initCRA();
		chrono.stop();
		tim += chrono;
		for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
			Ivect V(n),R(n) ;
			for (Ivect::iterator it = V.begin() ; it != V.end() ; ++it) {
				integer::random_lessthan<Unsigned>(*it,l) ;
			}
			ReductVect iteration(V);
			chrono.clear(); chrono.start();
			rns.cra(R,iteration);
			chrono.stop();
			tim += chrono ;
			if (!std::equal(R.begin(),R.end(),V.begin())) {
				std::cerr << "*** Givaro CRT failed " << (Unsigned?"positive":"general") << "***" << std::endl;
				std::cerr << R << std::endl << "expecting " << std::endl << V << std::endl;
			}
			// else
			// std::cerr << "ok" << std::endl;

		}
		std::cout << "GivCRT :" << tim << std::endl;
	}

	{ /*  do givaro fixed  */
		// Init RNS
		typedef LinBox::Modular<double> ModularField ;
		tim.clear();
		LinBox::RNSfixed<Unsigned> rns(l) ;
		chrono.clear() ; chrono.start() ;
		rns.initCRA();
		chrono.stop();
		tim += chrono;
		for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
			Ivect V(n),R(n) ;
			for (Ivect::iterator it = V.begin() ; it != V.end() ; ++it) {
				integer::random_lessthan<Unsigned>(*it,l) ;
			}
			ReductVect iteration(V);
			chrono.clear(); chrono.start();
			rns.cra(R,iteration);
			chrono.stop();
			tim += chrono ;
			if (!std::equal(R.begin(),R.end(),V.begin())) {
				std::cerr << "*** givaro fixed failed " << (Unsigned?"positive":"general") << "***" << std::endl;
				std::cerr << R << std::endl << "expecting " << std::endl << V << std::endl;
			}
			// else
			// std::cerr << "ok" << std::endl;

		}
		std::cout << "Giv CRT Fixed :" << tim << std::endl;
	}

#if 1 /*  IML */
#ifdef __LINBOX_HAVE_IML
	{ /*  do iml cra */
		typedef LinBox::Modular<double> ModularField ;
		tim.clear();

		/* Init RNS */
		chrono.clear() ; chrono.start() ;
		long basislen = 0 ;
		IML::Double primesize;
		integer product ;
		primesize = pow(2,22);
		product = pow((integer)2,l);
		// mpz_init(maxi); mpz_init(mp_maxInter);
		// comment les trouver ?

		IML::FiniteField ** RNS = IML::findRNS(primesize,product.get_mpz(),&basislen);
		IML::FiniteField * liftbasis = RNS[0] ; // findLiftbasisSmall(n, maxi, &basislen);
		IML::FiniteField * cmbasis   = RNS[1] ; // combBasis(basislen,basis);
		mpz_t mp_prod ;
		IML::FiniteField * bdcoeffs = NULL ;
		IML::Double * Vp = IML_XMALLOC(IML::Double,n*basislen);
		if (!Unsigned) {
			mpz_init(mp_prod);
			IML::basisProd(basislen,liftbasis,mp_prod);
			bdcoeffs =  IML::repBound(basislen, liftbasis, cmbasis) ;
		}
		chrono.stop();
		tim += chrono;

		/*  loop m times */
		for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
			/*  init result */
			Ivect V(n),R(n) ;
			for (Ivect::iterator it = V.begin() ; it != V.end() ; ++it) {
				integer::random_lessthan<Unsigned>(*it,l) ;
			}
			ReductVect iteration(V);
			for (size_t j = 0 ; j < (size_t)basislen ; ++j) {
				std::vector<double> G ;
				iteration(G,ModularField((integer)liftbasis[j]));
				for (size_t k = 0 ; k < (size_t)n ; ++k)
					Vp[j+k*basislen] = G[k] ;
			}

			/*  CRA */

			// fooooooooooooooor
			if (!Unsigned) {
				for (size_t j = 0 ; j < (size_t)n ; ++j)
					IML::ChineseRemainderPos(basislen, liftbasis, cmbasis, Vp+j, R[j].get_mpz());
			}
			else {
				for (size_t j = 0 ; j < (size_t)n ; ++j)
					IML::ChineseRemainder(basislen, mp_prod, liftbasis, cmbasis, bdcoeffs, Vp+j, R[j].get_mpz()) ;

			}
			IML_XFREE(cmbasis);
			IML_XFREE(liftbasis);
			if (!Unsigned) {
				mpz_clear(mp_prod);
				IML_XFREE(bdcoeffs);
			}
			/*  END */
		}
	}
#endif // __LINBOX_HAVE_IML
#endif

	/*  do ntl cra */
	// Init primes
	for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
	}
	return EXIT_SUCCESS ;
}

int main(int ac, char** av)
{
	static index_t m = 10 ;
	static index_t l = 200 ;
	static index_t n = 10 ;
	LinBox::PlotData<index_t>  Data(n,m);
	bench_cra<true>(n,m,l,Data);
	bench_cra<false>(n,m,l,Data);
	return EXIT_SUCCESS ;
}
