/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2011 LinBox
 * Written by BB <brice.boyer@imag.fr>
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


// #define _LB_LOG2 0.6931471807
#define _LB_LOG2 0.30102999566398119521


// typedef std::vector<LinBox::integer> Ivect ;
// typedef std::vector<double>          Fvect ;

struct ReductVect {
	std::vector<Integer>              & _v ;
	// Fvect                             & _w ;
	mutable std::vector<double>         _r ;

	ReductVect(std::vector<Integer> &v) :
		_v(v)
   	{
		_r.resize(v.size());
	} ;

	template<class Field>
	std::vector<typename Field::Element> &
	operator()(std::vector<typename Field::Element> & r, const Field & F) const
	{
		typedef typename Field::Element          Element;
		typedef std::vector<Element>               Fvect;
		typedef typename Fvect::iterator            Iter;


		//! @todo LinBox hom or magic here ?
		std::vector<Integer>::iterator j = _v.begin();

		for (Iter i = r.begin() ; i != r.end() ; ++i,++j) {
			F.init(*i,*j);
		}
		return r ;
	}

#if 0
	template<class Field>
	std::vector<typename Field::Element>::iterator &
	operator()(std::vector<typename Field::Element>::iterator & r, const Field & F) const
	{
		typedef typename Field::Element          Element;
		typedef std::vector<Element>               Fvect;
		typedef typename Fvect::iterator            Iter;

		Fvect  _r (_v.size()) ;
		//! @todo LinBox hom or magic here ?
		std::vector<Integer>::iterator j = _v.begin();
		for (Iter i = _r.begin() ; i != _r.end() ; ++i,++j) {
			F.init(*i,*j);
		}
		return r=_r.begin() ;
	}

#endif
	std::vector<double>::iterator &
	operator()(std::vector<double>::iterator & r, const LinBox::Modular<double> & F) const
	{
		// std::vector<double>  _r (_v.size()) ;
		//! @todo LinBox hom or magic here ?
		std::vector<Integer>::iterator j = _v.begin();
		typedef std::vector<double>::iterator Iter ;
		for (Iter i = _r.begin() ; i != _r.end() ; ++i,++j) {
			F.init(*i,*j);
		}
		// std::cout << _r << std::endl;
		return r=_r.begin() ;
	}

};

template<class Field>
struct ReductPoint {
	Integer & _v ;
	typedef typename Field::Element Element;
	ReductPoint(Integer &v) :
		_v(v)
   	{} ;
	Element & operator()(Element & r, const Field & F)
	{
		F.init(r,_v);
		return r ;
	}
};

namespace LinBox
{
	template<class Element> struct CRATemporaryVectorTrait<ReductVect ,Element> {
		typedef typename std::vector<double>::iterator Type_t;
	};
}


/*! Bench CRA.
 * @param n size of vector to reconstruct (>1)
 * @param m number of vectors to reconstruct
 * @param l size of the integers
 * @tparam Unsigned use >=0 random integers or not.
 * @param Data collects timings
 */
template<bool Unsigned>
int bench_cra(int n, int m, unsigned int l
	      , LinBox::PlotData<index_t> & Data)
{

	Timer tim, chrono ;
	typedef std::vector<LinBox::integer> Ivect ;
	{ /* LinBox CRA */
		tim.clear();
		typedef LinBox::Modular<double> ModularField ;
		typedef LinBox::FullMultipFixedCRA< ModularField > CRAbase ;
		size_t PrimeSize = 22;
		double logV = l*_LB_LOG2 ;
		if (!Unsigned) logV += _LB_LOG2 ;
		// std::cout << logV << std::endl;
		LinBox::RandomPrimeIterator genprime( PrimeSize );
		for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
			// create the vector to reconstruct
			Ivect V(n),R(n) ;
			for (Ivect::iterator it = V.begin() ; it != V.end() ; ++it) {
				Integer::random_lessthan<Unsigned>(*it,l) ;
			}
			LinBox::ChineseRemainder<  CRAbase >  cra( std::pair<size_t,double>(n, logV) );
			ReductVect iteration(V);
			chrono.clear(); chrono.start();
			Ivect::iterator Rit = R.begin();
			cra(Rit, iteration, genprime);
			chrono.stop();
			tim += chrono ;
			if (!std::equal(R.begin(),R.end(),V.begin())) {
				std::cerr << "*** LinBox CRA failed " << (Unsigned?"positive":"general") << " ***" << std::endl;
				std::cerr << R << std::endl << "expecting " << std::endl << V << std::endl;
			}
			else
				std::cerr << "ok" << std::endl;
		}
		std::cout << "LinBox CRA :" << tim << std::endl;
	}

	{/*  do givaro crt */
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
				Integer::random_lessthan<Unsigned>(*it,l) ;
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
else
				std::cerr << "ok" << std::endl;

		}
		std::cout << "GivCRT :" << tim << std::endl;
	}

	{/*  do givaro fixed  */
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
				Integer::random_lessthan<Unsigned>(*it,l) ;
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
			else
				std::cerr << "ok" << std::endl;

		}
		std::cout << "Giv CRT Fixed :" << tim << std::endl;
	}

	/*  do iml cra/mixed */
	// Init RNS
	for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
	}
	/*  do ntl cra */
	// Init primes
	for (size_t i = 0 ; i < (size_t) m ; ++i) { // repeat m times
	}
	return EXIT_SUCCESS ;
}

int main(int ac, char** av)
{
	static size_t m = 10 ;
	static size_t l = 10 ;
	static size_t n = 10 ;
	LinBox::PlotData<index_t>  Data(n,m);
	bench_cra<true>(n,m,l,Data);
	bench_cra<false>(n,m,l,Data);
	return EXIT_SUCCESS ;
}
