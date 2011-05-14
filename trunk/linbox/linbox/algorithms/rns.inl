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

/*! @file algorithms/rns.inl
 *  @ingroup algorithms
 * @brief <b>R</b>esidue <b>N</b>umber <b>S</b>ystem implementation.
 * NO DOC
 */

#ifndef __LINBOX_algorithms_rns_INL
#define __LINBOX_algorithms_rns_INL

#include <set>
#include "linbox/util/debug.h"

namespace LinBox
{
	/* Constructor */
	template<bool Unsigned>
	RNS<Unsigned>::RNS(unsigned long l, unsigned long ps) :
		_bit_(l),  _ps_(ps)
	{
		linbox_check(ps<30); // oupa, mais moins qu'un unsigned long
		unsigned int nb_primes = (unsigned int) std::ceil(double(l)/double(ps));
		// integer maxint = Integer::pow(2,l); XXX je veux faire ça !!!!!!
		integer maxint = pow((integer)2,(Unsigned?l:l+1));
		_primes_.resize(nb_primes);
		std::set<unsigned long> primeset ;
		size_t lg      = 0 ;
		int tries      = 0 ;
		int penalty    = 0 ;
		integer curint = 1 ;
		while(true) { // this loop cannot be infinite
			while(tries < 3) { // if we fail 3 times to insert, we don't have enough primes to sample.
				if (curint>maxint)
					break;
				RandomPrimeIterator genprimes( (unsigned int) (_ps_+penalty) );
				unsigned long p = genprimes.randomPrime() ;
				++genprimes;
				primeset.insert(p);
				if (lg < primeset.size()) {
					++lg ;
					Integer::mulin(curint,p) ;
				}
				else { // not inserted
					++ tries ;
				}
				//! @todo if log2(maxint/curint)<ps use smaller genprime.
			}
			if (curint>maxint)
				break;
			++penalty;  // no try : penalty.
			tries = 0 ; // restarting
		}
		_maxint_ = curint ;
		if (!Unsigned)
			Integer::div(_midint_,_maxint_,2);
		_size_   = lg;
		_primes_.resize(lg);
		std::set<unsigned long>::iterator pset = primeset.begin();
		Fvect::iterator pvec = _primes_.begin();
		for ( ; pvec != _primes_.end() ; ++pset, ++pvec) { // dumping set to vect
			*pvec = *pset ;
		}

		// std::cout << "primes are : " << _primes_ << std::endl;

	}

	template<bool Unsigned>
	void
	RNS<Unsigned>::initCRA()
	{
		Fvect::iterator pvec = _primes_.begin();
		_PrimeDoms_.resize( _size_ );

		typename Domains::iterator  i = _PrimeDoms_.begin();
		for(; i != _PrimeDoms_.end(); ++i, ++pvec ) {
			*i = Field( (double)*pvec  );
		}
		CRTSystem CRT( _PrimeDoms_ );
		_CRT_ = CRT  ;
		return ;
	}

	template<bool Unsigned>
	void
	RNS<Unsigned>::cra(integer & result, const std::vector<double> & residues)
	{
		Elements Moduli( _size_ );
		typename Elements::iterator e = Moduli.begin();
		std::vector<double>::const_iterator r = residues.begin();
		for (; e != Moduli.end() ; ++e, ++r) {
			*e = *r ;
		}
		_CRT_.RnsToRing( result, Moduli );
		linbox_check(result >=0);
		linbox_check(result < _maxint_ );
		if (!Unsigned)
			if (result>_midint_) {
			       	Integer::subin(result,_maxint_);
				linbox_check(result<=_midint_);
			}
		return ;
	}

	template<bool Unsigned>
	template<typename Function>
	void
	RNS<Unsigned>::cra(Ivect & result, Function & unitCRA)
	{
		std::vector<std::vector<double> > residues(_size_);
		for (int i = 0 ; i < (int)_size_ ; ++i) {
			residues[i].resize(result.size());
			unitCRA(residues[i],_PrimeDoms_[i]); // creates residue list
		}

		for (size_t i = 0 ; i < result.size() ; ++i) {
			Elements Moduli( _size_ );
			typename Elements::iterator e = Moduli.begin();
			std::vector<std::vector<double> >::iterator r = residues.begin();
			for (; e != Moduli.end() ; ++e,++r) {
				*e = (*r)[i] ;
			}
			_CRT_.RnsToRing( result[i], Moduli );
			// std::cout << result[i] << '<' << _maxint_ << std::endl;
			linbox_check(result[i] >=0);
			linbox_check(result [i]< _maxint_ );

			if (!Unsigned)
				if (result[i]>_midint_) {
				       	Integer::subin(result[i],_maxint_);
					linbox_check(result[i]<=_midint_);
				}
		}
		return ;
	}

}


namespace LinBox
{
	/* Constructor */
	template<bool Unsigned>
	RNSfixed<Unsigned>::RNSfixed(unsigned long l, unsigned long ps) :
		_bit_(l),  _ps_(ps)
	{
		linbox_check(ps<30); // oupa, mais moins qu'un unsigned long
		unsigned int nb_primes = (unsigned int) std::ceil(double(l)/double(ps));
		// integer maxint = Integer::pow(2,l); XXX je veux faire ça !!!!!!
		integer maxint = pow((integer)2,(Unsigned?l:l+1));
		// std::cout << "target max int : " << maxint << std::endl;
		_primes_.resize(nb_primes);
		std::set<unsigned long> primeset ;
		size_t lg      = 0 ;
		int tries      = 0 ;
		int penalty    = 0 ;
		integer curint = 1 ;
		while(true) { // this loop cannot be infinite
			while(tries < 3) { // if we fail 3 times to insert, we don't have enough primes to sample.
				if (curint>maxint)
					break;
				RandomPrimeIterator genprimes((unsigned int) (_ps_+penalty) );
				unsigned long p = genprimes.randomPrime() ;
				++genprimes;
				primeset.insert(p);
				if (lg < primeset.size()) {
					++lg ;
					Integer::mulin(curint,p) ;
				}
				else { // not inserted
					++ tries ;
				}
				//! @todo if log2(maxint/curint)<ps use smaller genprime.
			}
			if (curint>maxint)
				break;
			++penalty;  // no try : penalty.
			tries = 0 ; // restarting
		}
		// std::cout << "got max int : " << curint << std::endl;
		_maxint_ = curint;
		if (!Unsigned)
			Integer::div(_midint_,_maxint_,2);
		_size_ = lg;
		_primes_.resize(lg);
		std::set<unsigned long>::iterator pset = primeset.begin();
		Fvect::iterator pvec = _primes_.begin();
		for ( ; pvec != _primes_.end() ; ++pset, ++pvec) { // dumping set to vect
			*pvec = *pset ;
		}
		// std::cout << "primes are : " << _primes_ << std::endl;


	}

	template<bool Unsigned>
	void
	RNSfixed<Unsigned>::initCRA()
	{
		Fvect::iterator pvec = _primes_.begin();
		_Primes_.resize( _size_ );

		typename Prime_t::iterator  i = _Primes_.begin();
		for(; i != _Primes_.end(); ++i, ++pvec ) {
			*i = *pvec ;
		}
		CRTSystemFixed CRT(_Primes_);
		_CRT_ = CRT;
		return ;
	}

	template<bool Unsigned>
	void
	RNSfixed<Unsigned>::cra(integer & result, const std::vector<double> & residues)
	{
		Prime_t Moduli( _size_ );
		typename Prime_t::iterator e = Moduli.begin();
		std::vector<double>::const_iterator r = residues.begin();
		for (; e != Moduli.end() ; ++e, ++r) {
			*e = *r ;
		}
		_CRT_.RnsToRing( result, Moduli );
	linbox_check(result >=0);
		linbox_check(result < _maxint_ );

		if (!Unsigned)
			if (result>_midint_) {
				Integer::subin(result,_maxint_);
					linbox_check(result<=_midint_);
			}
		return ;
	}

	template<bool Unsigned>
	template<class Function>
	void
	RNSfixed<Unsigned>::cra(Ivect & result, Function & unitCRA)
	{
		std::vector<std::vector<double> > residues(_size_);
		for (size_t i = 0 ; i < _size_ ; ++i) {
			residues[i].resize(result.size());
			unitCRA(residues[i],Modular<double>(_Primes_[i]));
		}

		for (size_t i = 0 ; i < result.size() ; ++i) {
			Prime_t Moduli( _size_ );
			typename Prime_t::iterator e = Moduli.begin();
			std::vector<std::vector<double > >::iterator r = residues.begin();
			for (; e != Moduli.end() ; ++e,++r) {
				*e = (*r)[i] ;
			}
			_CRT_.RnsToRing( result[i], Moduli );
	linbox_check(result [i]>=0);
		linbox_check(result [i]< _maxint_ );

			if (!Unsigned)
				if (result[i]>_midint_){
					Integer::subin(result[i],_maxint_);
					linbox_check(result[i]<=_midint_);
				}
		}
		return ;
	}

}

#endif // __LINBOX_algorithms_rns_INL
