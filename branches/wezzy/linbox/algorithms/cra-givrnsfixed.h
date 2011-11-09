/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/cra-full-multip.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Time-stamp: <05 Apr 11 10:49:43 Jean-Guillaume.Dumas@imag.fr>
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*!@file algorithms/cra-givrnsfixed.h
 * @ingroup algorithms
 * @brief NO DOC
 */

#ifndef __LINBOX_cra_givrnsfix_H
#define __LINBOX_cra_givrnsfix_H

#include <stdlib.h>
#include <givaro/givrnsfixed.h>

namespace LinBox
{

	/*! NO DOC...
	 * @ingroup CRA
	 * @bib
	 */
	template<class Domain_Type>
	struct GivaroRnsFixedCRA : public ::Givaro::RNSsystemFixed<integer> {
		typedef Domain_Type			                      Domain;
		typedef typename Domain::Element           DomainElement;

		typedef ::Givaro::RNSsystemFixed<integer>			Father_t;
		typedef GivaroRnsFixedCRA<Domain> 		          Self_t;

		const size_t				nbloops;
		size_t					iterationnumber;

		std::vector< std::vector< integer > > 	residues;
		integer _product;
		integer _midprod;

	public:
		GivaroRnsFixedCRA(const std::vector<integer>& primes)
				: Father_t(primes),
				  nbloops(primes.size()),
				  iterationnumber(0),
				  _product(1)
		{
			for(size_t i=0; i<primes.size(); ++i)
				_product *= primes[i];
			::Givaro::Integer::div(_midprod,_product,2);
		}

		Integer& getModulus(Integer& m)
		{
			return m=_product;
		}

		template<template<class> class Vect>
		Vect<Integer>& getResidue(Vect<Integer>& r)
		{
			return result(r);
		}

		template< template<class, class> class Vect, template <class> class Alloc>
		void initialize (const Domain& D, const Vect<DomainElement, Alloc<DomainElement> >& e)
		{
			residues.resize(e.size());
			this->progress(D,e);
		}

		template< template<class, class> class Vect, template <class> class Alloc>
		void progress (const Domain& D, const Vect<DomainElement, Alloc<DomainElement> >& e)
		{
			++iterationnumber;
			typename Vect<DomainElement, Alloc<DomainElement> >::const_iterator eit=e.begin();
			std::vector<std::vector< Integer > >::iterator rit = residues.begin();

			for( ; eit != e.end(); ++eit, ++rit) {
				Integer tmp;
				D.convert(tmp, *eit);
				rit->push_back(tmp);
			}

		}

		template<template<class, class> class Vect, template <class> class Alloc>
		Vect<Integer, Alloc<Integer> >& result (Vect<Integer, Alloc<Integer> > &d)
		{
			d.resize(0);
			for(std::vector<std::vector< Integer > >::const_iterator rit = residues.begin(); rit != residues.end(); ++rit) {
				Integer tmp;
				RnsToRing(tmp, *rit);
				linbox_check(tmp>=0);
				linbox_check(tmp<_product);
				if (tmp>_midprod)
					tmp -= _product ;
				linbox_check(tmp<=_midprod);
				d.push_back(tmp);
			}
			return d;
		}

		bool terminated()
		{
			return iterationnumber >= nbloops;
		}

		bool noncoprime(const Integer& i) const
		{
			return false;
		}




	};

}


#endif //__LINBOX_cra_givrnsfix_H
