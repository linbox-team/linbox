/* Copyright (C)  LinBox
 * author: B. David Saunders and Zhendong Wan
 * parallelized for BOINC computing by Bryan Youse
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

/*! @file algorithms/cra-kaapi.h
 * @ingroup algorithms
 * @ingroup CRA
 * @brief NO DOC
 */

#ifndef __LINBOX_cra_kaapi_H
#define __LINBOX_cra_kaapi_H

#include <vector>
#include <cstdlib>
#include <utility>

#include "linbox/util/timer.h"
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"

#include "linbox/kaapi/communicate.h"

namespace LinBox
{

	template<class Function, class Domain >
	struct Residue {

		Function *_f;
		Residue() :
			_f(0)
		{}
		Residue(Function& ff) :
			_f(&ff)
		{}
		Residue(const Residue<Function,Domain>&r ) :
			_f(r._f)
		{}

		typename Domain::Element operator()(Domain D)
		{
			typename Domain::Element d;
			D.init(d);
			return (*_f)(d,D);
		}
	};

	/*! No doc.
	 * @ingroup CRA
	 */
	template<class CRABase>
	struct ChineseRemainder {

		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
	protected:
		CRABase Builder_;

	public:
		template<class Param>
		ChineseRemainder(const Param& b) :
			Builder_(b)
		{}

		/**
		 * CRA loop subroutine.
		 * @ingroup CRA
		 * given a function and a prime, this returns the residue by applying given function
		 * this must be thread safe and communicable
		 * \param primeiter Prime Integer
		 * \param Iteration function used to compute residue
		 * \param[out] res the residue
		 * \return the residue
		 */
		template<class Int, class Function, class PrimeIterator>
		Int& operator() (Int& res, Function& Iteration, PrimeIterator& primeiter)
		{

			double start = Util::WallTimer::gettime();

			++primeiter;
			Domain D(*primeiter);
			DomainElement r; D.init(r);
			Builder_.initialize( D, Iteration(r, D) );

			// the task used to extract residue
			Residue<Function,Domain> residue(Iteration);

			size_t nb_primes = 4;
			size_t nb_done=0;

			while( ! Builder_.terminated() )
			{
				Domain domains [nb_primes];
				DomainElement domainelements[nb_primes];

				// generate the array of domain
				for (size_t i=0;i < nb_primes;i++) {
					do {
						++primeiter;
					} while( Builder_.noncoprime(*primeiter) ) ;
					domains[i]=Domain(*primeiter);
				}

				// recursively call the send function
				a1::transform(
					      domains,
					      domains+nb_primes,
					      domainelements,
					      residue
					     );

				// when it's done, analyze the result
				for(size_t i=0;i<nb_primes;++i) {
					Builder_.progress(domains[i], domainelements[i]);
				}

				nb_done+=nb_primes;
				nb_primes=nb_done/2;
			}

			std::cout << "TIME=" << Util::WallTimer::gettime()-start << std::endl;
			return Builder_.result(res);
		}

		template<class Int, template <class T> class Vect, class Function, class PrimeIterator>
		Vect<Int> & operator() (Vect<Int>& res, Function& Iteration, PrimeIterator& primeiter)
		{

			++primeiter;
			Domain D(*primeiter);
			Vect<DomainElement> r;
			Builder_.initialize( D, Iteration(r, D) );

			while( ! Builder_.terminated() ) {
				++primeiter; while(Builder_.noncoprime(*primeiter) ) ++primeiter;
				Domain D(*primeiter);
				Vect<DomainElement> r;
				Builder_.progress( D, Iteration(r, D) );
			}
			return Builder_.result(res);
		}

	};

}

/*
 * marshalling operator,
 * WARNING: those are dummy ones, the real ones are *really* hard to implement because of BlackBoxes
 * whose interface is hidden and does not make it easy to be communicable
 */

template<class Function, class Domain >
a1::OStream& operator<<( a1::OStream& out, const LinBox::Residue<Function, Domain>&  )
{
	return out ;
}

template<class Function, class Domain >
a1::IStream& operator>>( a1::IStream& in,  LinBox::Residue<Function, Domain>&  )
{
	return in;
}
#endif //__LINBOX_cra_kaapi_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

