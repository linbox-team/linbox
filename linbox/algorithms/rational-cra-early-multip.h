/* Copyright (C) 2010 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __LINBOX_rational_early_multip_cra_H
#define __LINBOX_rational_early_multip_cra_H

#include "givaro/zring.h"
#include "linbox/algorithms/rational-cra-early-single.h"
#include "linbox/algorithms/rational-cra-full-multip.h"

namespace LinBox
{

	struct EarlyMultipRatCRA : public EarlySingleRatCRA, public FullMultipRatCRA {
        using SingleParent = EarlySingleRatCRA;
        using MultiParent = FullMultipRatCRA;
		typedef EarlyMultipRatCRA	Self_t;
	protected:
		// Random coefficients for a linear combination
		// of the elements to be reconstructed
		std::vector< unsigned long >      	randv;

		void initialize (const Integer& D, const Integer& e) {return;}; // DON'T TOUCH
		void progress (const Integer & D, const Integer & e) {return;};
        template <class Domain>
		void initialize (const Domain& D, const typename Domain::Element& e){return;};
		Integer& result(Integer& d){return d;};
        template <class Domain>
		void progress (const Domain& D, const typename Domain::Element& e){return;};






	public:


		EarlyMultipRatCRA(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD) :
			SingleParent(EARLY), MultiParent()
		{ }

		//!init
		template<class Domain, template<class, class> class Vect, template <class> class Alloc>
		void initialize (const Domain& D, const Vect <typename Domain::Element, Alloc<typename Domain::Element> >& e)
		{
			// Random coefficients for a linear combination
			// of the elements to be reconstructed
			srand48(BaseTimer::seed());
			randv. resize ( e.size() );
			for ( std::vector<unsigned long>::iterator int_p = randv. begin();
			      int_p != randv. end(); ++ int_p)
				*int_p = ((unsigned long)lrand48()) % 20000;

			typename Domain::Element z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			SingleParent::initialize(D,dot(z, D, e, randv) );
			MultiParent::initialize(D, e);
		}

        template <class Domain>
		void initialize (const Domain& D, const BlasVector<Domain>& e)
		{
			// Random coefficients for a linear combination
			// of the elements to be reconstructed
			srand48(BaseTimer::seed());
			randv. resize ( e.size() );
			for ( std::vector<unsigned long>::iterator int_p = randv. begin();
			      int_p != randv. end(); ++ int_p)
				*int_p = ((unsigned long)lrand48()) % 20000;

			typename Domain::Element z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			SingleParent::initialize(D,dot(z, D, e, randv) );
			MultiParent::initialize(D, e);
		}

		//!progress
		template<class Domain, template<class,class> class Vect, template <class> class Alloc>
		void progress (const Domain& D, const Vect<typename Domain::Element, Alloc<typename Domain::Element> >& e)
		{
			typename Domain::Element z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			SingleParent::progress(D, dot(z, D, e, randv));
			MultiParent::progress(D, e);
		}

        template <class Domain>
		void progress (const Domain& D, const BlasVector<Domain>& e)
		{
			typename Domain::Element z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			SingleParent::progress(D, dot(z, D, e, randv));
			MultiParent::progress(D, e);
		}

		//!result
		template<template<class, class> class Vect, template <class> class Alloc>
		Vect<Integer, Alloc<Integer> >& result(Vect<Integer, Alloc<Integer> >& num, Integer& den)
		{
			return MultiParent::result(num, den);
		}

		BlasVector<Givaro::ZRing<Integer> >& result(BlasVector<Givaro::ZRing<Integer>>& num, Givaro::ZRing<Integer>::Element& den)
		{
			return MultiParent::result(num, den);
		}

		//!tools
		bool terminated()
		{
			return SingleParent::terminated();
		}

		bool noncoprime(const Integer& i) const
		{
			return SingleParent::noncoprime(i);
		}

	protected:

		template <class Domain, template<class, class> class Vect1, template <class> class Alloc, class Vect2>
		typename Domain::Element& dot (typename Domain::Element& z, const Domain& D,
				    const Vect1<typename Domain::Element, Alloc<typename Domain::Element> >& v1, const Vect2& v2)
		{

			D.assign(z,D.zero); typename Domain::Element tmp;
			typename Vect1<typename Domain::Element, Alloc<typename Domain::Element> >::const_iterator v1_p;
			typename Vect2::const_iterator v2_p;
			for (v1_p  = v1. begin(), v2_p = v2. begin();
			     v1_p != v1. end();
			     ++ v1_p, ++ v2_p)
				D.axpyin(z, (*v1_p), D.init(tmp, (*v2_p)));
#if 0
			commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "v: " << v2 << std::endl;
			commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "z: " << z << std::endl;
#endif
			return z;
		}

		template <class Domain, class Vect2>
		typename Domain::Element& dot (typename Domain::Element& z, const Domain& D,
				    const BlasVector<Domain >& v1, const Vect2& v2)
		{

			D.assign(z,D.zero); typename Domain::Element tmp;
			typename BlasVector<Domain >::const_iterator v1_p;
			typename Vect2::const_iterator v2_p;
			for (v1_p  = v1. begin(), v2_p = v2. begin();
			     v1_p != v1. end();
			     ++ v1_p, ++ v2_p)
				D.axpyin(z, (*v1_p), D.init(tmp, (*v2_p)));
#if 0
			commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "v: " << v2 << std::endl;
			commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "z: " << z << std::endl;
#endif
			return z;
		}

	};
}

#endif //__LINBOX_rational_early_multip_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
