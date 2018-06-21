/* Copyright (C) 2007 LinBox
 * Written by JG Dumas
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

/*! @file algorithms/cra-early-multip.h
 * @ingroup algorithms
 * @brief vector Chinese remaindering with early termination
 */

#ifndef __LINBOX_cra_early_multip_H
#define __LINBOX_cra_early_multip_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

#include "linbox/algorithms/cra-single.h"
#include "linbox/algorithms/cra-full-multip.h"


namespace LinBox
{

	/*!  @brief NO DOC
	 * @ingroup CRA
	 * @brief vector CRA with early termination and a PRNG reference
	 */
	template<class Domain_Type, class RandIter_Type>
	struct EarlyMultipCRARGen : public EarlySingleCRA<Domain_Type>, public FullMultipCRA<Domain_Type> {
        // Could be much faster
        // - do not compute twice the product of moduli
        // - reconstruct one element of e until Early Termination,
        //   then only, try a random linear combination.

		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef EarlyMultipCRARGen<Domain,RandIter_Type> Self_t;
        using SingleParent = EarlySingleCRA<Domain>;
        using MultiParent = FullMultipCRA<Domain>;

        const unsigned short LINEAR_COMB_RBITS = 15; // # of random bits in linear combination

	protected:
        // PRNG
        RandIter_Type * const rgen_;

		// Random coefficients for a linear combination
		// of the elements to be reconstructed
		std::vector< unsigned long >      	randv_;

        /** @brief Chooses a new random linear combination
         * @ param size  is the sie of the vector
         * XXX: Does not change any stored moduli. Generally should be called
         * only in initialization.
         */
        inline void choose_linear_comb (unsigned long size) {
			// Random coefficients for a linear combination
			// of the elements to be reconstructed
			randv_. resize ( size );
			for ( auto int_p = randv_. begin(); int_p != randv_. end(); ++ int_p)
                *int_p = (*rgen_)() >> LINEAR_COMB_RBITS;
        }

		Integer& result(Integer &d) { std::cout << "should not be called" << std::endl; return d ;} ; // DON'T TOUCH
	public:

		EarlyMultipCRARGen( RandIter_Type * rgen,
                            const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD) :
			SingleParent(EARLY), MultiParent(),
            rgen_ (rgen)
		{}

		Integer& getModulus(Integer& m)
		{
			SingleParent::getModulus(m);
			return m;
		}
		Integer& getResidue(Integer& m)
		{
			SingleParent::getResidue(m);
			return m;
		}

		template<template<class T> class Vect>
		Vect<Integer>& getResidue(Vect<Integer>& m)
		{
			MultiParent::getResidue(m);
			return m;
		}

		//! Init
        template <class DomainOrInteger, class Vect>
        void initialize (const DomainOrInteger& D, const Vect& e) {
            choose_linear_comb(e.size());
			SingleParent::initialize(D, dot(D, e, randv_));
			MultiParent::initialize(D, e);
		}

		//! Progress
        template <class DomainOrInteger, class Vect>
        void progress (const DomainOrInteger& D, const Vect& e)
		{

			SingleParent::progress(D, dot(D, e, randv_));
			MultiParent::progress(D, e);
		}

		//! Result
		template<class Vect>
                Vect& result(Vect& d)
		{
			return MultiParent::result(d);
		}

		//! terminate
		bool terminated()
		{
			return SingleParent::terminated();
		}

		bool noncoprime(const Integer& i) const
		{
			return SingleParent::noncoprime(i);
		}

		bool changeVector()
		{
            choose_linear_comb(randv_.size());
			std::vector<Integer> e(randv_.size());
            // TODO clean this up, decouple with early and full classes
			/* clear CRAEarlySingle; */
			SingleParent::occurency_ = 0;
			SingleParent::nextM_ = 1UL;
			SingleParent::primeProd_ = 1UL;
			SingleParent::residue_ = 0;

			/* Computation of residue_ */
			std::vector< LazyProduct >::iterator _mod_it = MultiParent::RadixPrimeProd_.begin();// list of prime products
			std::vector< std::vector<Integer> >::iterator _tab_it = MultiParent::RadixResidues_.begin();// list of residues as vectors of size 1
			std::vector< bool >::iterator    _occ_it = MultiParent::RadixOccupancy_.begin();//flags of occupied fields
			int prev_shelf=0, shelf = 0;
			for (;_occ_it != MultiParent::RadixOccupancy_.end(); ++_mod_it, ++_tab_it, ++_occ_it ) {
				++shelf;
				if (*_occ_it) {
					Integer D = _mod_it->operator()();
					std::vector<Integer> e_v(randv_.size());
					e_v = *_tab_it;
                    auto z = dot(D, e_v, randv_);
					Integer prev_residue_ = SingleParent::residue_;
					SingleParent::progress(D,z);
					if (prev_residue_ == SingleParent::residue_ )
						SingleParent::occurency_ = SingleParent::occurency_ +  (shelf - prev_shelf);
					if ( SingleParent::terminated() ) {
						return true;
					}
					prev_shelf = shelf;
				}
			}
			return false;
		}

	protected:

        template <typename DomElem, typename Elt2>
        inline void axpyin (DomElem& out, const DomElem& a, const Elt2& b, const Domain& D) {
            DomElem tmp;
            D.axpyin(out, a, D.init(tmp, b));
        }

        template <typename Elt2>
        inline void axpyin (Integer& out, const Integer& a, const Elt2& b, const Integer& D) {
            out = (out + a*b) % D;
        }

        template <typename DomainOrInteger, typename Vect1, typename Vect2>
        typename Vect1::value_type
        dot (const DomainOrInteger& D, const Vect1& v1, const Vect2& v2) {
            typename Vect1::value_type z = 0;
            auto v1_p = v1.begin();
            auto v2_p = v2.begin();
            for (; v1_p != v1.end(); ++v1_p, ++v2_p) {
                axpyin(z, *v1_p, *v2_p, D);
            }
            return z;
        }

	};


	/*!  @brief NO DOC
	 * @ingroup CRA
	 * @brief vector CRA with early termination which creates its own PRNG object
	 */
	template<class Domain_Type, class RandIter_Type=Givaro::GivRandom>
	struct EarlyMultipCRA : public EarlyMultipCRARGen<Domain_Type,RandIter_Type> {
        using Parent_t = EarlyMultipCRARGen<Domain_Type,RandIter_Type>;
        using Domain = typename Parent_t::Domain;
        using DomainElement = typename Parent_t::DomainElement;

	protected:
        // PRNG
        RandIter_Type actual_rgen_;

	public:

        template <typename... RGArgs>
		EarlyMultipCRA(const unsigned long EARLY,
                       RGArgs&&... rgargs) :
            Parent_t(&actual_rgen_, EARLY),
            actual_rgen_(std::forward<RGArgs>(rgargs)...)
		{}

        EarlyMultipCRA() :
            Parent_t(&actual_rgen_)
        {}
    };
}
#endif //__LINBOX_cra_early_multip_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
