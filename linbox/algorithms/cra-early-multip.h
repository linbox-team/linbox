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
	template<class SingleCRA_Type, class RandIter_Type>
	struct EarlyMultipCRARGen : public SingleCRA_Type, public FullMultipCRA {
        // Could be much faster
        // - do not compute twice the product of moduli
        // - reconstruct one element of e until Early Termination,
        //   then only, try a random linear combination.

		using Self_t = EarlyMultipCRARGen<SingleCRA_Type,RandIter_Type>;
        using SingleParent = SingleCRA_Type;
        using MultiParent = FullMultipCRA;

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

        template <typename... Args>
		EarlyMultipCRARGen(RandIter_Type * rgen, Args&&... args) :
			SingleParent(std::forward<Args>(args)...),
            MultiParent(),
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

        /** @brief Update the early termination vector and values.
         * @return true iff the method is now terminated.
         */
		bool changeVector()
		{
            /* Pick a new linear combination */
            choose_linear_comb(randv_.size());

            /* re-compute residues for SingleParent */
            bool inited = false;

            for (auto shelf_it = MultiParent::shelves_begin();
                 shelf_it != MultiParent::shelves_end();
                 ++shelf_it)
            {
                auto newval = dot(shelf_it->mod(), shelf_it->residue, randv_);
                if (shelf_it->occupied) {
                    if (! inited) {
                        SingleParent::initialize(shelf_it->mod(), newval);
                        inited = true;
                    }
                    else {
                        SingleParent::progress(shelf_it->mod(), newval,
                                               shelf_it->count);
                    }
                    if (SingleParent::terminated()) return true;
                }
            }

            return false;
		}

	protected:

        template <typename Domain, typename Elt2>
        static inline void axpyin (typename Domain::Element& out, const typename Domain::Element& a, const Elt2& b, const Domain& D) {
            typename Domain::Element tmp;
            D.axpyin(out, a, D.init(tmp, b));
        }

        template <typename Elt2>
        static inline void axpyin (Integer& out, const Integer& a, const Elt2& b, const Integer& D) {
            out = (out + a*b) % D;
        }

        template <typename DomainOrInteger, typename Vect1, typename Vect2>
        static typename Vect1::value_type
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


	/** @brief vector CRA with early termination which creates its own PRNG object
	 * @ingroup CRA
	 */
	template<class RandIter_Type=Givaro::GivRandom>
	struct EarlyMultipCRA : public EarlyMultipCRARGen<EarlySingleCRA,RandIter_Type> {
        using Parent_t = EarlyMultipCRARGen<EarlySingleCRA,RandIter_Type>;

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


	/** @brief vector CRA with guaranteed error probability which creates its own PRNG object
	 * @ingroup CRA
	 */
	template<class RandIter_Type=Givaro::GivRandom>
	struct ProbMultipCRA : public EarlyMultipCRARGen<ProbSingleCRA,RandIter_Type> {
        using Parent_t = EarlyMultipCRARGen<ProbSingleCRA,RandIter_Type>;

	protected:
        // PRNG
        RandIter_Type actual_rgen_;

	public:

        /** @brief Creates a new probabiistic CRA object.
         *
         * @param bitbound  Upper bound on the number of bits in the result
         * @param failprob  Upper bound on the probability of failure
         * @param rgargs  Arguments to PRNG constructor
         */
        template <typename... RGArgs>
		ProbMultipCRA(size_t bitbound, double failprob,
                      RGArgs&&... rgargs) :
            Parent_t(&actual_rgen_, bitbound, failprob),
            actual_rgen_(std::forward<RGArgs>(rgargs)...)
		{}

        /** @brief Creates a new probabiistic CRA object.
         *
         * @param bitbound  Upper bound on the number of bits in the result
         */
        ProbMultipCRA(size_t bitbound) :
            Parent_t(&actual_rgen_, bitbound)
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
