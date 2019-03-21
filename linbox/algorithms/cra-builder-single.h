/* Copyright (C) 2007 LinBox
 * written by JG Dumas
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

/*! @file algorithms/cra-builder-single.h
 * @ingroup algorithms
 * @brief Chinese remaindering of a single value
 */


#ifndef __LINBOX_cra_single_H
#define __LINBOX_cra_single_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-builder-full-multip.h"
#include <vector>
#include <array>
#include <utility>

namespace LinBox
{
	/** @brief Lower bound on number of b-bit primes
	 * @ingroup CRA
	 */
	inline uint64_t primes_count(size_t pbits)
	{
		static std::array<uint64_t, 36> pctable = {{
			0, 0, 2, 2, 2, 5, 7, 13, 23, 43, 75, 137, 255,
			464, 872, 1612, 3030, 5709, 10749, 20390, 38635,
			73586, 140336, 268216, 513708, 985818, 1894120,
			3645744, 7027290, 13561907, 26207278, 50697537,
			98182656, 190335585, 369323305, 717267168,
		}}; // source: http://oeis.org/A162145/list
		static double C = 3. / (10. * log(2.));
		if (pbits < pctable.size()) {
			return pctable[pbits];
		}
		else if (pbits <= 71) {
			// source: Rosser & Schoenfeld
			return ceil((C / (pbits - 1)) * pow(2., pbits));
		}
		else {
			return std::numeric_limits<uint64_t>::max();
		}
	}

	/**  @brief Abstract base class for CRA builders
	 *
	 *   Subclasses must implement the termination functionality.
	 *
	 * @ingroup CRA
	 *
	 */
	template<class Domain_Type>
	struct CRABuilderSingleBase {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef CRABuilderSingleBase<Domain> Self_t;

	protected:
		// PrimeProd*nextM_ is the modulus
		Integer 	primeProd_;
		Integer		nextM_;
		Integer 	residue_; 	// remainder to be reconstructed

#ifdef _LB_CRATIMING
		mutable Timer tInit, tIteration, tImaging, tIRecon, tOther;
		mutable CRATimer totalTime;
        mutable size_t IterCounter_;
#endif

		/** @brief Update the residue and check whether it changed
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 *
		 * The initialize function should be called at least once before
		 * calling this one.
		 *
		 * @param D	The modulus of the new image
		 * @param e	The residue modulo D
		 * @returns	true iff the residue changed with this update
		 */
		bool progress_check (const Integer& D, const Integer& e)
		{
			// Precondition : initialize has been called once before
#ifdef _LB_CRATIMING
			tIRecon.clear();
			tIRecon.start();
            ++IterCounter_;
#endif
			bool was_updated = false;
			primeProd_ *= nextM_;
			nextM_ =D;
			Integer u0 = residue_ % D;//0
			Integer u1 = e % D;//e
			Integer m0 = primeProd_;//1
			if (u0 != u1) {
				was_updated = true;
				inv(m0, m0, D); // res <-- m0^{-1} mod m1//1
				u0 = u1 - u0;           // tmp <-- (u1-u0)//e
				u0 *= m0;       // res <-- (u1-u0)( m0^{-1} mod m1 )//e
				u0 %= D;
				Integer tmp(u0);//e
				if (u0 < 0)
					tmp += D;//e+D
				else
					tmp -= D;//e-D
				if (absCompare(u0,tmp) > 0) u0 = tmp;
				u0 *= primeProd_;          // res <-- (u1-u0)( m0^{-1} mod m1 ) m0       and res <= (m0m1-m0)
				residue_ += u0;   // res <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0  and res <  m0m1
			}
#ifdef _LB_CRATIMING
			tIRecon.stop();
			totalTime.ttIRecon += tIRecon;
#endif
			return was_updated;
		}

		/** @brief Update the residue and check whether it changed
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 *
		 * The initialize function should be called at least once before
		 * calling this one.
		 *
		 * @param D	The modulus of the new image
		 * @param e	The residue modulo D
		 * @returns	true iff the residue changed with this update
		 */
		bool progress_check (const Domain& D, const DomainElement& e)
		{
			// Precondition : initialize has been called once before
#ifdef _LB_CRATIMING
			tIRecon.clear();
			tIRecon.start();
            ++IterCounter_;
#endif
			bool was_updated = false;
			primeProd_ *= nextM_;
			D.characteristic( nextM_ );

			DomainElement u0;
			if (! D.areEqual( D.init(u0, residue_), e)) {
				was_updated = true;

				D.negin(u0);       	// u0  <-- -u0
				D.addin(u0, e);    	// u0  <-- e-u0

				DomainElement m0;
				D.init(m0, primeProd_);
				D.invin(m0);       	// m0  <-- m0^{-1} mod nextM_
				D.mulin(u0, m0);   	// u0  <-- (e-u0)( m0^{-1} mod nextM_ )

				Integer res;
				D.convert(res, u0);	// res <-- (e-u0)( m0^{-1} mod nextM_ )
				// and res < nextM_

				Integer tmp(res);
				tmp -= nextM_;
				if (absCompare(res,tmp)>0) res = tmp; // Normalize

				res *= primeProd_;	// res <-- (e-u0)( m0^{-1} mod nextM_ ) m0
				// and res <= (m0.nextM_-m0)

				residue_ += res;	// <-- u0 + (e-u0)( m0^{-1} mod nextM_ ) m0
				// and res <  m0.nextM_
			}
#ifdef _LB_CRATIMING
			tIRecon.stop();
			totalTime.ttIRecon += tIRecon;
#endif
			return was_updated;
		}

	public:

		CRABuilderSingleBase() :
			primeProd_(1U),
			nextM_(1U)
#ifdef _LB_CRATIMING
            , IterCounter_(0)
#endif
		{
#ifdef _LB_CRATIMING
			clearTimers();
			totalTime.clear();
#endif
		}

		/** @brief Initialize the CRA with the first residue.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 * This function must be called just once. Subsequent calls
		 * should be made to the progress() function.
		 *
		 * @param D	The modulus
		 * @param e	The residue
		 */
		void initialize (const Integer& D, const Integer& e)
		{
#ifdef _LB_CRATIMING
			tInit.clear();
			tInit.start();
            IterCounter_=1;
#endif
			primeProd_ = D;
			nextM_ = 1U;
			residue_ = e;
#ifdef _LB_CRATIMING
			tInit.stop();
			totalTime.ttInit += tInit;
#endif
		}

		/** @brief Initialize the CRA with the first residue.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 * This function must be called just once. Subsequent calls
		 * should be made to the progress() function.
		 *
		 * @param D	The modulus
		 * @param e	The residue
		 */
		void initialize (const Domain& D, const DomainElement& e)
		{
#ifdef _LB_CRATIMING
			tInit.clear();
			tInit.start();
            IterCounter_=1;
#endif
			D.characteristic( primeProd_ );
			nextM_ = 1U;
			D.convert( residue_, e);
#ifdef _LB_CRATIMING
			tInit.stop();
			totalTime.ttInit += tInit;
#endif
		}

		/** @brief Gets the result recovered so far.
		 *
		 * (This is the same as getResidue.)
		 */
		Integer& result(Integer& d)
		{
			return d=residue_;
		}

		/** @brief Gets the result recovered so far.
		 */
		Integer& getResidue(Integer& r )
		{
			return r= residue_;
		}

		/** @brief Gets the modulus of the result recovered so far.
		 */
		Integer& getModulus(Integer& m)
		{

#ifdef _LB_CRATIMING
			tOther.clear();
			tOther.start();
#endif
			m = primeProd_ * nextM_;
#ifdef _LB_CRATIMING
			tOther.stop();
			totalTime.ttOther += tOther;
#endif
			return m;
		}

		/** @brief Checks whether i is co-prime to the modulus so far.
		 *
		 * @return	true iff i shares a common factor with the modulus
		 */
		bool noncoprime(const Integer& i) const
		{
			Integer g;
			return ( (gcd(g, i, nextM_) != 1) || (gcd(g, i, primeProd_) != 1) );
		}

		/** @brief Returns a lower bound on the number of bits in the modulus.
		 */
		decltype(Integer().bitsize()) modbits() const
		{
			return primeProd_.bitsize() + nextM_.bitsize() - 1;
		}

		virtual ~CRABuilderSingleBase() {}

#ifdef _LB_CRATIMING
		void clearTimers() const
		{
			tInit.clear();
			//tIteration.clear();
			//tImaging.clear();
			tIRecon.clear();
			tOther.clear();
		}
	public:
		inline std::ostream& printTime(const Timer& timer, const char* title, std::ostream& os, const char* pref = "") const
		{
			if (timer.count() > 0) {
				os << pref << title;
				for (int i=strlen(title)+strlen(pref); i<28; i++)
					os << ' ';
				return os << timer << std::endl;
			}
			else
				return os;
		}

		inline std::ostream& printCRATime(const CRATimer& timer, const char* title, std::ostream& os) const
		{
			printTime(timer.ttInit, " Init: ", os, title);
			//printTime(timer.ttImaging, "Imaging", os, title);
			//printTime(timer.ttIteration, "Iteration", os, title);
			printTime(timer.ttIRecon, " Integer reconstruction: ", os, title);
			printTime(timer.ttOther, " Other: ", os, title);
			return os;
		}

		std::ostream& reportTimes(std::ostream& os) const
		{
            os <<  "CRA Iterations:" << IterCounter_ << "\n";
			printCRATime(totalTime, "CRA Time", os);
			return os;
		}
#endif

	};


	/**  @brief Heuristic Chinese Remaindering with early termination
	 *
	 * This approach stops as soon as the modulus doesn't changed for some
	 * fixed number of steps in a row.
	 *
	 * @ingroup CRA
	 *
	 */
	template<class Domain_Type>
	struct CRABuilderEarlySingle :public CRABuilderSingleBase<Domain_Type> {
		typedef CRABuilderSingleBase<Domain_Type>	Base;
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef CRABuilderEarlySingle<Domain> Self_t;

		const unsigned int    EARLY_TERM_THRESHOLD;

	protected:
		unsigned int    occurency_;	// number of equalities

	public:

		/** @brief Creates a new heuristic CRA object.
		 *
		 * @param	EARLY how many unchanging iterations until termination.
		 */
		CRABuilderEarlySingle(const size_t EARLY=LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD) :
			EARLY_TERM_THRESHOLD((unsigned)EARLY-1),
			occurency_(0U)
		{ }

		/** @brief Initialize the CRA with the first residue.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 * This function must be called just once. Subsequent calls
		 * should be made to the progress() function.
		 *
		 * Either the types of D and e should both be Integer,
		 * or D is the domain type (e.g., Modular<double>) and
		 * e is the element type (e.g., double)
		 *
		 * @param D	The modulus
		 * @param e	The residue
		 */
		template <typename ModType, typename ResType>
		void initialize (const ModType& D, const ResType& e)
		{
			Base::initialize(D,e);
			occurency_ = 1;
		}

		/** @brief Update the residue and termination condition.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 *
		 * The initialize function must be called at least once before
		 * calling this one.
		 *
		 * Either the types of D and e should both be Integer,
		 * or D is the domain type (e.g., Modular<double>) and
		 * e is the element type (e.g., double)
		 *
		 * @param D	The modulus of the new image
		 * @param e	The residue modulo D
		 */
		template <typename ModType, typename ResType>
		void progress (const ModType& D, const ResType& e)
		{
			// Precondition : initialize has been called once before
			// linbox_check(occurency_ > 0);
			if (Base::progress_check(D,e))
				occurency_ = 1;
			else
				occurency_ ++;
		}

		/** @brief Checks whether the CRA is (heuristically) finished.
		 *
		 * @return true iff the early termination condition has been reached.
		 */
		bool terminated() const
		{
			return occurency_ > EARLY_TERM_THRESHOLD;
		}
	};


	/**  @brief Chinese Remaindering with guaranteed probability bound and early termination.
	 *
	 * This strategy terminates when the probability of failure is below a
	 * certain threshold.
	 *
	 * @ingroup CRA
	 *
	 */
	template<class Domain_Type>
	struct CRABuilderProbSingle :public CRABuilderSingleBase<Domain_Type> {
		typedef CRABuilderSingleBase<Domain_Type>	Base;
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef CRABuilderProbSingle<Domain> Self_t;

	protected:
		const size_t bitbound_;
		const double failbound_;
		double curfailprob_; // the probability the result right now is incorrect

		size_t mod_bitsize(const Integer& D) const { return D.bitsize(); }

		template <typename Field>
		size_t mod_bitsize(const Field& D) const {
			Integer p;
			D.characteristic(p);
			return p.bitsize();
		}

	public:
		/** @brief Creates a new probabilistic CRA object.
		 *
		 * @param	bitbound  An upper bound on the number of bits in the result.
		 * @param	failprob  An upper bound on the probability of failure.
		 */
		CRABuilderProbSingle(const size_t bitbound, const double failprob = LINBOX_DEFAULT_FAILURE_PROBABILITY) :
			bitbound_(bitbound),
			failbound_(failprob),
			curfailprob_(-1.)
		{
		}

		/** @brief Initialize the CRA with the first residue.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 * This function must be called just once. Subsequent calls
		 * should be made to the progress() function.
		 *
		 * Either the types of D and e should both be Integer,
		 * or D is the domain type (e.g., Modular<double>) and
		 * e is the element type (e.g., double)
		 *
		 * @param D	The modulus
		 * @param e	The residue
		 */
		template <typename ModType, typename ResType>
		void initialize (const ModType& D, const ResType& e)
		{
			Base::initialize(D,e);
			curfailprob_ = 1.;
		}

		/** @brief Update the residue and termination condition.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 *
		 * The initialize function must be called at least once before
		 * calling this one.
		 *
		 * Either the types of D and e should both be Integer,
		 * or D is the domain type (e.g., Modular<double>) and
		 * e is the element type (e.g., double)
		 *
		 * @param D	The modulus of the new image
		 * @param e	The residue modulo D
		 */
		template <typename ModType, typename ResType>
		void progress (const ModType& D, const ResType& e)
		{
			// Precondition : initialize has been called once before
			// linbox_check(curfailprob_ >= 0.);
			if (Base::progress_check(D,e)) {
				curfailprob_ = (Base::modbits() > bitbound_) ? 0. : 1.;
			}
			else {
				// failure iff (failed so far) AND (this image fooled you)
				// = curfailprob * (# bad primes in range) / (total # primes in range)
				// <= curfailprob * (floor((bitbound_ - modbits + 1)/(pbits - 1)) / (# pbits primes)
				auto mbits = Base::modbits();
				if (mbits > bitbound_)
					curfailprob_ = 0.;
				else {
					auto pbits = mod_bitsize(D);
					double failupdate = double((bitbound_ + 1 - Base::modbits()) / (pbits - 1))
							    / primes_count(pbits);
					if (failupdate < 1.)
						curfailprob_ *= failupdate;
				}
			}
		}

		/** @brief Checks whether the CRA is (heuristically) finished.
		 *
		 * @return true iff the early termination condition has been reached.
		 */
		bool terminated() const
		{
			return curfailprob_ <= failbound_;
		}
	};


	/**  @brief Chinese Remaindering with full precision and no chance of failure.
	 *
	 * @ingroup CRA
	 *
	 */
	template<class Domain_Type>
	struct CRABuilderFullSingle :public CRABuilderFullMultip<Domain_Type>{
		typedef CRABuilderSingleBase<Domain_Type>	Base;
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef CRABuilderFullSingle<Domain> Self_t;
        using MultiParent = CRABuilderFullMultip<Domain>;

	public:
		/** @brief Creates a new deterministic CRA object.
		 *
		 * @param	bitbound  An upper bound on the number of bits in the result.
		 * @param	failprob  An upper bound on the probability of failure.
		 */
		CRABuilderFullSingle(const size_t bitbound) :
            MultiParent((double)bitbound)
		{
		}

		/** @brief Initialize the CRA with the first residue.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 * This function must be called just once. Subsequent calls
		 * should be made to the progress() function.
		 *
		 * Either the types of D and e should both be Integer,
		 * or D is the domain type (e.g., Modular<double>) and
		 * e is the element type (e.g., double)
		 *
		 * @param D	The modulus
		 * @param e	The residue
		 */
		template <typename ModType, typename ResType>
		void initialize (const ModType& D, const ResType& e)
		{
            std::array<ResType,1> v {e};
            MultiParent::initialize(D,v);
		}

		/** @brief Update the residue and termination condition.
		 *
		 * The eventually-recovered number will be congruent to e modulo D.
		 *
		 * The initialize function must be called at least once before
		 * calling this one.
		 *
		 * Either the types of D and e should both be Integer,
		 * or D is the domain type (e.g., Modular<double>) and
		 * e is the element type (e.g., double)
		 *
		 * @param D	The modulus of the new image
		 * @param e	The residue modulo D
		 */
		template <typename ModType, typename ResType>
		void progress (const ModType& D, const ResType& e)
		{
			// Precondition : initialize has been called once before
            std::array<ResType,1> v {e};
            MultiParent::progress(D, v);
		}

		/** @brief Gets the result recovered so far.
		 */
		inline const Integer& getResidue() const
		{
            return MultiParent::getResidue().front();
		}

		/** @brief Gets the result recovered so far.
		 */
		inline Integer& getResidue(Integer& r) const
		{
			return r = getResidue();
		}

		/** @brief alias for getResidue
		 */
		inline Integer& result(Integer& r) const
		{
			return r = getResidue();
		}
	};

}

#endif //__LINBOX_cra_single_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
