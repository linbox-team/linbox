/* linbox/algorithms/cra-domain-sequential.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Time-stamp: <28 Oct 22 18:29:45 Jean-Guillaume.Dumas@imag.fr>
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

/*! @file algorithms/cra-domain-sequential.h
 * @brief Sequencial version of \ref CRA
 * @ingroup CRA
 */

#ifndef __LINBOX_sequential_cra_H
#define __LINBOX_sequential_cra_H
#include "linbox/linbox-config.h"
#include "linbox/util/timer.h"
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include "linbox/vector/blas-vector.h"
#include <utility>
#include <stdlib.h>
#include "linbox/util/commentator.h"
#include "linbox/util/commentator.h"

namespace LinBox
{

        /// No doc.
        /// @ingroup CRA
	template<class CRABase>
	struct ChineseRemainderSequential {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;

	public:
		const int MAXSKIP = 1000;
		const int MAXNONCOPRIME = 1000;

	protected:
		CRABase Builder_;
		int ngood_ = 0;
		int nbad_ = 0;
		int nskip_ = 0;

		/** \brief Helper class to sample unique primes.
		*/
		template <class PrimeIterator, bool is_unique = PrimeIterator::UniqueSamplingTag::value>
		struct PrimeSampler {
			const ChineseRemainderSequential& outer_;
			PrimeIterator& primeiter_;

			PrimeSampler (const ChineseRemainderSequential& outer, PrimeIterator& primeiter) :
				outer_(outer), primeiter_(primeiter)
			{ }

			/*! \brief Returns the next coprime element from the iterator.
			 */
			decltype(*primeiter_) operator() () {
				if (outer_.ngood_ == 0) return *primeiter_;
				int coprime = 0;
				while (outer_.Builder_.noncoprime(*primeiter_)) {
					++primeiter_;
					++coprime;
					if (coprime > outer_.MAXNONCOPRIME) {
						commentator().report(Commentator::LEVEL_ALWAYS,INTERNAL_ERROR) << "you are running out of primes. " << outer_.iterCount() << " used and " << coprime << " coprime primes tried for a new one.";
						throw LinboxError("LinBox ERROR: ran out of primes in CRA\n");
					}
				}
				return *primeiter_;
			}
		};

		/** \brief Call this when a bad prime is skipped.
		 */
		void doskip() {
			commentator().report(Commentator::LEVEL_IMPORTANT,INTERNAL_WARNING) << "bad prime, skipping\n";
			++nbad_;
			if (++nskip_ > MAXSKIP) {
				commentator().report(Commentator::LEVEL_ALWAYS,INTERNAL_ERROR) << "you are running out of GOOD primes. " << ngood_ << " good primes and " << nbad_ << " bad primes with " << nskip_ << " skipped in a row.\n";
				throw LinboxError("LinBox ERROR: ran out of good primes in CRA\n");
			}
		}

		/** \brief Gets a prime from the iterator that is coprime to the curent modulus.
		 */
		template <class PrimeIterator>
		inline auto get_coprime(PrimeIterator& primeiter) const -> decltype(*primeiter) {
			return PrimeSampler<PrimeIterator>(*this, primeiter)();
		}

	public:
		/** \brief Pass-through constructor to create the underlying builder.
		 */
		template <typename... Args>
		ChineseRemainderSequential(Args&&... args) :
			Builder_(std::forward<Args>(args)...)
		{ }

		/** \brief How many iterations have been performed so far.
		 *
		 * (This used to be stored in the public field IterCounter.)
		 */
		int iterCount() const {
			return ngood_ + nbad_;
		}

            /** \brief The \ref CRA loop
             *
             * Given a function to generate residues \c mod a single prime,
             * this loop produces the residues resulting from the Chinese
             * remainder process on sufficiently many primes to meet the
             * termination condition.
			 *
             * \param[out] res  an integer
             *
             * \param Iteration  Function object of two arguments, \c
             * Iteration(r, F), given prime field \p F it sets \p r
			 * to the residue(s) and returns an IterationResult
			 * to indicate how to incorporate the new residue.
             * This loop may be parallelized.  \p
             * Iteration  must be reentrant, thread safe. For example, \p
             * Iteration may be returning the coefficients of the minimal
             * polynomial of a matrix \c mod \p F.
             *
             * \param primeiter  iterator for generating primes.
             */
		template<class ResultType, class Function, class PrimeIterator>
		ResultType& operator() (ResultType& res, Function& Iteration, PrimeIterator& primeiter)
            {
                commentator().start ("Givaro::Modular iteration", "mmcravit");
				(*this)(-1, res, Iteration, primeiter);
                commentator().stop ("done", NULL, "mmcravit");
				return res;
            }

            /** \brief Run the CRA loop a certain number of times.
			 *
			 * This runs the CRA loop up to k times, or until termination
			 * if k is negative.
			 *
			 * \param k  maximum number of iterations, or run until termination
			 * if k is negative.
			 *
             * \param[out] res  an integer
             *
             * \param Iteration  Function object of two arguments, \c
             * Iteration(r, F), given prime field \p F it sets \p r
			 * to the residue(s) and returns an IterationResult
			 * to indicate how to incorporate the new residue.
             * This loop may be parallelized.  \p
             * Iteration  must be reentrant, thread safe. For example, \p
             * Iteration may be returning the coefficients of the minimal
             * polynomial of a matrix \c mod \p F.
             *
             * \param primeiter  iterator for generating primes.
             */
		template<class ResultType, class Function, class PrimeIterator>
		bool operator() (int k, ResultType& res, Function& Iteration, PrimeIterator& primeiter)
            {
				while (k != 0 && ngood_ == 0) {
					--k;
					Domain D(*primeiter);
                    commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
					++primeiter;
					auto r = CRAResidue<ResultType,Function>::create(D);
#ifdef __LB_CRA_TIMING__
                    Timer chrono; chrono.start();
#endif
					if (Iteration(r,D) == IterationResult::SKIP) {
						doskip();
					}
					else {
						++ngood_;
						Builder_.initialize(D,r);
					}
#ifdef __LB_CRA_TIMING__
                    chrono.stop();
                    std::clog << "1st iter : " << chrono << std::endl;
#endif
				}

				while (k != 0 && ! Builder_.terminated()) {
					--k;
					Domain D(get_coprime(primeiter));
                    commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
					++primeiter;
					auto r = CRAResidue<ResultType,Function>::create(D);

					switch (Iteration(r, D)) {
					case IterationResult::CONTINUE:
						++ngood_;
						Builder_.progress(D, r);
						break;
					case IterationResult::SKIP:
						doskip();
						break;
					case IterationResult::RESTART:
						commentator().report(Commentator::LEVEL_IMPORTANT,INTERNAL_WARNING) << "previous primes were bad; restarting\n";
						nbad_ += ngood_;
						ngood_ = 1;
						Builder_.initialize(D, r);
						break;
					}
				}

                Builder_.result(res);
				return ngood_ > 0 && Builder_.terminated();
            }

		template<class Param>
		bool changeFactor(const Param& p)
            {
                return Builder_.changeFactor(p);
            }

		template<class Param>
		Param& getFactor(Param& p)
            {
                return Builder_.getFactor(p);
            }

		bool changePreconditioner(const Integer& f, const Integer& m=Integer(1))
            {
                return Builder_.changePreconditioner(f,m);
            }

		Integer& getModulus(Integer& m)
            {
                Builder_.getModulus(m);
                return m;
            }

		Integer& getResidue(Integer& m)
            {
                Builder_.getResidue(m);
                return m;
            }

		Integer& result(Integer& m)
            {
                Builder_.result(m);
                return m;
            }

		template<class Int, template <class, class> class Vect, template <class> class Alloc >
		Vect<Int, Alloc<Int> >& result(Vect<Int, Alloc<Int> >& m)
            {
                Builder_.result(m);
                return m;
            }

#ifdef __LB_CRA_TIMING__
		inline std::ostream& reportTimes(std::ostream& os)
            {
                os <<  "Iterations:" << iterCount() << "\n";
                Builder_.reportTimes(os);
                return os;
            }
#endif

	};

	/** \brief Helper class to sample unique primes.
	 *
	 * This is the specialization for prime iterators that are already
	 * guaranteed to return unique primes (so that no checking is necessary).
	*/
	template <class CRABase>
	template <class PrimeIterator>
	struct ChineseRemainderSequential<CRABase>::PrimeSampler<PrimeIterator,true> {
		PrimeIterator& primeiter_;

		PrimeSampler (const ChineseRemainderSequential<CRABase>&, PrimeIterator& primeiter) :
			primeiter_(primeiter)
		{ }

		/*! \brief Returns the next coprime element from the iterator.
			*/
		decltype(*primeiter_) operator() () {
			return *primeiter_;
		}
	};

#ifdef __LB_CRA_TIMING__
	class CRATimer {
	public:
		mutable Timer ttInit, ttIRecon, /* ttImaging, ttIteration,*/ ttOther;
		void clear() const
            {
                ttInit.clear();
                ttIRecon.clear();
                    //ttImaging.clear();
                    //ttIteration.clear();
                ttOther.clear();
            }
	};
#endif

}

#endif //__LINBOX_sequential_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
