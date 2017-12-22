/* linbox/algorithms/cra-domain-seq.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Time-stamp: <04 Dec 17 14:59:50 Jean-Guillaume.Dumas@imag.fr>
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

/*! @file algorithms/cra-domain-seq.h
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

namespace LinBox
{

	template<class Function, class Field> struct CRATemporaryVectorTrait {
		typedef BlasVector<Field> Type_t;
	};

	/** \brief Glorified typedef for the CRA type based on the result type.
	 */
	template <typename ResultType, typename Domain>
	struct CRAResidueType {
		using Type = typename ResultType::template rebind<Domain>::other;
	};

	template <typename Domain>
	struct CRAResidueType<Integer, Domain> {
		using Type = typename Domain::Element;
	};

        /// No doc.
        /// @ingroup CRA
	template<class CRABase>
	struct ChineseRemainderSeq {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
	protected:
		CRABase Builder_;

	public:
		int IterCounter;

		/** \brief Pass-through constructor to create the underlying builder.
		 */
		template <typename... Args>
		ChineseRemainderSeq(Args&&... args) :
			Builder_(std::forward<Args>(args)...),
			IterCounter(0)
		{ }

            /** \brief The \ref CRA loop
             *
             * Given a function to generate residues \c mod a single prime,
             * this loop produces the residues resulting from the Chinese
             * remainder process on sufficiently many primes to meet the
             * termination condition.
             *
             * \param Iteration  Function object of two arguments, \c
             * Iteration(r, F), given prime field \p F it outputs
             * residue(s) \p r. This loop may be parallelized.  \p
             * Iteration  must be reentrant, thread safe. For example, \p
             * Iteration may be returning the coefficients of the minimal
             * polynomial of a matrix \c mod \p F.
             *
             * @warning  We won't detect bad primes.
             *
             * \param[out] res  an integer
             *
             * \param primeiter  iterator for generating primes.
             */
		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter)
            {
                commentator().start ("Givaro::Modular iteration", "mmcrait");
				(*this)(-1, res, Iteration, primeiter);
                commentator().stop ("mmcrait");
				return res;
			}

		template<class Iterator, class Function, class PrimeIterator>
		Iterator& operator() (Iterator& res, Function& Iteration, PrimeIterator& primeiter)
            {
                commentator().start ("Givaro::Modular vectorized iteration", "mmcravit");
				(*this)(-1, res, Iteration, primeiter);
                commentator().stop ("done", NULL, "mmcravit");
				return res;
            }

            /*
             *progress for k iterations
             */
		template<class ResultType, class Function, class PrimeIterator>
		bool operator() (int k, ResultType& res, Function& Iteration, PrimeIterator& primeiter)
            {
				using ResidueType = typename CRAResidueType<ResultType,Domain>::Type;
                if ((IterCounter ==0) && (k !=0)) {
                    --k;
                    ++IterCounter;
                    Domain D(*primeiter);
                    commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
                    ++primeiter;
					ResidueType r;
					// D.init(r); TODO is this needed???
#ifdef _LB_CRATIMING
                    Timer chrono; chrono.start();
#endif
                    Builder_.initialize( D, Iteration(r, D) );
#ifdef _LB_CRATIMING
                    chrono.stop();
                    std::clog << "1st iter : " << chrono << std::endl;
#endif
                }

                int coprime =0, nbprimes=0;
                int maxnoncoprime = 1000;

				while (k != 0 && ! Builder_.terminated()) {
					--k;
                    ++IterCounter;

                    while(Builder_.noncoprime(*primeiter) ) {
                        ++primeiter;
                        ++coprime;
                        if (coprime > maxnoncoprime) {
                            commentator().report(Commentator::LEVEL_ALWAYS,INTERNAL_ERROR) << "you are running out of primes. " << nbprimes << " used and " << maxnoncoprime << " coprime primes tried for a new one.";
                            return true;//term TODO why true, shouldn't it be false?
                        }
                    }

                    coprime =0;
                    Domain D(*primeiter);
                    commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
                    ++primeiter; ++nbprimes;

					ResidueType r;
					// D.init(r); TODO is this needed???
                    Builder_.progress( D, Iteration(r, D) );
                }
                Builder_.result(res);
				return Builder_.terminated();
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

#ifdef _LB_CRATIMING
		inline std::ostream& reportTimes(std::ostream& os)
            {
                os <<  "Iterations:" << IterCounter << "\n";
                Builder_.reportTimes(os);
                return os;
            }
#endif

	};

#ifdef _LB_CRATIMING
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
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
