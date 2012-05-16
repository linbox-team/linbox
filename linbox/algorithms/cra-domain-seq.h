/* linbox/algorithms/cra-domain-seq.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Time-stamp: <01 Apr 11 15:47:48 Jean-Guillaume.Dumas@imag.fr>
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
#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>
#include "linbox/util/commentator.h"

//$define _LB_CRATIMING

namespace LinBox
{

	template<class Function, class Element> struct CRATemporaryVectorTrait {
		typedef std::vector<Element> Type_t;
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

		template<class Param>
		ChineseRemainderSeq(const Param& b) :
			Builder_(b)
		{
			IterCounter=0;
		}
		ChineseRemainderSeq(const CRABase& b) :
			Builder_(b)
		{
			IterCounter=0;
		}

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
		 * \param primeiter  iterator for generating primes.
		 *
		 * \param[out] res  an integer
		 */
		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter)
		{
			commentator().start ("Modular iteration", "mmcrait");
			if (IterCounter==0) {
				++IterCounter;
				Domain D(*primeiter);
				std::ostream& report = commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
				report << "With prime " << *primeiter << std::endl;
				//commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
				++primeiter;
				DomainElement r; D.init(r);
				Builder_.initialize( D, Iteration(r, D) );
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			while( ! Builder_.terminated() ) {
				++IterCounter;
				while(Builder_.noncoprime(*primeiter) ) {
					++primeiter;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						return Builder_.result(res);
					}
				}
				coprime =0;
				Domain D(*primeiter);
				commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
                                ++primeiter;
				DomainElement r; D.init(r);
				Builder_.progress( D, Iteration(r, D) );
			}
			commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << IterCounter << " primes." << std::endl;
			return Builder_.result(res);
		}

		/*
		 * progress for k>=0 iterations
		 * run until terminated if k < 0
		 */
		template<class Function, class PrimeIterator>
		bool operator() (const int k, Integer& res, Function& Iteration, PrimeIterator& primeiter)
		{

			int i=0;
			if ((IterCounter ==0) && (k !=0)) {
				++i;
				++IterCounter;
				Domain D(*primeiter);
				commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
				++primeiter;
				DomainElement r; D.init(r);
				Builder_.initialize( D, Iteration(r, D) );
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			while ((k <0) || (i < k)) {
				if (Builder_.terminated()) break;
				++i;
				while(Builder_.noncoprime(*primeiter)) {
					++primeiter;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						return true ;//term
					}
				}
				coprime =0;

				Domain D(*primeiter);
				commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
				++primeiter;
				DomainElement r; D.init(r);
				Builder_.progress( D, Iteration(r, D) );
			}
			Builder_.result(res);
			if (Builder_.terminated() ) return true;
			else return false;
		}


		template<class Iterator, class Function, class PrimeIterator>
		Iterator& operator() (Iterator& res, Function& Iteration, PrimeIterator& primeiter)
		{
			commentator().start ("Modular vectorized iteration", "mmcravit");

			if (IterCounter==0) {
				Domain D(*primeiter);
				commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
				++primeiter;
				typename CRATemporaryVectorTrait<Function, DomainElement>::Type_t r;
				Builder_.initialize( D, Iteration(r, D) );
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			while( ! Builder_.terminated() ) {
				++IterCounter;
				while(Builder_.noncoprime(*primeiter) ) {
					++primeiter;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						return Builder_.result(res);
					}
				}

                                Domain D(*primeiter);
				commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
                                ++primeiter;
				typename CRATemporaryVectorTrait<Function, DomainElement>::Type_t r;
				Builder_.progress( D, Iteration(r, D) );
			}
			commentator().stop ("done", NULL, "mmcravit");
			return Builder_.result(res);
		}

		/*
		 *progress for k iterations
		 */

		template<class Iterator, class Function, class PrimeIterator>
		bool operator() (const int k, Iterator& res, Function& Iteration, PrimeIterator& primeiter)
		{

			int i=0;
			if ((IterCounter ==0) && (k !=0)) {
				++i;
				++IterCounter;
				Domain D(*primeiter);
				commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
				++primeiter;
				typename CRATemporaryVectorTrait<Function, DomainElement>::Type_t r;
				Builder_.initialize( D, Iteration(r, D) );
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			while( (k <0 ) || (i < k)) {
				if (Builder_.terminated()) break;
				++i;
				++IterCounter;

				while(Builder_.noncoprime(*primeiter) ) {
					++primeiter;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are runnig out of primes. " << maxnoncoprime << " coprime primes found";
						return true;//term
					}
				}

				coprime =0;
				Domain D(*primeiter);
				++primeiter;

				typename CRATemporaryVectorTrait<Function, DomainElement>::Type_t r;
				Builder_.progress( D, Iteration(r, D) );
			}
			Builder_.result(res);
			if (Builder_.terminated()) return true;
			else return false;
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

#undef _LB_CRATIMING

#endif //__LINBOX_sequential_cra_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

