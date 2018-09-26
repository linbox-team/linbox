/* linbox/algorithms/cra-domain-omp.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Naive parallel chinese remaindering
 * Launch NN iterations in parallel, where NN=omp_get_max_threads()
 * Then synchronization and termintation test.
 * Time-stamp: <13 Mar 12 13:49:58 Jean-Guillaume.Dumas@imag.fr>
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

/*! @file algorithms/cra-domain-omp.h
 * @brief Parallel (OMP) version of \ref CRA
 * @ingroup CRA
 */

#ifndef __LINBOX_omp_cra_H
#define __LINBOX_omp_cra_H
// commentator is not thread safe
#define DISABLE_COMMENTATOR
#include <omp.h>
#include <set>
#include "linbox/algorithms/cra-domain-seq.h"

#include "linbox/algorithms/rational-cra.h"

namespace LinBox
{

	template<class CRABase>
	struct ChineseRemainderOMP : public ChineseRemainderSeq<CRABase> {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
		typedef ChineseRemainderSeq<CRABase>    Father_t;

		template<class Param>
		ChineseRemainderOMP(const Param& b) :
			Father_t(b)
		{}

		ChineseRemainderOMP(const CRABase& b) :
			Father_t(b)
		{}

		template <class ResultType, class Function, class PrimeIterator>
		ResultType& operator() (ResultType& res, Function& Iteration, PrimeIterator& primeiter)
		{
			using ResidueType = typename CRAResidue<ResultType,Function>::template ResidueType<Domain>;
			size_t NN = omp_get_max_threads();
			//std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);

			std::vector<Domain> ROUNDdomains; ROUNDdomains.reserve(NN);
			std::vector<ResidueType> ROUNDresidues; ROUNDresidues.reserve(NN);
			std::vector<IterationResult> ROUNDresults(NN);
			std::set<Integer> coprimeset;

			while (! this->Builder_.terminated()) {
				ROUNDdomains.clear();
				ROUNDresidues.clear();
				coprimeset.clear();

				while (coprimeset.size() < NN) {
					coprimeset.emplace(this->get_coprime(primeiter));
					++primeiter;
				}

				for(auto coprimesetiter = coprimeset.cbegin(); coprimesetiter != coprimeset.cend(); ++coprimesetiter) {
					// std::cerr << "With prime: " << *coprimesetiter << std::endl;
					ROUNDdomains.emplace_back(*coprimesetiter);
					ROUNDresidues.emplace_back(CRAResidue<ResultType,Function>::create(ROUNDdomains.back()));
				}

#pragma omp parallel for
				for(size_t i=0;i<NN;++i) {
					ROUNDresults[i] = Iteration(ROUNDresidues[i], ROUNDdomains[i]);
				}
#pragma omp barrier

				// if any thread says RESTART, then all CONTINUEs become SKIPs
				bool anyrestart = false;
				for (auto res : ROUNDresults) {
					if (res == IterationResult::RESTART) anyrestart = true;
				}
				if (anyrestart) {
					this->nbad_ += this->ngood_;
					this->ngood_ = 0;
				}

				for (size_t i = 0; i < NN; ++i) {
					if (ROUNDresults[i] == IterationResult::SKIP) {
						this->doskip();
					}
					else if (anyrestart && ROUNDresults[i] == IterationResult::CONTINUE) {
						// commentator should indicate that this prime is bad
						++this->nbad_;
					}
					else if (this->ngood_ == 0) {
						this->ngood_ = 1;
						this->Builder_.initialize(ROUNDdomains[i], ROUNDresidues[i]);
					}
					else {
						++this->ngood_;
						this->Builder_.progress(ROUNDdomains[i], ROUNDresidues[i]);
					}
				}
				//std::cerr << "Computed: " << iterCount() << " primes." << std::endl;
			}

			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << this->iterCount() << " primes." << std::endl;
			return this->Builder_.result(res);
		}
	};





	template<class CRABase>
	struct ChineseRemainderRatOMP : public RationalRemainder<CRABase> {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
		typedef RationalRemainder<CRABase>    Father_t;

		template<class Param>
		ChineseRemainderRatOMP(const Param& b) :
			Father_t(b)
		{}

		ChineseRemainderRatOMP(const CRABase& b) :
			Father_t(b)
		{}

		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
			//! @bug why why why ???
			/** erreur: ‘omp_get_max_threads’ has not been declared
			 * ../linbox/algorithms/cra-domain-omp.h:152:16: note: suggested alternative:
			 * /usr/lib/gcc/x86_64-linux-gnu/4.6/include/omp.h:64:12: note:   ‘Givaro::omp_get_max_threads’
			 */
			size_t NN = omp_get_max_threads();
			//std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);

			int coprime =0;
			int maxnoncoprime = 1000;

			if (this->IterCounter==0) {
				std::set<Integer> coprimeset;
				while(coprimeset.size() < NN) {
					++primeiter;
					while(this->Builder_.noncoprime(*primeiter) ) {
						++primeiter;
						++coprime;
						if (coprime > maxnoncoprime) {
							std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
							return this->Builder_.result(res);
						}
					}
					coprime =0;
					coprimeset.insert(*primeiter);
				}
				std::vector<Domain> ROUNDdomains; ROUNDdomains.reserve(NN);
				std::vector<DomainElement> ROUNDresidues(NN);
				typename std::vector<DomainElement>::iterator resit=ROUNDresidues.begin();
				for(std::set<Integer>::const_iterator coprimesetiter = coprimeset.begin(); coprimesetiter != coprimeset.end(); ++coprimesetiter,++resit) {
					// std::cerr << "With prime: " << *coprimesetiter << std::endl;
					ROUNDdomains.push_back( Domain(*coprimesetiter) );
					ROUNDdomains.back().init( *resit );
				}

#pragma omp parallel for
				for(size_t i=0;i<NN;++i) {

#pragma omp critical
					Iteration(ROUNDresidues[i], ROUNDdomains[i]);

				}
#pragma omp barrier
				++this->IterCounter;
				this->Builder_.initialize( ROUNDdomains[0],ROUNDresidues[0]);
				for(size_t i=1;i<NN;++i) {
					++this->IterCounter;
					this->Builder_.progress( ROUNDdomains[i],ROUNDresidues[i]);
				}
				// commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
			}

			while( ! this->Builder_.terminated() ) {
				//std::cerr << "Computed: " << this->IterCounter << " primes." << std::endl;
				std::set<Integer> coprimeset;
				while(coprimeset.size() < NN) {
					++primeiter;
					while(this->Builder_.noncoprime(*primeiter) ) {
						++primeiter;
						++coprime;
						if (coprime > maxnoncoprime) {
							std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
							return this->Builder_.result(res);
						}
					}
					coprime =0;
					coprimeset.insert(*primeiter);
				}
				std::vector<Domain> ROUNDdomains; ROUNDdomains.reserve(NN);
				std::vector<DomainElement> ROUNDresidues(NN);
				typename std::vector<DomainElement>::iterator resit=ROUNDresidues.begin();
				for(std::set<Integer>::const_iterator coprimesetiter = coprimeset.begin(); coprimesetiter != coprimeset.end(); ++coprimesetiter,++resit) {
					// std::cerr << "With prime: " << *coprimesetiter << std::endl;
					ROUNDdomains.push_back( Domain(*coprimesetiter) );
					ROUNDdomains.back().init( *resit );
				}

#pragma omp parallel for 
				for(size_t i=0;i<NN;++i) {

//#pragma omp critical
					Iteration(ROUNDresidues[i], ROUNDdomains[i]);

				}
#pragma omp barrier
				for(size_t i=0;i<NN;++i) {
					++this->IterCounter;
					this->Builder_.progress( ROUNDdomains[i],ROUNDresidues[i]);
				}
			}
			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << this->IterCounter << " primes." << std::endl;
			return this->Builder_.result(res);
		}


#if 1
		template<class Container, class Function, class PrimeIterator>
		Container& operator()  (Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{

			size_t NN = 8*omp_get_max_threads();
			std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (omp_get_max_threads() == 1) return Father_t::operator()(res, den,Iteration,primeiter);

			int coprime =0;

			long IterCounter=0;

			if (IterCounter==0) {
				std::set<Integer> coprimeset;
				while(coprimeset.size() < NN) {
					++primeiter;
					while(this->Builder_.noncoprime(*primeiter) ) {
						++primeiter;
						++coprime;

					}
					coprime =0;
					coprimeset.insert(*primeiter);
				}
				std::vector<Domain> ROUNDdomains; ROUNDdomains.reserve(NN);

				std::vector<DomainElement> ROUNDresidues(NN); 
//				typename std::vector<DomainElement>::iterator resit=ROUNDresidues.begin();

				for(std::set<Integer>::const_iterator coprimesetiter = coprimeset.begin(); coprimesetiter != coprimeset.end(); ++coprimesetiter) {

					ROUNDdomains.push_back( Domain(*coprimesetiter) );
				}


Iteration(ROUNDresidues[0], ROUNDdomains[0]);
				++IterCounter;
				this->Builder_.initialize( ROUNDdomains[0],ROUNDresidues[0]);
#pragma omp parallel for schedule(dynamic)
				for(size_t i=1;i<NN;++i) {

					Iteration(ROUNDresidues[i], ROUNDdomains[i]);
//					++IterCounter;
#pragma omp critical(ROUNDresidues)
					this->Builder_.progress( ROUNDdomains[i],ROUNDresidues[i]);
				}

				// commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
			}



			while( ! this->Builder_.terminated() ) {

				std::set<Integer> coprimeset;
				while(coprimeset.size() < NN) {
					++primeiter;
					while(this->Builder_.noncoprime(*primeiter) ) {
						++primeiter;
						++coprime;
					}

					coprime = 0;
					coprimeset.insert(*primeiter);
				}

				std::vector<Domain> ROUNDdomains; ROUNDdomains.reserve(NN);
				std::vector<DomainElement> ROUNDresidues(NN); 
//				typename std::vector<DomainElement>::iterator resit=ROUNDresidues.begin();

				for(std::set<Integer>::const_iterator coprimesetiter = coprimeset.begin(); coprimesetiter != coprimeset.end(); ++coprimesetiter) {

					ROUNDdomains.push_back( Domain(*coprimesetiter) );

				}



#pragma omp parallel for schedule(dynamic)
		 		for(size_t i=0;i<NN;++i) {

					Iteration(ROUNDresidues[i], ROUNDdomains[i]);
					//++IterCounter;
/*
if(omp_in_parallel())  std::cerr << "Thread("<<omp_get_thread_num()<<")Begin parallel executing >>>>>>>>>>> ("<<i<<")"<< std::endl;
else std::cerr << "Thread("<<omp_get_thread_num()<<") >>>>>>>>>>> ("<<i<<")"<< std::endl;
*/
#pragma omp critical(ROUNDresidues)
					this->Builder_.progress( ROUNDdomains[i],ROUNDresidues[i]);

/*
if(omp_in_parallel())  std::cerr << "Thread("<<omp_get_thread_num()<<")end parallel executing <<<<<<<<<<<<< ("<<i<<")"<< std::endl;
else std::cerr << "Thread("<<omp_get_thread_num()<<") <<<<<<<<<<<<< ("<<i<<")"<< std::endl;
*/

				}







			}

			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << IterCounter << " primes." << std::endl;
			return this->Builder_.result(res,den);
		}
#endif


	};



}

#endif //__LINBOX_omp_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
































































 






