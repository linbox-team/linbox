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
#include "linbox/randiter/random-prime.h"
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

		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter)
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

		template<class Container, class Function, class PrimeIterator>
		Container& operator() (Container& res, Function& Iteration, PrimeIterator& primeiter)
		{
			typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
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
				std::vector<ElementContainer> ROUNDresidues(NN);
				typename std::vector<ElementContainer>::iterator resit=ROUNDresidues.begin();
				for(std::set<Integer>::const_iterator coprimesetiter = coprimeset.begin(); coprimesetiter != coprimeset.end(); ++coprimesetiter,++resit) {
					// std::cerr << "With prime: " << *coprimesetiter << std::endl;
					ROUNDdomains.push_back( Domain(*coprimesetiter) );
					//                     for(typename ElementContainer::iterator reselit = resit->begin();
					//                         reselit != resit->end(); ++reselit)
					//                         ROUNDdomains.back().init( *reselit );
				}
                
#pragma omp parallel for
				for(size_t i=0;i<NN;++i) {
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
				std::vector<ElementContainer> ROUNDresidues(NN);
				typename std::vector<ElementContainer>::iterator resit=ROUNDresidues.begin();
				for(std::set<Integer>::const_iterator coprimesetiter = coprimeset.begin(); coprimesetiter != coprimeset.end(); ++coprimesetiter,++resit) {
					// std::cerr << "With prime: " << *coprimesetiter << std::endl;
					ROUNDdomains.push_back( Domain(*coprimesetiter) );
					//                     for(typename ElementContainer::iterator reselit = resit->begin();
					//                         reselit != resit->end(); ++reselit)
					//                         ROUNDdomains.back().init( *reselit );
				}
                
#pragma omp parallel for
				for(size_t i=0;i<NN;++i) {
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
            int Tile = 8; //A magic number to boost the performance for unknown reason
			typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
			int NN = Tile*omp_get_max_threads();
			std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (omp_get_max_threads() == 1) return Father_t::operator()(res, den,Iteration,primeiter);

            std::set<int> coprimeset;
            std::vector<ElementContainer> ROUNDresidues;ROUNDresidues.resize(omp_get_max_threads());
            std::vector<Domain> ROUNDdomains;ROUNDdomains.resize(omp_get_max_threads());
            std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;

            for(auto j=0;j<omp_get_max_threads();j++){
                LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> m_primeiter( j,omp_get_max_threads(),25);
                m_primeiters.push_back(m_primeiter);
            }

			while( ! this->Builder_.terminated() ) {
                
#pragma omp parallel for num_threads(NN/Tile)
                for(auto j=0;j<NN/Tile;j++)
                    {

#pragma omp task                        
                        {

                            for(auto i=0; i<NN; ){
                                // Avoid unnecessary computation so as to terminate as early as possible
                                if( this->Builder_.terminated() ) break; 
                              ++m_primeiters[ omp_get_thread_num()]; while(this->Builder_.noncoprime(*m_primeiters[ omp_get_thread_num()]) && coprimeset.find(*m_primeiters[omp_get_thread_num()])!=coprimeset.end()) ++m_primeiters[ omp_get_thread_num()];     

                                ROUNDdomains[ omp_get_thread_num()] = Domain(*m_primeiters[ omp_get_thread_num()]);
                                
                                Iteration(ROUNDresidues[ omp_get_thread_num()], ROUNDdomains[ omp_get_thread_num()]);
                           
#pragma omp critical
                                if(coprimeset.size()>0){
//                                    if(coprimeset.find(*m_primeiters[ omp_get_thread_num()])==coprimeset.end()){
                                        this->Builder_.progress( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);
                                        i++;                                        
                                        coprimeset.insert(*m_primeiters[ omp_get_thread_num()]);
//                                    }
                                }else{
                                    
                                    this->Builder_.initialize( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);
                                    i++;
                                    coprimeset.insert(*m_primeiters[ omp_get_thread_num()]);
                                }

                            }
                        }

                    }
                
			}
            
			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << IterCounter << " primes." << std::endl;
			return this->Builder_.result(res,den);
		}
        
#else
		template<class Container, class Function, class PrimeIterator>
		Container& operator()  (Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
            int Tile = 8;
			typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
			auto NN = Tile*omp_get_max_threads();
			std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (omp_get_max_threads() == 1) return Father_t::operator()(res, den,Iteration,primeiter);
            
            std::set<int> coprimeset;

            
			while( ! this->Builder_.terminated() ) {
                
                
                
#pragma omp parallel for num_threads(NN/Tile) 
                for(auto j=0;j<NN/Tile;j++)
                    {                        
                        
#pragma omp task                        
                        {
                            ElementContainer	ROUNDres;
                            Domain			ROUNDdom;
                            PrimeIterator m_primeiter; ++m_primeiter;
                            
                            for(auto i=0; i<NN; ){
                                // Avoid unnecessary computation so as to terminate as early as possible
                                if( this->Builder_.terminated() ) break; 
                                                                
                                while(this->Builder_.noncoprime(*m_primeiter) ) ++m_primeiter;
                                ROUNDdom = Domain(*m_primeiter);
                                Iteration(ROUNDres, ROUNDdom);
                                
                                
#pragma omp critical
                                if(coprimeset.size()>0){
                                    
                                    if(coprimeset.find(*m_primeiter)==coprimeset.end()){
                                        this->Builder_.progress( ROUNDdom, ROUNDres);
                                        i++;
                                        coprimeset.insert(*m_primeiter);
                                        
                                    }
                                    
                                }else{
                                    
                                    this->Builder_.initialize( ROUNDdom, ROUNDres);
                                    i++;
                                    coprimeset.insert(*m_primeiter);
                                }
                                
                            }
                        }
                        
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
