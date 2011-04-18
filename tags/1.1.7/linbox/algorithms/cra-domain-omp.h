/* linbox/algorithms/cra-domain-omp.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Naive parallel chinese remaindering
 * Launch NN iterations in parallel, where NN=omp_get_max_threads()
 * Then synchronization and termintation test.
 * Time-stamp: <16 Jul 10 16:59:17 Jean-Guillaume.Dumas@imag.fr> 
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#ifndef __LINBOX_omp_cra_H
#define __LINBOX_omp_cra_H
// commentator is not thread safe
#define DISABLE_COMMENTATOR
#include "linbox/algorithms/cra-domain-seq.h"
#include <set>
#include <vector>

namespace LinBox 
{

    template<class CRABase>
    struct ChineseRemainderOMP : public ChineseRemainderSeq<CRABase> {
        typedef typename CRABase::Domain	Domain;
        typedef typename CRABase::DomainElement	DomainElement;
        typedef ChineseRemainderSeq<CRABase>    Father_t;

        template<class Param>
        ChineseRemainderOMP(const Param& b) : Father_t(b) {}
        
	ChineseRemainderOMP(const CRABase& b) : Father_t(b) {}
        
        template<class Function, class PrimeIterator>
        Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter) {
            size_t NN = omp_get_max_threads();
            std::cerr << "Blocs: " << NN << " iterations." << std::endl;
// commentator.start ("Parallel OMP Modular iteration", "mmcrait");
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
                for(long i=0;i<NN;++i) {
                    Iteration(ROUNDresidues[i], ROUNDdomains[i]);
                }
#pragma omp barrier
                ++this->IterCounter;
            	this->Builder_.initialize( ROUNDdomains[0],ROUNDresidues[0]);
                for(size_t i=1;i<NN;++i) {
                    ++this->IterCounter;
                    this->Builder_.progress( ROUNDdomains[i],ROUNDresidues[i]);
                }
// commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
            }

            while( ! this->Builder_.terminated() ) {
std::cerr << "Computed: " << this->IterCounter << " primes." << std::endl;
                size_t NN = omp_get_max_threads();
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
                for(long i=0;i<NN;++i) {
                    Iteration(ROUNDresidues[i], ROUNDdomains[i]);
                }
#pragma omp barrier
                for(size_t i=0;i<NN;++i) {
                    ++this->IterCounter;
                    this->Builder_.progress( ROUNDdomains[i],ROUNDresidues[i]);
                }
            }
// commentator.stop ("done", NULL, "mmcrait");
std::cerr << "Used: " << this->IterCounter << " primes." << std::endl;
            return this->Builder_.result(res);
        }

	template<class Container, class Function, class PrimeIterator>
        Container& operator() (Container& res, Function& Iteration, PrimeIterator& primeiter) {
            typedef typename CRATemporaryVectorTrait<Function, DomainElement>::Type_t ElementContainer;
            size_t NN = omp_get_max_threads();
            std::cerr << "Blocs: " << NN << " iterations." << std::endl;
// commentator.start ("Parallel OMP Modular iteration", "mmcrait");
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
// commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
            }

            while( ! this->Builder_.terminated() ) {
std::cerr << "Computed: " << this->IterCounter << " primes." << std::endl;
                size_t NN = omp_get_max_threads();
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
// commentator.stop ("done", NULL, "mmcrait");
std::cerr << "Used: " << this->IterCounter << " primes." << std::endl;
            return this->Builder_.result(res);
        }


				
    };
}

#endif //__LINBOX_omp_cra_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
