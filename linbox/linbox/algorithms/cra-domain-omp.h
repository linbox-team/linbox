/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/cra-domain-omp.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Naive parallel chinese remaindering
 * Launch NN iterations in parallel, where NN=omp_get_max_threads()
 * Then synchronization and termintation test.
 * Time-stamp: <30 Mar 10 15:27:09 Jean-Guillaume.Dumas@imag.fr> 
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
#ifndef __LINBOX_OMP_CRA_H
#define __LINBOX_OMP_CRA_H
// commentator is not thread safe
#define DISABLE_COMMENTATOR
#include "linbox/algorithms/cra-domain-seq.h"

namespace LinBox {

    template<class CRABase>
    struct ChineseRemainderOMP : public ChineseRemainderSeq<CRABase> {
        typedef typename CRABase::Domain	Domain;
        typedef typename CRABase::DomainElement	DomainElement;
        template<class Param>
        ChineseRemainderOMP(const Param& b) : ChineseRemainderSeq<CRABase>(b) {}
        
	ChineseRemainderOMP(const CRABase& b) : ChineseRemainderSeq<CRABase>(b) {}
        
        template<class Function, class PrimeIterator>
        Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter) {
// commentator.start ("Parallel OMP Modular iteration", "mmcrait");
	    if (this->IterCounter==0) {
		++primeiter; 
		++this->IterCounter;
            	Domain D(*primeiter); 
// commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
            	DomainElement r; D.init(r);
            	this->Builder_.initialize( D, Iteration(r, D) );
            }

	    int coprime =0;
	    int maxnoncoprime = 1000;
std::cerr << "Blocs: " << omp_get_max_threads() << " iterations." << std::endl;

            while( ! this->Builder_.terminated() ) {
                size_t NN = omp_get_max_threads();
                std::vector<DomainElement> ROUNDresidues(NN);
                std::vector<Domain> ROUNDdomains(NN);
                for(size_t i=0;i<NN;++i) {
                    ++primeiter; ++this->IterCounter;
                    while(this->Builder_.noncoprime(*primeiter) ) {
			++primeiter;
			++coprime;
			if (coprime > maxnoncoprime) {
                            std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
				return this->Builder_.result(res);
			}
                    }
                    coprime =0;
                    ROUNDdomains[i] = Domain(*primeiter); 
// commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime[" << i << "]: " << *primeiter << std::endl;
                    ROUNDdomains[i].init(ROUNDresidues[i]);
                }
#pragma omp parallel for
                for(size_t i=0;i<NN;++i) {
                    Iteration(ROUNDresidues[i], ROUNDdomains[i]);
                }

#pragma omp barrier
                for(size_t i=0;i<NN;++i) {
                    this->Builder_.progress( ROUNDdomains[i],ROUNDresidues[i]);
                }
            }
// commentator.stop ("done", NULL, "mmcrait");
std::cerr << "Used: " << this->IterCounter << " primes." << std::endl;
            return this->Builder_.result(res);
        }

				
    };
}

#endif
