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
            int Tile = 8; //A magic number to boost the performance for unknown reason
			size_t NN = omp_get_max_threads();
			//std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);

            this->para_compute(Tile,NN,Iteration);

			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << this->IterCounter << " primes." << std::endl;
			return this->Builder_.result(res);
		}

		template<class Container, class Function, class PrimeIterator>
		Container& operator() (Container& res, Function& Iteration, PrimeIterator& primeiter)
		{
            int Tile = 8; //A magic number to boost the performance for unknown reason

			size_t NN = omp_get_max_threads();
			//std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);

            this->para_compute(Tile,NN,Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << this->IterCounter << " primes." << std::endl;
			return this->Builder_.result(res);
		}

        template<class Function>
        void para_compute(int Tile, int NN, Function Iteration){ 			
                    typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
                    std::set<int> coprimeset;
                    std::vector<ElementContainer> ROUNDresidues;ROUNDresidues.resize(omp_get_max_threads());
                    std::vector<Domain> ROUNDdomains;ROUNDdomains.resize(omp_get_max_threads());
                    std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;
                    //std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag>> m_primeiters;
                    
                    for(auto j=0;j<omp_get_max_threads();j++){
                        LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> m_primeiter( j, omp_get_max_threads());
                        //LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag> m_primeiter( j, omp_get_max_threads());
                        m_primeiters.push_back(m_primeiter);

                    }

                    if(omp_get_max_threads()>50){

                        early_termination_compute_task( (this->Builder_), NN,  Tile, m_primeiters, coprimeset, Iteration,  ROUNDdomains, ROUNDresidues);
                    
                    }else{
                        
                        compute_task( (this->Builder_), NN,  Tile, m_primeiters, coprimeset, Iteration,  ROUNDdomains, ROUNDresidues);
                        
                    }
        }

		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
			//! @bug why why why ???
			/** erreur: ‘omp_get_max_threads’ has not been declared
			 * ../linbox/algorithms/cra-domain-omp.h:152:16: note: suggested alternative:
			 * /usr/lib/gcc/x86_64-linux-gnu/4.6/include/omp.h:64:12: note:   ‘Givaro::omp_get_max_threads’
			 */
            int Tile = 8; //A magic number to boost the performance for unknown reason
			size_t NN = omp_get_max_threads();
			//std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res, den, Iteration,primeiter);

			int coprime =0;
			int maxnoncoprime = 1000;

            this->para_compute(Tile,NN,Iteration);			
			
			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << this->IterCounter << " primes." << std::endl;
			return this->Builder_.result(res,den);
		}
		

        template< class Function, class PrimeIterator, class Domain, class ElementContainer>
        void solve_with_prime(std::vector<PrimeIterator>& m_primeiters, std::set<int>& coprimeset, Function& Iteration, std::vector<Domain>& ROUNDdomains, std::vector<ElementContainer>& ROUNDresidues){

            ++m_primeiters[ omp_get_thread_num()];

#pragma omp critical
            while(this->Builder_.noncoprime(*m_primeiters[ omp_get_thread_num()]) &&
                  coprimeset.find(*m_primeiters[omp_get_thread_num()])!=coprimeset.end())
                ++m_primeiters[ omp_get_thread_num()];
            
            ROUNDdomains[ omp_get_thread_num()] = Domain(*m_primeiters[ omp_get_thread_num()]);

            Iteration(ROUNDresidues[ omp_get_thread_num()], ROUNDdomains[ omp_get_thread_num()]);

        }
        
        
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, int NN, int Tile, std::vector<PrimeIterator>& m_primeiters, std::set<int>& coprimeset, Function& Iteration, std::vector<Domain>& ROUNDdomains, std::vector<ElementContainer>& ROUNDresidues)
        {
        

			while( ! this->Builder_.terminated() ) {
#pragma omp parallel num_threads(NN/Tile)  
                //for(auto j=0;j<NN/Tile;j++)
                {
  
//#pragma omp task                        
                    {
                            // std::cout<<"Coucou thread_num = "<<omp_get_thread_num()<<std::endl;
                        for(auto i=0; i<Tile; ){
                            
                            solve_with_prime(m_primeiters, coprimeset, Iteration, ROUNDdomains, ROUNDresidues);
                            
#pragma omp critical
                            if(coprimeset.size()>0){
                                //if(coprimeset.find(*m_primeiters[ omp_get_thread_num()])==coprimeset.end()){
                                this->Builder_.progress( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);
                                i++;                                        
                                coprimeset.insert(*m_primeiters[ omp_get_thread_num()]);
                                //}
                            }else{
                                
                                this->Builder_.initialize( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);
                                i++;
                                coprimeset.insert(*m_primeiters[ omp_get_thread_num()]);
                            }
                            
                        }
                    }

                }

            }
            
        }
        
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void early_termination_compute_task(pFunc& pF, int NN, int Tile, std::vector<PrimeIterator>& m_primeiters, std::set<int>& coprimeset, Function& Iteration, std::vector<Domain>& ROUNDdomains, std::vector<ElementContainer>& ROUNDresidues)
        {
            
			while( ! this->Builder_.terminated() ) {
                
#pragma omp parallel num_threads(NN/Tile)
                //for(auto j=0;j<NN/Tile;j++)
                {
                    
#pragma omp task
                    {
                        
                        for(auto i=0; i<Tile; ){
                            // Avoid unnecessary computation so as to terminate as early as possible
                            if( this->Builder_.terminated() ) break;

                            solve_with_prime(m_primeiters, coprimeset, Iteration, ROUNDdomains, ROUNDresidues);
                            
#pragma omp critical
                            if(coprimeset.size()>0){
                                //if(coprimeset.find(*m_primeiters[ omp_get_thread_num()])==coprimeset.end()){
                                this->Builder_.progress( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);
                                i++;                                        
                                coprimeset.insert(*m_primeiters[ omp_get_thread_num()]);
                                //}
                            }else{
                                
                                this->Builder_.initialize( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);
                                i++;
                                coprimeset.insert(*m_primeiters[ omp_get_thread_num()]);
                            }
                            
                        }
                    }
                    
                }
            }
            
        }

        
		template<class Container, class Function, class PrimeIterator>
		Container& operator()  (Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
            int Tile = 8; //A magic number to boost the performance for unknown reason

			int NN = Tile*omp_get_max_threads();//Maybe replace omp_get_max_threads() with user required number ?
			std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (omp_get_max_threads() == 1) return Father_t::operator()(res, den, Iteration,primeiter);
            
            this->para_compute(Tile,NN,Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << IterCounter << " primes." << std::endl;

			return this->Builder_.result(res,den);

		}

        
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
