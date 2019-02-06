/* linbox/algorithms/cra-domain-omp.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Updated by Hongguang Zhu <zhuhongguang2014@gmail.com>
 *
 *
 * Naive parallel chinese remaindering
 * Launch NN iterations in parallel, where NN=omp_get_num_threads()
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
        
        double B;
        
		template<class Param>
		ChineseRemainderOMP(const Param& b) :
			Father_t(b), B(b)
		{}
        
		ChineseRemainderOMP(const CRABase& b) :
			Father_t(b)
		{}
        
		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter)
		{
			//! @bug why why why ???
			/** erreur: ‘omp_get_num_threads’ has not been declared
			 * ../linbox/algorithms/cra-domain-omp.h:152:16: note: suggested alternative:
			 * /usr/lib/gcc/x86_64-linux-gnu/4.6/include/omp.h:64:12: note:   ‘Givaro::omp_get_num_threads’
			 */
            
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);
            
            this->para_compute(Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res);
		}
        
		template<class Container, class Function, class PrimeIterator>
		Container& operator() (Container& res, Function& Iteration, PrimeIterator& primeiter)
		{
            
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);
            
            this->para_compute(Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res);
		}
        
        
		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
			//! @bug why why why ???
			/** erreur: ‘omp_get_num_threads’ has not been declared
			 * ../linbox/algorithms/cra-domain-omp.h:152:16: note: suggested alternative:
			 * /usr/lib/gcc/x86_64-linux-gnu/4.6/include/omp.h:64:12: note:   ‘Givaro::omp_get_num_threads’
			 */
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res, den, Iteration,primeiter);
            
            this->para_compute(Iteration);			
			
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res,den);

		}
		
        template<class Function>
        void para_compute(Function Iteration){

            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; }

#if 1       //----------------------------------task based paladin--------------------------------------------
            typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
            //std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;
            std::vector<int> m_primeiters;
            LinBox::PrimeIterator<LinBox::IteratorCategories::DeterministicTag> gen;
            long Niter=std::ceil(1.442695040889*B/(double)(gen.getBits()-1));
            m_primeiters.reserve(Niter);
            std::vector<CRABase> vBuilders;
            vBuilders.reserve(NN);
            for(auto j=0;j<Niter;j++){
                ++gen;
                while(this->Builder_.noncoprime(*gen) )
                    ++gen;

                m_primeiters.push_back(*gen);
                
            }

            for(auto j=0;j<NN;j++){
               
                CRABase Builder_(B);
                vBuilders.push_back(Builder_);
                
            }
#else
            typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
            //std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;
            std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag>> m_primeiters;
           
            m_primeiters.reserve(NN);
            std::vector<CRABase> vBuilders;
            vBuilders.reserve(NN);


            for(auto j=0;j<NN;j++){
                //LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> m_primeiter( j, NN);
                LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag> m_primeiter( j, NN);
                m_primeiters.push_back(m_primeiter);
                
                CRABase Builder_(B);
                vBuilders.push_back(Builder_);
                
            }
            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
#endif

            
            
            std::vector<ElementContainer> ROUNDresidues;ROUNDresidues.resize(Niter);
            std::vector<Domain> ROUNDdomains;ROUNDdomains.resize(Niter);

            compute_task( (this->Builder_), m_primeiters, Iteration,  ROUNDdomains,
                          ROUNDresidues, vBuilders);
            
            
        }


#if 1       //----------------------------------task based paladin--------------------------------------------
        template< class Function, class Domain, class ElementContainer>
        void solve_with_prime(int m_primeiters,
                              Function& Iteration, Domain& ROUNDdomains,
                              ElementContainer& ROUNDresidues, CRABase& vBuilders)
        {
            
            ROUNDdomains = Domain(m_primeiters);

            Iteration(ROUNDresidues, ROUNDdomains);
            
        }
#else
        template< class Function, class PrimeIterator, class Domain, class ElementContainer>
        void solve_with_prime(PrimeIterator& m_primeiters,
                              Function& Iteration, Domain& ROUNDdomains,
                              ElementContainer& ROUNDresidues, CRABase& vBuilders)
        {

            ++m_primeiters;

            while(vBuilders.noncoprime(*m_primeiters) )
                ++m_primeiters;

            
            ROUNDdomains = Domain(*m_primeiters);

            Iteration(ROUNDresidues, ROUNDdomains);
            
        }
#endif



#if 1 //=============================================Task based ===========================================
#if 1 // ---------------------------- Paladin impl ----------------------------
/*
#define TASK(M, I)                             \
    PRAGMA_OMP_IMPL(omp task M)           \
    {I;}

#define FORBLOCK1D(iter, m, Helper, Args...)                                       \
    { FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param > iter(m, Helper); \
      for(iter.initialize(); !iter.isTerminated(); ++iter){ {Args;}  } }

#define FOR1D(i, m, Helper, Args...)                             \
    FORBLOCK1D(_internal_iterator, m, Helper,                           \
        TASK( , \
             {for(auto i=_internal_iterator.begin(); i!=_internal_iterator.end(); ++i) \
                 { Args; } });)                                         \
        WAIT;
        
// overload of SPLITTER
#define splitting_0() FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads>()
#define splitting_1(a) FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads>(a)
#define splitting_2(a,c) FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,c>(a)
#define splitting_3(a,b,c) FFLAS::ParSeqHelper::Parallel<b,c>(a)

// parallel region
#define PAR_BLOCK  PRAGMA_OMP_IMPL(omp parallel)   \
    PRAGMA_OMP_IMPL(omp single)
*/

        template<class pFunc, class Function, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<int>& m_primeiters,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {
            
            LinBox::PrimeIterator<LinBox::IteratorCategories::DeterministicTag> gen;
            long Niter=std::ceil(1.442695040889*B/(double)(gen.getBits()-1));
            
//std::cout<<"============================= Niter:= "<<Niter<<" ================================"<<std::endl;
/*@Try PlanB as follows:
*  As the Niter is known then first generate Niter prime numbers and put them into a vector then start paladin impl
*/


            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 
            
//            FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Grain> H;
splitting_3(H,FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Grain);


                PAR_BLOCK{

	                FOR1D(i,Niter,H,
{
//std::cout<<"i = "<<i<<" == "<<"  _internal_iterator.blockindex()  : "<<"_internal_iterator.blockindex()%NN = "<<_internal_iterator.blockindex()%NN<<std::endl;
   std::cout << "Threads: " << NUM_THREADS << ", max: " << MAX_THREADS << std::endl;
   std::cerr << "OMP: " << __FFLASFFPACK_USE_OPENMP << ", max " << omp_get_max_threads() << std::endl;
                      solve_with_prime(
                            m_primeiters[i],
                            Iteration,
                            ROUNDdomains[i],
                            ROUNDresidues[i],
                            vBuilders[_internal_iterator.blockindex()%NN]
                        );
}
                    )

                }



            this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

            for(auto j=1;j<Niter;j++)
                {

                        this->Builder_.progress( ROUNDdomains[j], ROUNDresidues[j]);

                }
           
        }
#else // ---------------------------- OMP impl ----------------------------
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<PrimeIterator>& m_primeiters,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {
            
            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
            int NN;
#pragma omp parallel
#pragma omp single 
            NN=NUM_THREADS;

solve_with_prime(m_primeiters[0], Iteration, ROUNDdomains[0], ROUNDresidues[0], vBuilders[0]);
this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

#pragma omp parallel num_threads(NN)
#pragma omp single
	        for(auto j=1;j<Niter;j++){
#pragma omp task
                solve_with_prime(m_primeiters[omp_get_thread_num()], Iteration, ROUNDdomains[j], ROUNDresidues[j], vBuilders[omp_get_thread_num()]);
                
                }


//            this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

            for(auto j=1;j<Niter;j++)
                {

                        this->Builder_.progress( ROUNDdomains[j], ROUNDresidues[j]);

                }
           
        }
#endif  // --------------------------------------------------------------------------

#else   //=============================================For loop based ===========================================
#if 0   // ---------------------------- Paladin impl ----------------------------
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<PrimeIterator>& m_primeiters,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {
            
            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 
            
            FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Grain> H;

            if(NN>Niter){
SYNCH_GROUP(
	            PARFORBLOCK1D(k,Niter,H,
                    solve_with_prime(m_primeiters[k], Iteration, ROUNDdomains[k], ROUNDresidues[k], vBuilders[k]);
                );
)

            }else{
SYNCH_GROUP(
	            PARFORBLOCK1D(k,NN,H,
	                for(auto j=k*(Niter/NN);j<(k+1)*(Niter/NN);j++)
                        solve_with_prime(m_primeiters[k], Iteration, ROUNDdomains[j], ROUNDresidues[j], vBuilders[k]);
                );
)
                if(Niter%NN>0){
SYNCH_GROUP(
	                PARFORBLOCK1D(k,Niter%NN,H,
                        solve_with_prime(m_primeiters[k], Iteration, ROUNDdomains[k+NN*(Niter/NN)], ROUNDresidues[k+NN*(Niter/NN)], vBuilders[k]);
                    );
)
                }

            }


            this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

            for(auto j=1;j<Niter;j++)
                {

                        this->Builder_.progress( ROUNDdomains[j], ROUNDresidues[j]);

                }
           
        }
#else   // ---------------------------- OMP impl ----------------------------
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<PrimeIterator>& m_primeiters,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {
            
            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
            int NN;
#pragma omp parallel
#pragma omp single 
            NN=NUM_THREADS;

solve_with_prime(m_primeiters[0], Iteration, ROUNDdomains[0], ROUNDresidues[0], vBuilders[0]);
this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

#pragma omp parallel for schedule(dynamic,1) num_threads(NN)
	        for(auto j=1;j<Niter;j++){
                solve_with_prime(m_primeiters[omp_get_thread_num()], Iteration, ROUNDdomains[j], ROUNDresidues[j], vBuilders[omp_get_thread_num()]);
                
                }


//            this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

            for(auto j=1;j<Niter;j++)
                {

                        this->Builder_.progress( ROUNDdomains[j], ROUNDresidues[j]);

                }
           
        }
        

#endif  // ---------------------------------------------------------------
#endif  //=============================================================================================
        
		template<class Container, class Function, class PrimeIterator>
		Container& operator()  (Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
            
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res, den, Iteration,primeiter);
            
            this->para_compute(Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
            
			//return this->Builder_.result(res,den);
this->Builder_.result(res,den);
//std::cerr << ">>>>res: " << std::endl; for(long j=0;j<res.size();j++) std::cerr << res.getEntry(j) << std::endl; 
return res;
            
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
