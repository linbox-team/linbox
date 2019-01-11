/* Copyright (C) 2007 LinBox
 * Updated by Hongguang ZHU
 * Written by bds and zw
 * author: B. David Saunders and Zhendong Wan
 * parallelized for BOINC computing by Bryan Youse
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


#ifndef __LINBOX_cra_hybrid_H
#define __LINBOX_cra_hybrid_H

#define MPICH_IGNORE_CXX_SEEK //BB: ???
#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/rational-cra2.h"
#include "linbox/algorithms/rational-cra.h"
#include "linbox/util/mpicpp.h"


#include <unordered_set>
#include "linbox/randiter/random-prime.h"

#include "linbox/algorithms/cra-domain-omp.h"


namespace LinBox
{
    
	template<class CRABase>
	struct HybridChineseRemainder  {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
	protected:
		CRABase Builder_;
		Communicator* _commPtr;
		unsigned int _numprocs;
		double HB;//hadamard bound
        
	public:
		template<class Param>
		HybridChineseRemainder(const Param& b, Communicator *c) :
			Builder_(b), _commPtr(c), _numprocs(c->size())
			, HB(b)//Init with hadamard bound
		{}
        
		/** \brief The CRA loop.
		 *
		 * termination condition.
		 *
		 * \param Iteration  Function object of two arguments, \c
		 * Iteration(r, p), given prime \c p it outputs residue(s) \c
		 * r.  This loop may be parallelized.  \p Iteration must be
		 * reentrant, thread safe.  For example, \p Iteration may be
		 * returning the coefficients of the minimal polynomial of a
		 * matrix \c mod \p p.
		 @warning  we won't detect bad primes.
		 *
		 * \param primeg  RandIter object for generating primes.
		 * \param[out] res an integer
		 */
		template<class Function, class PrimeIterator>
		Integer & operator() (Integer& res, Function& Iteration, PrimeIterator& primeg)
		{
			//  defer to standard CRA loop if no parallel usage is desired
			if(_commPtr == 0 || _commPtr->size() == 1) {
				ChineseRemainder< CRABase > sequential(Builder_);
				return sequential(res, Iteration, primeg);
			}
			
			para_compute(res, Iteration, primeg);
			if(_commPtr->rank() == 0){
				return Builder_.result(res);
			}
			else{
                return res;
			}
		}
		
		template<class Function, class PrimeIterator>
		Integer & operator() (Integer& num, Integer& den, Function& Iteration, PrimeIterator& primeg)
		{

			//  defer to standard CRA loop if no parallel usage is desired
			if(_commPtr == 0 || _commPtr->size() == 1) {
				RationalRemainder< CRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}
			para_compute(num, Iteration, primeg);
			if(_commPtr->rank() == 0){
				return Builder_.result(num,den);
			}
			else{
                return num;
			}
		}


		template<class Vect, class Function, class PrimeIterator>
		Vect & operator() (Vect& num,  Integer& den, Function& Iteration, PrimeIterator& primeg)
		{
			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {

                RationalRemainder< CRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
                
			}
            para_compute(num, Iteration, primeg);
			if(_commPtr->rank() == 0){
				return Builder_.result(num,den);
			}
			else{
                return num;
			}
		}
		
		template<class Vect, class Function, class PrimeIterator>
		Vect & operator() (Vect& num, Function& Iteration, PrimeIterator& primeg)
		{
			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {

                ChineseRemainder< CRABase > sequential(Builder_);
				return sequential(num, Iteration, primeg);
                
			}
            para_compute(num, Iteration, primeg);
			if(_commPtr->rank() == 0){
				return Builder_.result(num);
			}
			else{
                return num;
			}
		}
#if 0
		template<class Function, class PrimeIterator>
		void  para_compute( Integer& res, Function& Iteration, PrimeIterator& primeg)
		{    
			int procs = _commPtr->size();
			int process = _commPtr->rank();
            std::unordered_set<int> prime_used;
			//  parent process
			if(process == 0 ){

				//  create an array to store primes
				int primes[procs - 1];
				DomainElement r;
				//  send each child process a new prime to work on
				for(int i=1; i<procs; i++){
					++primeg; 
					while(Builder_.noncoprime(*primeg) || prime_used.find(*primeg) != prime_used.end()  ) ++primeg;
					prime_used.insert(*primeg);
					primes[i - 1] = *primeg;
					_commPtr->send(primes[i - 1], i);
				}
				bool first_time = true;
				int poison_pills_left = procs - 1;
				//  loop until all execution is complete
				while( poison_pills_left > 0 ){
					int idle_process = 0;
					//  receive sub-answers from child procs
					_commPtr->recv(r, MPI_ANY_SOURCE);
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
					Domain D(primes[idle_process - 1]);
					//  assimilate results
					if(first_time){
						Builder_.initialize(D, r);
						first_time = false;
					}
					else
						Builder_.progress( D, r );
					//  queue a new prime if applicable
					if(! Builder_.terminated()){
						++primeg; 
						while(Builder_.noncoprime(*primeg) || prime_used.find(*primeg) != prime_used.end()  ) ++primeg;
    					prime_used.insert(*primeg);
						primes[idle_process - 1] = *primeg;
					}
					//  otherwise, queue a poison pill
					else{
						primes[idle_process - 1] = 0;
						poison_pills_left--;
					}
					//  send the prime or poison pill
					_commPtr->send(primes[idle_process - 1], idle_process);
				}  // end while

			}  // end if(parent process)
			//  child processes
			else{

				int pp;
				while(true){
					//  receive the prime to work on, stop
					//  if signaled a zero
					_commPtr->recv(pp, 0);
					if(pp == 0)
						break;
					Domain D(pp);
					DomainElement r; D.init(r);
					Iteration(r, D);
					//Comm->buffer_attach(rr);
					// send the results
					_commPtr->send(r, 0);
				}

			}
		}
#endif
		template<class Vect, class Function, class PrimeIterator>
		void  para_compute( Vect& num, Function& Iteration, PrimeIterator& primeg)
		{    
            
            Domain D(*primeg);
            BlasVector<Domain> r(D);
            Timer chrono;

			//  parent propcess
			if(_commPtr->rank() == 0){
               
                master_process_task(Iteration, D, r);

			}
			//  child process
			else{

                worker_process_task(Iteration, r);

			}

		}
		
///////////////////////////////////////////////////////////<-----------------Need to be multithreaded----------------
        
        template< class Function, class PrimeIterator, class Domain, class ElementContainer>
        void solve_with_prime(std::vector<PrimeIterator>& m_primeiters, std::set<int>& coprimeset,
                              Function& Iteration, std::vector<Domain>& ROUNDdomains,
                              ElementContainer& ROUNDresidues, //std::vector<ElementContainer>& ROUNDresidues, 
                              std::vector<CRABase>& vBuilders)
        {
            
            ++m_primeiters[ omp_get_thread_num()];

            while(vBuilders[ omp_get_thread_num()].noncoprime(*m_primeiters[ omp_get_thread_num()]) )
                ++m_primeiters[ omp_get_thread_num()];

            
            ROUNDdomains[ omp_get_thread_num()] = Domain(*m_primeiters[ omp_get_thread_num()]);
            
            Iteration(ROUNDresidues, ROUNDdomains[ omp_get_thread_num()]);
// std::cout<<" Thread("<<omp_get_thread_num()<<") : Worker("<<_commPtr->rank()<<") is inserting "<< *m_primeiters[ omp_get_thread_num()] << std::endl;
            ROUNDresidues.push_back(*m_primeiters[ omp_get_thread_num()]);//<------------------------
            
        }
        
        
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<PrimeIterator>& m_primeiters, std::set<int>& coprimeset,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders, size_t Niter)
        {
            
			size_t NN = Niter;//<=================Later modif required-------------
std::cout<<" Worker("<<_commPtr->rank()<<") found "<<omp_get_num_procs()<<" cores"<< std::endl;
   
#pragma omp parallel 
#pragma omp single
NN=omp_get_num_threads();

std::cout<<" Worker("<<_commPtr->rank()<<") has "<<NN<<" threads"<< std::endl;


#pragma omp parallel for num_threads(NN) schedule(dynamic,1)
            for(auto j=0;j<Niter;j++)
                {
  
                    solve_with_prime(m_primeiters, coprimeset, Iteration, ROUNDdomains, ROUNDresidues[j], vBuilders);

                }
 
            
        }


        template<class Vect, class Function>
        void worker_process_task(Function& Iteration,  Vect &r)
        {
           
            int pp=0;
            LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>   gen(_commPtr->rank(),_commPtr->size());
            //LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag>   gen(_commPtr->rank(),_commPtr->size());

            _commPtr->recv(pp, 0);

            if(pp!=0){
                std::unordered_set<int> prime_used;
std::cout<<" Worker("<<_commPtr->rank()<<") received Niter:="<<pp<< std::endl;
//<------------------------------could be replaced with OpenMP impl---------------------------
#if 1
			size_t NN = pp;//<=================Later modif required-------------
#pragma omp parallel 
#pragma omp single
NN=omp_get_num_threads();
            
            std::set<int> coprimeset;
            std::vector<BlasVector<Domain>> ROUNDresidues;ROUNDresidues.resize(pp);//<---------Niter dependant----------
            std::vector<Domain> ROUNDdomains;ROUNDdomains.resize(NN);
            std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;//<---------to be replace with generated primes from masked generator for a given process, no longer thread dependant which may cause duplication of some primes --------------------------------------------------------------------------------------------------
            std::vector<CRABase> vBuilders;vBuilders.resize(NN);
            
            for(auto j=0u;j<NN;j++){
                LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> m_primeiter( j, NN);
                
                CRABase Builder_(HB);vBuilders.push_back(Builder_);
                
                m_primeiters.push_back(m_primeiter);
                
            }
          
            compute_task( (this->Builder_), m_primeiters, coprimeset, Iteration,  ROUNDdomains,
                          ROUNDresidues, vBuilders, pp);	


                for(long i=0; i<pp; i++){

                    _commPtr->send(ROUNDresidues[i].begin(), ROUNDresidues[i].end(), 0, 0);

                 }

#endif
//------------------------------could be replaced with OpenMP impl--------------------------->

            };std::cout<<" Worker("<<_commPtr->rank()<<") Finished !!!!! "<< std::endl;

        }
///////////////////////////////////////////////////////////-----------------Need to be multithreaded---------------->

        template<class Vect>
        void compute_state_comm(int *primes, Vect &r, int &pp, int &idle_process, int &poison_pills_left)
        {
            
            idle_process = 0;
            
            r.resize (r.size()+1);

            //receive the beginnin and end of a vector in heapspace
            _commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0);
          
            //Dind out which process sent the solution and the coresponding prime number
            idle_process = (_commPtr->get_stat()).MPI_SOURCE;
            
            //Update the number of iterations for the next step
            poison_pills_left--;//poison_pills_left-=primes[idle_process - 1];
            
           
            //Store the corresponding prime number
            pp = r[r.size()-1];

            //Restructure the vector without added prime number
            r.resize (r.size()-1);
            
            
        }
        template<class Vect>
        void master_compute(int *primes, Vect &r, int Niter)
        {

            int poison_pills_left = Niter;//int poison_pills_left = _commPtr->size() - 1;
            int pp;
            int idle_process = 0;

            while(poison_pills_left > 0 ){
               
                compute_state_comm(primes, r, pp, idle_process, poison_pills_left);

                Domain D(pp);
                
                Builder_.progress(D, r);
                
                //primes[idle_process - 1] = (Builder_.terminated()) ? 1:0;

            }

        }
        
        template<class Vect, class Function>
        void master_init(int *primes, Function& Iteration, Domain &D, Vect &r, int &Niter)
        {
			int procs = _commPtr->size();

			LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>   gen(_commPtr->rank(),_commPtr->size());
            Niter=std::ceil(1.442695040889*HB/(double)(gen.getBits()-1));
std::cout<<" >>>>>>>>>>>>> Niter:="<< Niter << std::endl;
            //Compute nb of tasks ought to be realized for each process

            if(Niter<(procs-1)){

                for(long i=1; i<Niter+1; i++){
                    primes[i - 1] = 1;
                    _commPtr->send(primes[i - 1], i);             

                }
                for(long i=Niter+1; i<procs; i++){
                    primes[i - 1] = 0;
                    _commPtr->send(primes[i - 1], i);

                }

                }else{
                for(long i=1; i<Niter%(procs-1)+1; i++){
                    primes[i - 1] = Niter/(procs-1)+1;
                    _commPtr->send(primes[i - 1], i);

                }
                for(long i=Niter%(procs-1)+1; i<procs; i++){
                    primes[i - 1] = Niter/(procs-1);
                    _commPtr->send(primes[i - 1], i);
         
                }
            }

            
            //Initialize the buider and the receiver vector r
            Builder_.initialize( D, Iteration(r, D) );
        }
        
        template<class Vect, class Function>
        void master_process_task(Function& Iteration, Domain &D, Vect &r)
        {
            int primes[_commPtr->size() - 1];
            int Niter = 0;
            master_init(primes, Iteration, D, r, Niter);
            
            master_compute(primes, r, Niter);
   
        }



    };
    
}

#undef MPICH_IGNORE_CXX_SEEK
#endif // __LINBOX_cra_mpi_H
// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
