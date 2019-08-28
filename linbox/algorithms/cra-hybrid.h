/* Copyright (C) 2007 LinBox
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

#pragma once

#include <unordered_set>
#include <utility>
#include <vector>

#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/rational-cra.h"
#include "linbox/algorithms/rational-cra-var-prec.h"
#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/mpicpp.h"
#include "linbox/util/timer.h"

#include <unordered_set>
#include "linbox/randiter/random-prime.h"

#include "linbox/algorithms/cra-domain-omp.h"

#if defined(__LINBOX_HAVE_MPI)
namespace LinBox
{
    /* @Warning /!\
     * The MPI standard imposes that if the following program is launched with only 2 processes, multithreading
     * will not be available.
     * Multithreading is only active if launched with more than 2 processes while each process is mapped to
     * one node.
     * To use multithreading, user needs to set the environment useing export for the number of cores/threads on
     * each node if each process is mapped to a node however this only works for one run so this will not work
     * for multiple executions where only the -t option is needed to set the number of threads
     *
     * For the execution on several nodes, each node needs to have at least 3 processes mapped otherwise no
     * multithreading will be available according to the MPI standards
     */
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

		int getNiter(){
		    return std::ceil(1.442695040889*HB/(double)(LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>(0,_commPtr->size()).getBits()-1));
		}

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
				RationalChineseRemainder< CRABase > sequential(Builder_);
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
			//  then proceed sequentially (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {

                RationalChineseRemainder< CRABase > sequential(Builder_);
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
			//  then proceed sequentially (without parallel)
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


        template< class Function, class ElementContainer>
        void solve_with_prime(int m_primeiter,
                              Function& Iteration,
                              ElementContainer& VECTORresidues
                              )
        {
            Domain D(m_primeiter);
            //@warning: The commentator within the following function's function call hierarchy will raise a crash if the whole program is not compiled with all commentator disabled as the commentator is not thread safe
            Iteration(VECTORresidues, D );
            VECTORresidues.push_back(m_primeiter);
        }


        template<class pFunc, class Function, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<int>& m_primeiters,
                          Function& Iteration,
                          std::vector<ElementContainer>& VECTORresidues, size_t Ntask)
        {
            //@fixme: cannot use export OMP_NUM_THREADS to set the desired number of threads for more than one iterations and only omp_set_num_thread() can be used for this purpose, whick takes the option -t
#if 0
            PAR_BLOCK{
                auto sp=SPLITTER(NUM_THREADS,FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Threads);
                SYNCH_GROUP({
                    FORBLOCK1D(iter, Ntask, sp,{
                        TASK(MODE(CONSTREFERENCE(m_primeiters,Iteration,VECTORresidues)),{
                                for(auto j=iter.begin(); j!=iter.end(); ++j)
                                {
                                    solve_with_prime(m_primeiters[j], Iteration, VECTORresidues[j]);
                                }
                         })
                     });
                });
            }
#else
            PAR_BLOCK{
                auto sp=SPLITTER(NUM_THREADS,FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Threads);
                SYNCH_GROUP({
                    FOR1D(iter, Ntask, sp,MODE(CONSTREFERENCE(m_primeiters,Iteration,VECTORresidues)),{
                                    solve_with_prime(m_primeiters[iter], Iteration, VECTORresidues[iter]);
                     });
                });
            }
#endif
        }


        template<class Vect, class Function>
        void worker_process_task(Function& Iteration,  Vect &r)
        {
            int Ntask=0;
            //LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>   gen(_commPtr->rank(),_commPtr->size());
            LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag>   gen(_commPtr->rank(),_commPtr->size());
            ++gen;
            _commPtr->recv(Ntask, 0);

            if(Ntask!=0){
                std::unordered_set<int> prime_used;

                std::vector<BlasVector<Domain>> VECTORresidues;VECTORresidues.resize(Ntask);
                std::vector<int> m_primeiters;m_primeiters.reserve(Ntask);

                for(auto j=0;j<Ntask;j++){
                    ++gen;
                    while(this->Builder_.noncoprime(*gen) || prime_used.find(*gen)!=prime_used.end())
                        ++gen;
                    m_primeiters.push_back(*gen);
                    prime_used.insert(*gen);
                }

                compute_task( (this->Builder_), m_primeiters, Iteration, VECTORresidues, Ntask);

                for(long i=0; i<Ntask; i++){
                    _commPtr->send(VECTORresidues[i].begin(), VECTORresidues[i].end(), 0, 0);
                 }


            };

        }

	template<class Vect, class Function, class PrimeIterator>
		void  para_compute( Vect& num, Function& Iteration, PrimeIterator& primeg)
		{

            Domain D(*primeg);
            BlasVector<Domain> r(D);

#ifdef __Detailed_Time_Measurement
            double starttime;
            double endtime;
#endif

			//  parent propcess
			if(_commPtr->rank() == 0){

                master_process_task(Iteration, D, r);

			}
			//  child process
			else{
#ifdef __Detailed_Time_Measurement
		        chrono.start();
#endif
                worker_process_task(Iteration, r);
#ifdef __Detailed_Time_Measurement
		        chrono.stop();
                std::cout<<" process("<<_commPtr->rank()<<") used total CPU time (seconds): " << endtime-starttime<<std::endl;
#endif
			}

		}

        template<class Vect>
        void master_recv_residues(Vect &r, int &pp, int &Nrecv)
        {
            r.resize (r.size()+1);

           //receive the beginnin and end of a vector in heapspace
            _commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0);

            //Update the number of iterations for the next step
            Nrecv--;

            //Store the corresponding prime number
            pp = r[r.size()-1];

            //Restructure the vector without added prime number
            r.resize (r.size()-1);

        }

        template<class Vect>
        void master_compute(Vect &r)
        {

            int pp;

#ifdef __Detailed_Time_Measurement
            Timer chrono;
#endif
            int Nrecv=this->getNiter();

            while(Nrecv > 0 ){

                master_recv_residues(r, pp, Nrecv);

                Domain D(pp);

#ifdef __Detailed_Time_Measurement
		chrono.start();
#endif
                Builder_.progress(D, r);
#ifdef __Detailed_Time_Measurement
		chrono.stop();
		std::cout<<"Builder_.progress(D, r) in the manager process used CPU time (seconds): " <<chrono.usertime()<<std::endl;
#endif
            }

        }


        template<class Vect, class Function>
        void master_process_task(Function& Iteration, Domain &D, Vect &r)
        {
            int vNtask_per_proc[_commPtr->size() - 1];

            master_init(vNtask_per_proc, Iteration, D, r);

            master_compute(r);

        }

        template<class Vect, class Function>
        void master_init(int *vNtask_per_proc, Function& Iteration, Domain &D, Vect &r)
        {
			int procs = _commPtr->size();

            int Niter=this->getNiter();

            //Compute and send nb of tasks ought to be realized for each worker process
            if(Niter<(procs-1)){

                for(long i=1; i<Niter+1; i++){
                    vNtask_per_proc[i - 1] = 1;
                    _commPtr->send(vNtask_per_proc[i - 1], i);
                }
                for(long i=Niter+1; i<procs; i++){
                    vNtask_per_proc[i - 1] = 0;
                    _commPtr->send(vNtask_per_proc[i - 1], i);
                }

             }else{
                for(long i=1; i<Niter%(procs-1)+1; i++){
                    vNtask_per_proc[i - 1] = Niter/(procs-1)+1;
                    _commPtr->send(vNtask_per_proc[i - 1], i);
                }
                for(long i=Niter%(procs-1)+1; i<procs; i++){
                    vNtask_per_proc[i - 1] = Niter/(procs-1);
                    _commPtr->send(vNtask_per_proc[i - 1], i);
                }
            }

            //Initialize the buider and the receiver vector r
            Builder_.initialize( D, Iteration(r, D) );
        }

    };

}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
