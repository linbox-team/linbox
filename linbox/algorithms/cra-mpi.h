/* Copyright (C) 2007 LinBox
 * Written by bds and zw
 *
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


#ifndef __LINBOX_cra_mpi_H
#define __LINBOX_cra_mpi_H

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
/*
template <typename T > class chooseMPItype;
template <> struct chooseMPItype<unsigned int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED;};
template <> struct chooseMPItype<unsigned long long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG_LONG;};
template <> struct chooseMPItype<unsigned long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG;};
#include <gmp++/gmp++.h>
#include <string>
*/
#include <unordered_set>
#include <linbox/randiter/random-prime.h>

namespace LinBox
{

	template<class CRABase>
	struct MPIChineseRemainder  {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
	protected:
		CRABase Builder_;
		Communicator* _commPtr;
		unsigned int _numprocs;

	public:
		template<class Param>
		MPIChineseRemainder(const Param& b, Communicator *c) :
			Builder_(b), _commPtr(c), _numprocs(c->size())
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

			int procs = _commPtr->size();
			int process = _commPtr->rank();

			//  parent process
			if(process == 0 ){
				//  create an array to store primes
				int primes[procs - 1];
				DomainElement r;
				//  send each child process a new prime to work on
				for(int i=1; i<procs; i++){
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
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
				return Builder_.result(res);
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
				return res;
			}
		}

#if 0
		template<class V, class F, class P>
		V & operator() (V& res, F& it, P&primeg){ return res; }
#endif
		template<class Vect, class Function, class PrimeIterator>
		Vect & operator() (Vect& res, Function& Iteration, PrimeIterator& primeg)
		{
			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {
				ChineseRemainder< CRABase > sequential(Builder_);
				return sequential(res, Iteration, primeg);
			}

			int procs = _commPtr->size();
			int process = _commPtr->rank();
// 			std::vector<DomainElement> r;
			typename Rebind<Vect, Domain>::other r;

			//  parent propcess
			if(process == 0){
				int primes[procs - 1];
				Domain D(*primeg);
				//  for each slave process...
				for(int i=1; i<procs; i++){
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					//  fix the array of currently sent primes
					primes[i - 1] = *primeg;
					//  send the prime to a slave process
					_commPtr->send(primes[i - 1], i);
				}
				Builder_.initialize( D, Iteration(r, D) );
				int poison_pills_left = procs - 1;
				while(poison_pills_left > 0 ){
					int idle_process = 0;
					//  receive the beginnin and end of a vector in heapspace
					_commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0);
					//  determine which process sent answer
					//  and give them a new prime
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
					Domain D(primes[idle_process - 1]);
					Builder_.progress(D, r);
					//  if still working, queue a prime
					if(! Builder_.terminated()){
						++primeg;
						primes[idle_process - 1] = *primeg;
					}
					//  otherwise, queue a poison pill
					else{
						primes[idle_process - 1] = 0;
						poison_pills_left--;
					}
					//  send the prime or poison
					_commPtr->send(primes[idle_process - 1], idle_process);
				}  // while
				return Builder_.result(res);
			}
			//  child process
			else{
				int pp;
				//  get a prime, compute, send back start and end
				//  of heap addresses
				while(true){
					_commPtr->recv(pp, 0);
					if(pp == 0)
						break;
					Domain D(pp);
					Iteration(r, D);
					_commPtr->send(r.begin(), r.end(), 0, 0);
				}
				return res;
			}
		}
	};
        
        
        
        
	template<class RatCRABase>
	struct MPIratChineseRemainder  {
		typedef typename RatCRABase::Domain	Domain;
		typedef typename RatCRABase::DomainElement	DomainElement;
	protected:
		RatCRABase Builder_;
		Communicator* _commPtr;
		unsigned int _numprocs;
                
	public:
		template<class Param>
		MPIratChineseRemainder(const Param& b, Communicator *c) :
			Builder_(b), _commPtr(c), _numprocs(c->size())
		{}

		template<class Function, class PrimeIterator>
		Integer & operator() (Integer& num, Integer& den, Function& Iteration, PrimeIterator& primeg)
		{
			//  defer to standard CRA loop if no parallel usage is desired
			if(_commPtr == 0 || _commPtr->size() == 1) {
				RationalRemainder< RatCRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}

			int procs = _commPtr->size();
			int process = _commPtr->rank();

			//  parent process
			if(process == 0 ){
				//  create an array to store primes
				int primes[procs - 1];
				DomainElement r;
				//  send each child process a new prime to work on
				for(int i=1; i<procs; i++){
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
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
						Builder_.initialize( D, Iteration(r, D) );
						first_time = false;
					}
					else
						Builder_.progress( D, Iteration(r, D) );
					//  queue a new prime if applicable
					if(! Builder_.terminated()){
						++primeg;
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
				return Builder_.result(num,den);
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
				return num;
			}
		}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ROI Optimization<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		template<class Function, class PrimeIterator>
		BlasVector<Givaro::ZRing<Integer> > & operator() ( BlasVector<Givaro::ZRing<Integer> > & num, Integer& den, Function& Iteration, PrimeIterator& primeg)
		{

#if 0       //Defuat implementation
                        
			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {
				RationalRemainder< RatCRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}
                        
			int procs = _commPtr->size();
			int process = _commPtr->rank();
			//typename Rebind<BlasVector< Givaro::ZRing<Integer> > , Domain>::other r;
                        Domain D(*primeg);
                        BlasVector<Domain> r(D);
                        //Timer chrono;
                        
			//  parent propcess
			if(process == 0){
                std::unordered_set<int> prime_sent;
				int primes[procs - 1];
				//Domain D(*primeg);
				//  for each slave process...
				for(int i=1; i<procs; i++){
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) || prime_sent.find(*primeg) != prime_sent.end()) ++primeg;	//while(Builder_.noncoprime(*primeg)) ++primeg;
					//  fix the array of currently sent primes
					primes[i - 1] = *primeg;
                    prime_sent.insert(*primeg);
					//  send the prime to a slave process
					_commPtr->send(primes[i - 1], i);
				}  

				Builder_.initialize( D, Iteration(r, D) );
				int poison_pills_left = procs - 1;
				while(poison_pills_left > 0 ){
					int idle_process = 0;
					//  receive the beginnin and end of a vector in heapspace
					_commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0); 

					//  determine which process sent answer
					//  and give them a new prime
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
					Domain D(primes[idle_process - 1]);
                                        //chrono.start(); 
					Builder_.progress(D, r);
                                        //chrono.stop(); 
                                        //std::cout<<"Builder_.progress(D, r) in the manager process used CPU time (seconds): " <<chrono.usertime()<<std::endl;
                                        
					//  if still working, queue a prime
					if(! Builder_.terminated()){
						++primeg; while(Builder_.noncoprime(*primeg) || prime_sent.find(*primeg) != prime_sent.end()) ++primeg;
                        prime_sent.insert(*primeg);
						primes[idle_process - 1] = *primeg;
					}
					//  otherwise, queue a poison pill
					else{
						primes[idle_process - 1] = 0;
						poison_pills_left--;
					}
                                        
                                        
					//  send the prime or poison
					_commPtr->send(primes[idle_process - 1], idle_process);
				}  // while
                                
				return Builder_.result(num,den);
                                
			}
			//  child process
			else{
				int pp;
				//  get a prime, compute, send back start and end
				//  of heap addresses
				while(true){
					_commPtr->recv(pp, 0);
            
					if(pp == 0)
						break;
					Domain D(pp);
                                        //chrono.start();
					Iteration(r, D);
                                        //chrono.stop(); 
                                        //std::cout << "Iteration(r,D) in the worker process used CPU time (seconds):  " << chrono.usertime() << std::endl;
std::cout << pp << std::endl;

					_commPtr->send(r.begin(), r.end(), 0, 0); 
				}
			}
                        
#endif
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#if 0       //Using news prime number generation function to reduce MPI communication between manager and workers

			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {
				RationalRemainder< RatCRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}
                        
			int procs = _commPtr->size();
			int process = _commPtr->rank();
			//typename Rebind<BlasVector< Givaro::ZRing<Integer> > , Domain>::other r;
                        Domain D(*primeg);
                        BlasVector<Domain> r(D);
                        Timer chrono;
                        
			//  parent propcess
			if(process == 0){
                //std::unordered_set<int> prime_sent;
				int primes[procs - 1];
				//Domain D(*primeg);
				//  for each slave process...
				for(int i=1; i<procs; i++){
					//  generate a new prime
					//++primeg; while(Builder_.noncoprime(*primeg) /*|| prime_sent.find(*primeg) != prime_sent.end()*/) ++primeg;	//while(Builder_.noncoprime(*primeg)) ++primeg;
					//  fix the array of currently sent primes
					primes[i - 1] = 1;
                                      //  prime_sent.insert(*primeg);
					//  send the prime to a slave process
					_commPtr->send(primes[i - 1], i);
				}  

				Builder_.initialize( D, Iteration(r, D) );
				int poison_pills_left = procs - 1;
int pp;
				while(poison_pills_left > 0 ){
					int idle_process = 0;
r.resize (r.size()+1);
					//  receive the beginnin and end of a vector in heapspace
					_commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0); 

					//  determine which process sent answer
					//  and give them a new prime
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
if(primes[idle_process - 1]==0)  poison_pills_left--;
					//  send the prime or poison
					_commPtr->send(primes[idle_process - 1], idle_process);

//std::cout << "<<<< received r:  " << r << std::endl;
pp = r[r.size()-1];
r.resize (r.size()-1); 
//std::cout << "<<<< resized r:  " << r << std::endl;



if(!Builder_.noncoprime(pp)){

					Domain D(pp); //Domain D(primes[idle_process - 1]);

//std::cerr << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "<< std::endl;
					Builder_.progress(D, r);
//std::cerr << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< "<< std::endl;
}//END FOR :  if(Builder_.noncoprime(pp))

					if(Builder_.terminated()){
						primes[idle_process - 1] = 0;
						//poison_pills_left--;
					}                                      
                                        

				}  // while
                                
				return Builder_.result(num,den);
                                
			}
			//  child process
			else{
				int pp;
        LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>   gen(process,procs);  
				//  get a prime, compute, send back start and end
				//  of heap addresses
				while(true){
					_commPtr->recv(pp, 0);
//std::cout << "Proc("<<process<<") received tag: "<< pp << std::endl;
					if(pp == 0)
						break;
++gen; while(Builder_.noncoprime(*gen) ) ++gen;
std::cout << *gen << std::endl;
					Domain D(*gen); //Domain D(pp);                                       
					Iteration(r, D);

r.push_back(*gen);
//std::cout <<"Proc("<<process<<")  Sending r:  " << r <<" for prime "<< *gen << std::endl;
					_commPtr->send(r.begin(), r.end(), 0, 0); 
				}
			}
   
#endif
/////////////////////////////////////////BEGIN OF ROI////////////////////////////////////////////////////
#if 1       //Restruture using Bcast and Gather (>_<) still not work as it never terminates

			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {
				RationalRemainder< RatCRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}
                        
			int procs = _commPtr->size();
			int process = _commPtr->rank();
            Domain D(*primeg);
            BlasVector<Domain> r(D); BlasVector<Domain> rr(D);
int tag=1;
Builder_.initialize( D, Iteration(r, D) );
rr.resize ((num.size()+1)*(procs));
r.resize (num.size());

while(true){

MPI_Bcast(&tag, 1, MPI_INT, 0, MPI_COMM_WORLD);
if(tag==0) break;

			//  child process
			if(process != 0){ 

        LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>   gen(process,procs);  

++gen; while(Builder_.noncoprime(*gen) ) ++gen;
//std::cout << *gen << std::endl;
					Domain D(*gen); //Domain D(pp);                                       
					Iteration(r, D);

std::cout <<"Proc("<<process<<")  Sending r:  " << r <<" for prime "<< *gen << std::endl;
r.push_back(*gen);

			}


//if(process != 0) std::cout <<"Proc("<<process<<")  Sending r:  " << r << std::endl;

MPI_Gather(&r.front(), (num.size()+1), MPI_DOUBLE, &rr.front(), (num.size()+1), MPI_DOUBLE, 0,  MPI_COMM_WORLD);

//std::cout <<"Proc("<<process<<")  received rr is of size :"<<rr.size()<<std::endl; 



			//  parent propcess
			if(process == 0){

int pp;

					int idle_process = 0;
/*
std::cerr <<"Received rr : "<< std::endl;
 for(long iter=0; iter<(num.size()+1)*(procs-1)-1; iter++)  std::cerr <<rr[iter]<<"   "; 
std::cerr << std::endl;
*/
std::cerr <<"Proc("<<process<<")  received rr: "<<rr<<std::endl; 

    for(long iter=1; iter<procs; iter++){
std::cerr <<"Proc("<<process<<") >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "<< std::endl;

        for(long i=0; i<num.size()+1; i++) r[i] = rr[iter*(num.size()+1)+i];
        pp = r[num.size()];

        std::cerr << "<<<< received r:  " << r <<" for prime: "<<pp<< std::endl;

        r.resize (num.size()); 

        //std::cout << "<<<< resized r:  " << r << std::endl;

        if(!Builder_.noncoprime(pp)){
					        Domain D(pp); //Domain D(primes[idle_process - 1]);
					        Builder_.progress(D, r);
std::cerr <<"<><><><> Checked r:  " << r <<" for prime "<< pp <<" <><><><>"<< std::endl;
        }//END :  if(Builder_.noncoprime(pp))
std::cerr <<"Proc("<<process<<") <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< "<< std::endl;
//MPI_Recv(&r, 1, MPI_INT,  0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    if(Builder_.terminated()){
              tag = 0;                             
              break;
        }

    }//END for loop


    }//END : if(process == 0)

r.resize (num.size()+1);

   
}//END while(true)
return Builder_.result(num,den);
#endif

//////////////////////////////////////////END OF ROI//////////////////////////////////////////////////////

		}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END ROI Optimization>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
