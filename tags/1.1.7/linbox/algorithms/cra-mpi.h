/* Copyright (C) 2007 LinBox
 * Written by bds and zw
 *
 * author: B. David Saunders and Zhendong Wan
 * parallelized for BOINC computing by Bryan Youse
 *
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
#include "linbox/util/mpicpp.h"

namespace LinBox 
{

	template<class CRABase>
		struct MPIChineseRemainder {
			typedef typename CRABase::Domain	Domain;
			typedef typename CRABase::DomainElement	DomainElement;
			protected:
			CRABase Builder_;
			Communicator* _commPtr;
			unsigned int _numprocs;

			public:
			template<class Param>
				MPIChineseRemainder(const Param& b, Communicator *c) : Builder_(b), _commPtr(c), _numprocs(c->size()) {}

			/** \brief The CRA loop

			  termination condition.

			  \param F - Function object of two arguments, F(r, p), given prime p it outputs residue(s) r.
			  This loop may be parallelized.  F must be reentrant, thread safe.
			  For example, F may be returning the coefficients of the minimal polynomial of a matrix mod p.
			  Warning - we won't detect bad primes.

			  \param genprime - RandIter object for generating primes.
			  \param Comm - Pointer to Communicator object to delegate parallelism using MPI
			  \result res - an integer
			  */
			template<class Function, class PrimeIterator>
				Integer & operator() (Integer& res, Function& Iteration, PrimeIterator& primeg) {
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
						int idle_process = 0;
						bool first_time = true;
						int poison_pills_left = procs - 1;
						//  loop until all execution is complete
						while( poison_pills_left > 0 ){
							//  receive sub-answers from child procs
							_commPtr->recv(r, MPI_ANY_SOURCE);
							idle_process = (_commPtr->get_stat()).MPI_SOURCE;
							Domain D(primes[idle_process - 1]); 
							//  assimilate results
							if(first_time){
								Builder_.initialize(D, r);
								first_time = false;
							} else
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

			/*
			template<class V, class F, class P>
				V & operator() (V& res, F& it, P&primeg){ return res; }
				*/

			template<class Vect, class Function, class PrimeIterator>
				Vect & operator() (Vect& res, Function& Iteration, PrimeIterator& primeg) {
					//  if there is no communicator or if there is only one process,
					//  then proceed normally (without parallel)
					if(_commPtr == 0 || _commPtr->size() == 1) {
						ChineseRemainder< CRABase > sequential(Builder_);
						return sequential(res, Iteration, primeg);
					}

					int procs = _commPtr->size();
					int process = _commPtr->rank();
					std::vector<DomainElement> r;

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
						int idle_process = 0;
						int poison_pills_left = procs - 1;
						while(poison_pills_left > 0 ){
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
}

#undef MPICH_IGNORE_CXX_SEEK
#endif // __LINBOX_cra_mpi_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
