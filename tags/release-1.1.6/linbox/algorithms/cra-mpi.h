/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/* author: B. David Saunders and Zhendong Wan*/
// parallelized for BOINC computing by Bryan Youse
// ======================================================================= //
// Time-stamp: <15 Mar 07 17:41:24 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_CRA_MPI_H
#define __LINBOX_CRA_MPI_H

#define MPICH_IGNORE_CXX_SEEK
#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>
#include "linbox/algorithms/cra-domain.h"
#include "linbox/util/mpicpp.h"

namespace LinBox {
    
    template<class CRABase>
    struct MPIChineseRemainder {
        typedef typename CRABase::Domain	Domain;
        typedef typename CRABase::DomainElement	DomainElement;
    protected:
        CRABase Builder_;
        Communicator               *commPtr;
        unsigned int               numprocs;
        
    public:
        template<class Param>
        MPIChineseRemainder(const Param& b, Communicator *c = NULL) : Builder_(b), commPtr(c), numprocs(c->size()) {}

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
        Integer & operator() (Integer& res, Function& Iteration, PrimeIterator& primeg, Communicator *Comm) {
                //  defer to standard CRA loop if no parallel usage is desired
            if(Comm == 0 || Comm->size() == 1) {
                ChineseRemainder< CRABase > sequential(Builder_);
                return sequential(res, Iteration, primeg);
            }
            
            int procs = Comm->size();
            int process = Comm->rank();
            
                //  parent process
            if(process == 0 ){
                    //  create an array to store primes
                int primes[procs - 1];
                DomainElement r; 
                    //  send each child process a new prime to work on
                for(int i=1; i<procs; i++){
                    ++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
                    primes[i - 1] = *primeg;
                    Comm->send(primes[i - 1], i);
                }
                int idle_process = 0;
                bool first_time = true;
                int poison_pills_left = procs - 1;
                    //  loop until all execution is complete
                while( poison_pills_left > 0 ){
                        //  receive sub-answers from child procs
                    Comm->recv(r, MPI_ANY_SOURCE);
                    idle_process = (Comm->get_stat()).MPI_SOURCE;
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
                    Comm->send(primes[idle_process - 1], idle_process);
                }  // end while
                return Builder_.result(res);
            }  // end if(parent process)
                //  child processes
            else{
                int pp;
                while(true){ 
                        //  receive the prime to work on, stop
                        //  if signaled a zero
                    Comm->recv(pp, 0);
                    if(pp == 0)
                        break;
                    Domain D(pp);
                    DomainElement r; D.init(r);
                    Iteration(r, D);
                        //Comm->buffer_attach(rr);
                        // send the results
                    Comm->send(r, 0);
                }
                return res;
            }
        }

        
        template<template <class T> class Vect, class Function, class PrimeIterator>
        Vect<Integer> & operator() (Vect<Integer>& res, Function& Iteration, PrimeIterator& primeg, Communicator *Comm) {
                //  if there is no communicator or if there is only one process,
                //  then proceed normally (without parallel)
            if(Comm == 0 || Comm->size() == 1) {
                ChineseRemainder< CRABase > sequential(Builder_);
                return sequential(res, Iteration, primeg);
            }
            
            int procs = Comm->size();
            int process = Comm->rank();
            Vect<DomainElement> r;
            
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
                    Comm->send(primes[i - 1], i);
                }
                Builder_.initialize( D, Iteration(r, D) ); 
                int idle_process = 0;
                int poison_pills_left = procs - 1;
                while(poison_pills_left > 0 ){
                        //  receive the beginnin and end of a vector in heapspace
                    Comm->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0);
                        //  determine which process sent answer 
                        //  and give them a new prime
                    idle_process = (Comm->get_stat()).MPI_SOURCE;
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
                    Comm->send(primes[idle_process - 1], idle_process);
                }  // while
                return Builder_.result(res);	
            }
                //  child process
            else{
                int pp;
                    //  get a prime, compute, send back start and end
                    //  of heap addresses 
                while(true){
                    Comm->recv(pp, 0);
                    if(pp == 0)
                        break;
                    Domain D(pp);
                    Iteration(r, D);
                    Comm->send(r.begin(), r.end(), 0, 0);
                }
                return res;
            }
        }

    };

}

#endif
