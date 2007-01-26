/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/* author: B. David Saunders and Zhendong Wan*/
// parallelized for BOINC computing by Bryan Youse
// ======================================================================= //
// Time-stamp: <15 Jul 05 18:53:10 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_CRA_H
#define __LINBOX_CRA_H

#define MPICH_IGNORE_CXX_SEEK
#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>
//#include "backend_lib.h"
#include "linbox/util/mpicpp.h"

//#include <linbox/util/timer.h>

namespace LinBox {


        // Lazy computation of the product of the moduli
    struct LazyProduct : public std::vector< Integer > {
        typedef std::vector< Integer > Father_t;
    protected:
        bool                _tobecomputed;
    public:
        LazyProduct() : Father_t(), _tobecomputed(true) {}
        void initialize(const Integer& i) {
            _tobecomputed = false;
            this->resize(0);
            this->push_back(i);
        }
            
        bool mulin(const Integer& i) {
            if (this->size()) {
                if (i != this->back()) {
                    this->push_back( i );
                    return _tobecomputed = true;
                } else {
                    return _tobecomputed;
                }
            
            } else {
                this->push_back( i );
                return _tobecomputed = false;
            }
        }
      
        bool mulin(const LazyProduct& i) {
            this->insert(this->end(), i.begin(), i.end());
            return _tobecomputed = (this->size()>1);
        }
      
        Integer & operator() () {
            if (_tobecomputed) {
                Father_t::const_iterator iter = this->begin();
                Father_t::iterator       prod = this->begin();
                for(++iter; iter != this->end(); ++iter)
                    *prod *= *iter;
                this->resize(1);
                _tobecomputed = false;
            }
            return this->back();
        }

        bool noncoprime(const Integer& i) const {
            Integer g;
            for(Father_t::const_iterator iter = this->begin(); iter != this->end(); ++iter)
                if ( gcd(g,i,*iter) > 1) return true;
            return false;
        }   
       
        friend std::ostream& operator<< (std::ostream& o, const LazyProduct& C) {
            o << "{";
            for(Father_t::const_iterator refs = C.begin();
                refs != C.end() ;
                ++refs )
                o << (*refs) << " " ;
            return o << "}";
        }

    };

/* Warning, we won't detect bad primes */

    template<class Domain>
    struct ChineseRemainder {
    protected:
        typedef typename Domain::Element DomainElement;
        typedef unsigned long            BaseRing; // used anywhere?
      

        unsigned int               occurency;
        double                     dSizes0; // used anywhere?
        Integer             nextm;
        Integer               Modulo0;
        Integer               Table0;
        bool                       Occupation0; //used anywhere?
        Communicator               *commPtr;
        unsigned int               numprocs;
        std::vector< unsigned long >      randv;
        std::vector< double >           dSizes;
        std::vector< LazyProduct >      Modulo;
        std::vector< std::vector<Integer> >    Table;
        std::vector< bool >             Occupation;
        unsigned int             EARLY_TERM_THRESHOLD;
        double                LOGARITHMIC_UPPER_BOUND;
        unsigned int             step;

    public:

        ChineseRemainder(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD, const size_t n=1) : Modulo0(0) {
            initialize(EARLY, n, 0.0);
        }

		  ChineseRemainder(Communicator *c = NULL, const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD,
		                   const size_t n=1) : Modulo0(0) {
            commPtr = c;
            numprocs = commPtr->size();
            initialize(EARLY, n, 0.0);
		  }

        ChineseRemainder(const double BOUND) : Modulo0(0) {
            initialize(0, 0, BOUND);
        }

        virtual ~ChineseRemainder() {}
        

        void initialize(const unsigned long EARLY=1, const size_t n=1, const double BOUND=0.0) {
            step = 0;
            occurency = 0;
            Occupation0 = false;
            dSizes.resize(1); Modulo.resize(1); Table.resize(1); Occupation.resize(1);
            Modulo0 = 1;
            EARLY_TERM_THRESHOLD = EARLY;
            LOGARITHMIC_UPPER_BOUND = BOUND;
            if (EARLY) {
                srand48(BaseTimer::seed());
                std::vector<unsigned long>::iterator int_p;
                randv. resize (n);
                for (int_p = randv. begin(); 
                     int_p != randv. end(); ++ int_p) 
                    *int_p = ((unsigned long)lrand48()) % 20000;
            }
        }   

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
        template<class Function, class RandPrime>
        Integer & operator() (Integer& res, Function& Iteration, RandPrime& genprime, Communicator *Comm) {
		  		//  defer to standard CRA loop if no parallel usage is desired
            if(Comm == 0 || Comm->size() == 1)
					return this->operator()(res, Iteration, genprime);

            Integer p;

            int procs = Comm->size();
            int process = Comm->rank();

            if (EARLY_TERM_THRESHOLD) {
				    //  parent process
                if(process == 0 ){
					     //  create an array to store primes
				        int primes[procs - 1];
                    genprime.randomPrime(p);
                    DomainElement r; 
						  //  send each child process a new prime to work on
					     for(int i=1; i<procs; i++){
                     	genprime.randomPrime(p);
                    	   while(this->Early_noncoprimality(p) )
                    	     genprime.randomPrime(p);
								primes[i - 1] = p;
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
									this->First_Early_progress(D, r);
									first_time = false;
								}
								else
									this->Early_progress( D, r );
								//  queue a new prime if applicable
								if(! this->Early_terminated()){
								   genprime.randomPrime(p);
									primes[idle_process - 1] = p;
								}
								//  otherwise, queue a poison pill
								else{
									primes[idle_process - 1] = 0;
									poison_pills_left--;
								}
								//  send the prime or poison pill
								Comm->send(primes[idle_process - 1], idle_process);
						 	}  // end while
                  return this->Early_result(res);
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
            } else {
				    if(process == 0)
                while( ! this->Full_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Full_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    DomainElement r; D.init(r);
                    this->Full_progress( D, Iteration(r, D) );
                }
                return this->Full_result(res);
            }
        }


            /** \brief The CRA loop
            
            Given a function to generate residues mod a single prime, this loop produces the residues 
            resulting from the Chinese remainder process on sufficiently many primes to meet the 
            termination condition.
         
            \param F - Function object of two arguments, F(r, p), given prime p it outputs residue(s) r.
            This loop may be parallelized.  F must be reentrant, thread safe.
            For example, F may be returning the coefficients of the minimal polynomial of a matrix mod p.
            Warning - we won't detect bad primes.
         
            \param genprime - RandIter object for generating primes.
            \result res - an integer
            */
        template<class Function, class RandPrime>
        Integer & operator() (Integer& res, Function& Iteration, RandPrime& genprime) {
            Integer p;
            if (EARLY_TERM_THRESHOLD) {
                {
                    genprime.randomPrime(p);
                    Domain D(p); 
                    DomainElement r; D.init(r);
                    this->First_Early_progress( D, Iteration(r, D) );
                }
                while( ! this->Early_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Early_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    DomainElement r; D.init(r);
                    this->Early_progress( D, Iteration(r, D) );
                }
                return this->Early_result(res);
            } else {
                while( ! this->Full_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Full_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    DomainElement r; D.init(r);
                    this->Full_progress( D, Iteration(r, D) );
                }
                return this->Full_result(res);
            }
        }

        template<template <class T> class Vect, class Function, class RandPrime>
        Vect<Integer> & operator() (Vect<Integer>& res, Function& Iteration, RandPrime& genprime, Communicator *Comm) {
				//  if there is no communicator or if there is only one process,
				//  then proceed normally (without parallel)
		  		if(Comm == 0 || Comm->size() == 1)
					return this->operator()(res, Iteration, genprime);

            Integer p;
				int procs = Comm->size();
				int process = Comm->rank();
				Vect<DomainElement> r;

            if (EARLY_TERM_THRESHOLD) {
					  //  parent propcess
					  if(process == 0){
				        int primes[procs - 1];
					     genprime.randomPrime(p);
						  Domain D(p);
						  //  for each slave process...
						  for(int i=1; i<procs; i++){
						     //  generate a new prime
						     genprime.randomPrime(p);
							  while(this->Early_noncoprimality(p) )
							     genprime.randomPrime(p);
							  //  fix the array of currently sent primes
							  primes[i - 1] = p;	
							  //  send the prime to a slave process
							  Comm->send(primes[i - 1], i);
						  }
						  this->First_Early_progress( D, Iteration(r, D) ); 
						  int idle_process = 0;
						  int poison_pills_left = procs - 1;
						  while(poison_pills_left > 0 ){
						     //  receive the beginnin and end of a vector in heapspace
						     Comm->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0);
							  //  determine which process sent answer 
							  //  and give them a new prime
							  idle_process = (Comm->get_stat()).MPI_SOURCE;
							  Domain D(primes[idle_process - 1]);
							  this->Early_progress(D, r);
							  //  if still working, queue a prime
							  if(! this->Early_terminated()){
							     genprime.randomPrime(p);
								  primes[idle_process - 1] = p;
							  }
							  //  otherwise, queue a poison pill
							  else{
							     primes[idle_process - 1] = 0;
								  poison_pills_left--;
							  }
							  //  send the prime or poison
							  Comm->send(primes[idle_process - 1], idle_process);
						  }  // while
					     return this->Early_result(res);	
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
            } else {
				    if(process == 0)
                while( ! this->Full_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Full_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<DomainElement> r; 
                    this->Full_progress( D, Iteration(r, D) );
                }
                return this->Full_result(res);
            }
			}


        template<template <class T> class Vect, class Function, class RandPrime>
        Vect<Integer> & operator() (Vect<Integer>& res, Function& Iteration, RandPrime& genprime) {
            Integer p;
            if (EARLY_TERM_THRESHOLD) {
                {
                    genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<DomainElement> r; 
                    this->First_Early_progress( D, Iteration(r, D) );            
                }
                while( ! this->Early_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Early_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<DomainElement> r; 
                    this->Early_progress( D, Iteration(r, D) );
                }
                return this->Early_result(res);
            } else {
                while( ! this->Full_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Full_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<DomainElement> r; 
                    this->Full_progress( D, Iteration(r, D) );
                }
                return this->Full_result(res);
            }
        }

       
            /** \brief Function for adding a new prime and it's residue to the CRA process.
                \param D - A domain (e.g. a Field).  
                \param e - A residue, image in the domain D of the desired value.
            */
        template<class VectOrInt>
        void progress (const Domain& D, const VectOrInt& e) {
            if( ++step == 1) {
                if (EARLY_TERM_THRESHOLD)
                    First_Early_progress(D, e);
                else
                    Full_progress(D, e);
            } else {
                if (EARLY_TERM_THRESHOLD)
                    Early_progress(D, e);
                else
                    Full_progress(D, e);
            }
        }

   int steps() {return step;}        
      
            /** \brief result mod the lcm of the moduli.
            
            A value mod the lcm of the progress step moduli
            which agrees with each residue mod the corresponding modulus. 
            */
        template<class VectOrInt>
        VectOrInt& result (VectOrInt &d){
            return (EARLY_TERM_THRESHOLD? Early_result(d) : Full_result(d) );
        }
 
        bool terminated() { 
            if (EARLY_TERM_THRESHOLD) 
                return Early_terminated();
            else
                return Full_terminated();
        }


        bool noncoprime(const Integer& i) const { 
            if (EARLY_TERM_THRESHOLD) 
                return Early_noncoprimality(i);
            else
                return Full_noncoprimality(i);
        }

    protected:
      
        bool Early_noncoprimality(const Integer& i) const {
            Integer g;
            return (gcd(g, i, Modulo0) != 1) ;
        }
      
            
        bool Full_noncoprimality(const Integer& i) const {
            std::vector< LazyProduct >::const_iterator _mod_it = Modulo.begin();
            std::vector< bool >::const_iterator    _occ_it = Occupation.begin();
            for( ; _occ_it != Occupation.end() ; ++_mod_it, ++_occ_it)
                if ((*_occ_it) && (_mod_it->noncoprime(i))) return true;
            return false;
        }
      
        void Full_progress (const Domain& D, const DomainElement& e) {
            std::vector<DomainElement> Ve; Ve.push_back(e);
            Full_progress(D, Ve);
        }

            
        template<template<class T> class Vect>
        void Full_progress (const Domain& D, const Vect<DomainElement>& e) {         
            std::vector< double >::iterator  _dsz_it = dSizes.begin();
            std::vector< LazyProduct >::iterator _mod_it = Modulo.begin();
            std::vector< std::vector<Integer> >::iterator _tab_it = Table.begin();
            std::vector< bool >::iterator    _occ_it = Occupation.begin();
            std::vector<Integer> ri(e.size()); LazyProduct mi; double di;
            if (*_occ_it) {
                typename Vect<DomainElement>::const_iterator  e_it = e.begin();
                std::vector<Integer>::iterator       ri_it = ri.begin();
                std::vector<Integer>::const_iterator t0_it = _tab_it->begin();
                for( ; e_it != e.end(); ++e_it, ++ri_it, ++ t0_it)
                    fieldreconstruct(*ri_it, D, *e_it, *t0_it, (*_mod_it).operator()() );
                Integer tmp; D.characteristic(tmp);
                di = *_dsz_it + log(double(tmp)); 
                mi.mulin(tmp);
                mi.mulin(*_mod_it);
                *_occ_it = false;
            } else {
                Integer tmp; D.characteristic(tmp);
                _mod_it->initialize(tmp);
                *_dsz_it = log(double(tmp));
                typename Vect<DomainElement>::const_iterator e_it = e.begin();
                _tab_it->resize(e.size());
                std::vector<Integer>::iterator t0_it= _tab_it->begin();
                for( ; e_it != e.end(); ++e_it, ++ t0_it)
                    D.convert(*t0_it, *e_it);
                *_occ_it = true;
                return;
            }
         
            for(++_dsz_it, ++_mod_it, ++_tab_it, ++_occ_it ; _occ_it != Occupation.end() ; ++_dsz_it, ++_mod_it, ++_tab_it, ++_occ_it) {
                if (*_occ_it) {
                    std::vector<Integer>::iterator      ri_it = ri.begin();
                    std::vector<Integer>::const_iterator t_it= _tab_it->begin();
                    for( ; ri_it != ri.end(); ++ri_it, ++ t_it)
                        smallbigreconstruct(*ri_it, mi(), *t_it, _mod_it->operator()()); 
                    mi.mulin(*_mod_it);
                    di += *_dsz_it;
                    *_occ_it = false;
                } else {
                    *_dsz_it = di;
                    *_mod_it = mi;
                    *_tab_it = ri;
                    *_occ_it = true;
                    return;
                }
            }
                // Everything was occupied
            dSizes.push_back( di );
            Table.push_back( ri );
            Modulo.push_back( mi );
            Occupation.push_back ( true );
        }

        Integer& Full_result (Integer &d){
            std::vector<Integer> Vd; Vd.push_back(d);
            Full_result(Vd);
            return d=Vd.front();
        }
      
        template<template<class T> class Vect>
        Vect<Integer>& Full_result (Vect<Integer> &d){
            d.resize( (Table.front()).size() );
            std::vector< LazyProduct >::iterator          _mod_it = Modulo.begin();
            std::vector< std::vector< Integer > >::iterator _tab_it = Table.begin();
            std::vector< bool >::iterator                _occ_it = Occupation.begin();
            LazyProduct Product;
            for( ; _occ_it != Occupation.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
                if (*_occ_it) {
                    Product = *_mod_it;
                    std::vector<Integer>::iterator t0_it = d.begin();
                    std::vector<Integer>::iterator t_it = _tab_it->begin();
                    if (++_occ_it == Occupation.end()) {
                        for( ; t0_it != d.end(); ++t0_it, ++t_it)
                            normalize(*t0_it = *t_it, *t_it, _mod_it->operator()());
                        Modulo.resize(1);
                        return d;
                    } else {
                        for( ; t0_it != d.end(); ++t0_it, ++t_it)
                            *t0_it  = *t_it;
                        ++_mod_it; ++_tab_it; 
                        break;
                    }
                }
            }
            for( ; _occ_it != Occupation.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
                if (*_occ_it) {
                    std::vector<Integer>::iterator t0_it = d.begin();
                    std::vector<Integer>::const_iterator t_it = _tab_it->begin();
                    for( ; t0_it != d.end(); ++t0_it, ++t_it)
                        normalizesmallbigreconstruct(*t0_it, Product(), *t_it, _mod_it->operator()() );
                    Product.mulin(*_mod_it);
                }
            }
            Modulo.resize(1);
            Modulo.front() = Product;
            return d;
        }
      
        bool Full_terminated() {
            double logm(0.0);
            std::vector< double >::iterator _dsz_it = dSizes.begin();
            std::vector< bool >::iterator _occ_it = Occupation.begin();
            for( ; _occ_it != Occupation.end() ; ++_dsz_it, ++_occ_it) {
                if (*_occ_it) {
                    logm = *_dsz_it;
                    ++_dsz_it; ++_occ_it;
                    break;
                }
            }
            for( ; _occ_it != Occupation.end() ; ++_dsz_it, ++_occ_it) {
                if (*_occ_it) {
                    logm += *_dsz_it;
                }
            }
            return logm > LOGARITHMIC_UPPER_BOUND;
        }


        template<template<class T> class Vect>
        void Early_progress (const Domain& D, const Vect<DomainElement>& e) {
            DomainElement z;
                // Could be much faster
                // - do not compute twice the product of moduli
                // - reconstruct one element of e until Early Termination,
                //   then only, try a random linear combination.
            Early_progress(D, dot(z, D, e, randv));
            Full_progress(D, e);
        }

        template<template<class T> class Vect>
        void First_Early_progress (const Domain& D, const Vect<DomainElement>& e) {
            DomainElement z;
                // Could be much faster
                // - do not compute twice the product of moduli
                // - reconstruct one element of e until Early Termination,
                //   then only, try a random linear combination.
            First_Early_progress(D,dot(z, D, e, randv) );
            Full_progress(D, e);
        }

        virtual void Early_progress (const Domain& D, const DomainElement& e) {
				if(Modulo0 == 0){
               D.characteristic( Modulo0 );
               D.convert(Table0, e);
               nextm = 1;
				}
				else{
               Modulo0 *= nextm; D.characteristic( nextm );
               Early_normalized_fieldreconstruct(Table0, Modulo0, D, e, nextm);
				}
        }

        virtual void First_Early_progress (const Domain& D, const DomainElement& e) {
            D.characteristic( Modulo0 );
            D.convert(Table0, e);
            nextm = 1;
        }

        template<template<class T> class Vect>
        Vect<Integer>& Early_result(Vect<Integer>& d) {
            return Full_result(d);
        }      


        Integer& Early_result(Integer& d) {
            return d=Table0;
        }      

        bool Early_terminated() {
            return occurency>EARLY_TERM_THRESHOLD;
        }

        Integer& fieldreconstruct(Integer& res, const Domain& D1, const DomainElement& u1, const Integer& r0, const Integer& P0) {
            DomainElement u0, m0; D1.init(u0, r0);
            if (D1.areEqual(u1, u0))
                return res=r0;
            else
                return fieldreconstruct(res, D1, u1, u0, D1.init(m0, P0), r0, P0);
        }

        Integer& Early_normalized_fieldreconstruct(Integer& res, const Integer& P0, const Domain& D1, const DomainElement& u1, const Integer& M1) {
            DomainElement u0, m0; D1.init(u0, res);

//             D1. write( 
//                 D1. write( 
//                     D1.write(commentator.report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) 
//                              << " over ") 
//                     << " : " << res << " mod " << P0 << " is ", u0) 
//                 << " || ", u1) << " mod " << M1 << std::endl;
                

            if (D1.areEqual(u1, u0)) {
                ++occurency;
                return res;
            } else {
                occurency = 0;
                return normalized_fieldreconstruct(res, D1, u1, u0, D1.init(m0, P0), Integer(res), P0, M1);
            }
         
        }

        Integer& fieldreconstruct(Integer& res, const Domain& D1, const DomainElement& u1, DomainElement& u0, DomainElement& m0, const Integer& r0, const Integer& P0) {
                // u0 and m0 are modified
            D1.negin(u0);   // u0 <-- -u0
            D1.addin(u0,u1);   // u0 <-- u1-u0
            D1.invin(m0);   // m0 <-- m0^{-1} mod m1
            D1.mulin(u0, m0);   // u0 <-- (u1-u0)( m0^{-1} mod m1 )
            D1.convert(res, u0);// res <-- (u1-u0)( m0^{-1} mod m1 )         and res <  m1
            res *= P0;      // res <-- (u1-u0)( m0^{-1} mod m1 ) m0      and res <= (m0m1-m0)
            return res += r0;   // res <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0 and res <  m0m1
        }

        Integer& normalized_fieldreconstruct(Integer& res, const Domain& D1, const DomainElement& u1, DomainElement& u0, DomainElement& m0, const Integer& r0, const Integer& P0, const Integer& m1) {
                // u0 and m0 are modified
            D1.negin(u0);          // u0 <-- -u0
            D1.addin(u0,u1);       // u0 <-- u1-u0
            D1.invin(m0);          // m0 <-- m0^{-1} mod m1
            D1.mulin(u0, m0);      // u0 <-- (u1-u0)( m0^{-1} mod m1 )
            D1.convert(res, u0);   // res <-- (u1-u0)( m0^{-1} mod m1 )       and res <  m1
            Integer tmp(res);
            tmp -= m1;
            if (absCompare(res,tmp)>0) res = tmp;
            res *= P0;      // res <-- (u1-u0)( m0^{-1} mod m1 ) m0       and res <= (m0m1-m0)
            return res += r0;   // res <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0  and res <  m0m1
        }

      
        Integer& smallbigreconstruct(Integer& u1, const Integer& m1, const Integer& u0, const Integer& m0) {
            Integer tmp=u1;
            inv(u1, m0, m1);      // res <-- m0^{-1} mod m1
            tmp -= u0;            // tmp <-- (u1-u0)
            u1 *= tmp;            // res <-- (u1-u0)( m0^{-1} mod m1 )   
            u1 %= m1;      //                    res <  m1 
            u1 *= m0;             // res <-- (u1-u0)( m0^{-1} mod m1 ) m0      and res <= (m0m1-m0)
            return u1 += u0;      // res <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0 and res <  m0m1
        }

        Integer& normalize(Integer& u1, Integer& tmp, const Integer& m1) {
            if (u1 < 0)
                tmp += m1;
            else
                tmp -= m1;
            if (absCompare(u1,tmp) > 0) return u1 = tmp;
            else          return u1;
        }
        


        Integer& normalizesmallbigreconstruct(Integer& u1, const Integer& m1, const Integer& u0, const Integer& m0) {
            Integer tmp=u1;
            inv(u1, m0, m1);   // res <-- m0^{-1} mod m1
            tmp -= u0;      // tmp <-- (u1-u0)
            u1 *= tmp;      // res <-- (u1-u0)( m0^{-1} mod m1 )   
            u1 %= m1;
            normalize(u1, tmp=u1, m1);         // Normalization
            u1 *= m0;          // res <-- (u1-u0)( m0^{-1} mod m1 ) m0       and res <= (m0m1-m0)
            return u1 += u0;   // res <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0  and res <  m0m1
        }


   template <template<class T> class Vect1, class Vect2>
   DomainElement& dot (DomainElement& z, const Domain& D, const Vect1<DomainElement>& v1, const Vect2& v2){
            D.init(z,0); DomainElement tmp;
            typename Vect1<DomainElement>::const_iterator v1_p;
            typename Vect2::const_iterator v2_p;
            for (v1_p  = v1. begin(), v2_p = v2. begin(); 
                 v1_p != v1. end(); 
                 ++ v1_p, ++ v2_p)       
                D.axpyin(z, (*v1_p), D.init(tmp, (*v2_p)));

//             commentator.report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "v: " << v2 << std::endl;
//             commentator.report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "z: " << z << std::endl;
            return z;
   }

    };
}

#endif
