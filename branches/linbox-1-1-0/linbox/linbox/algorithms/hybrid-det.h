/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/hybrid-det.h
 * Copyright (C) 2005 Anna Urbanska
 *
 * Written by Anna Urbanska <aniau@astronet.pl>
 * Modified by JGD <Jean-Guillaume.Dumas@imag.fr>
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
#ifndef __HYBRID_DET_H
#define __HYBRID_DET_H

#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/solutions/methods.h"
#include "linbox/blackbox/dense.h"

#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/prime-stream.h"
#include "linbox/util/debug.h"

//#include "linbox/solutions/solve.h"
#include "linbox/field/gmp-rational.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/field/PID-integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/field/PIR-modular-int32.h"
#include "linbox/field/ntl-ZZ_p.h"
#include "linbox/field/PIR-ntl-ZZ_p.h"
// Namespace in which all LinBox library code resides

//#include "linbox/algorithms/cra.h"
#include "linbox/field/modular.h"
//#include "linbox/field/givaro-zpz.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

#include "linbox/solutions/det.h"



namespace LinBox {

       

         
        template <class Blackbox, class MyMethod>
        struct IntegerModularDetReduced {      
        private:
            
            const Blackbox &A;
            const MyMethod &M;
                /* contains the factor by which to divide det */
            integer beta;

            size_t factor;
            size_t iter_count;
            size_t iter_count2;
            typename Vector<PID_integer>::Dense moduli;
            
        public:

            typename Vector<PID_integer>::Dense primes;

            size_t iterations() {
                return iter_count;
            }
            size_t iterations2() {
                return iter_count2;
            }
            
                    
            
                //	    int iter_count;
                //	    int iter_count2;
                //	    Vector <PID_integer>:: Dense moduli;
                //	    Vector <PID_integer>:: Dense primes;
            
            IntegerModularDetReduced(const Blackbox& b, const MyMethod& n, const integer& divisor, const size_t& f)
                    : A(b), M(n), beta(divisor), factor(f) {                
                moduli.resize(factor);
                primes.resize(factor);

            }
            
            template<typename Field>
            typename Field::Element& operator()(typename Field::Element& d, const Field& F) {
                
                if (beta > 1) {
                    ++this->iter_count2;
                    if (iter_count2 < factor) {
                        Field D(primes[iter_count2]);
                        typename Field::Element pbeta;
                        typename Field::Element kbeta;
                        typename Field::Element current_moduli;
                        D.init(pbeta, beta);
                        D.init(current_moduli, moduli[iter_count2]);
                        D.div(kbeta, current_moduli,pbeta);
                        return d=kbeta;
                    }
                }
                
                typedef typename Blackbox::template rebind<Field>::other FBlackbox;
                FBlackbox * Ap;
                MatrixHom::map(Ap, A, F);
                det( d, *Ap, M);
                
                if (beta > 1) {
                    typename Field::Element y;
                    F.init(y,beta);
                    F.div(d,d,y);
                }

                delete Ap;
	
                if (iter_count < factor) {
                    moduli[iter_count] = d;
                }
                ++this->iter_count;

                return d;
            }

            void Beta(Integer& b) { beta = b; }

        };
	
	/** \brief Compute the determinant of A over the integers
	 *
	 * The determinant of a linear operator A, represented as a
	 * black box, is computed over the integers.
         * 
         * This variant is a hybrid between Last Invariant Factor and Chinese Remaindering
         * It performs several determinants mod primes
         * Then switches to the LIF method, producing a factor of det.
         * It then comes back to the CRA if necessary to compute
         * the remaining (usually small) factor of the determinant.
	 *
	 * @param d Field element into which to store the result
	 * @param A Black box of which to compute the determinant
         * @param tag explicit over the integers
	 * @param M may be a Method::BlasElimination (default) or a Method::Wiedemann.
         \ingroup solutions
        */
        template <class Blackbox, class MyMethod>
        typename Blackbox::Field::Element & lif_cra_det (typename Blackbox::Field::Element         &d,
                                                        const Blackbox                            &A,
                                                        const RingCategories::IntegerTag          &tag,
                                                        const MyMethod                            &M)
            {
                typedef Modular<double> myModular;
                commentator.start ("Integer Determinant - hybrid version ", "det");
                size_t myfactor=6;
                size_t early_counter = 0;

                Integer lif = 1;
                Integer bonus = 1;
                Integer beta = 1;
                d=1;

                RandomPrime genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));
                ChineseRemainder< myModular > cra(3UL);
                IntegerModularDetReduced<Blackbox,MyMethod> iteration(A, M, beta,myfactor);

                if (A.rowdim() < 50 ) {
                    cra(d,iteration,genprime);
                    commentator.stop ( "first step", NULL, "det");
                    commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                        << "Iterations done " << iteration.iterations() << "\n";
                } else { 
                    Integer p;
                    Integer res;
                    while ( early_counter < myfactor && !cra.terminated() ) {
                        genprime.randomPrime(p);
                        myModular D(p);
                        iteration.primes[early_counter] = p;
                            //		prime(p, early_counter);
                        myModular::Element r;
                        D.init(r,0);
                        cra.progress( D, iteration(r, D));
                        ++early_counter; 
                    }

                    cra.result(res);

                    if (early_counter < myfactor) {
                            /* determinant found */
                        commentator.stop ( "first step", NULL, "det");
                        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                            << "Iterations done " << iteration.iterations()<< "\n";
                        return d=res;
                    }
		
                    RandomPrime genprime1( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));
                    PID_integer ZZ;
                    RationalSolver < PID_integer , myModular, RandomPrime, DixonTraits > RSolver(ZZ, genprime1); 

                    LastInvariantFactor < PID_integer ,RationalSolver < PID_integer, Modular<double>, RandomPrime, DixonTraits > >  LIF(RSolver);
                    LIF.lastInvariantFactor_Bonus(lif, bonus, A);	
	
                    if (lif==0) {
                        d = 0;
                        commentator.stop ("done", NULL, "det");
                        return d;
                    }

                    if (bonus == 1) {
                        d = lif;
                        commentator.stop ("done", NULL, "det");
                        return d;
                    }
                        
                    beta = lif*bonus;
                    iteration.Beta(beta);

                    RandomPrime genprime2( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
                    ChineseRemainder< Modular<double> > cra2(3UL);
                    Integer k = 1;

                    early_counter = 0;
                    while ( early_counter < myfactor && !cra2.terminated() ) {
                        myModular D(iteration.primes[early_counter]);
                        myModular::Element r;
                        D.init(r,0);
                        cra2.progress( D, iteration(r, D));
                        ++early_counter; 
                    }

                    if (early_counter < myfactor) {
                            /* determinant found */
                        k = cra2.result(res);
                    }
                    else {
                            /* enter the cra loop */
                        cra2(k,iteration, genprime);
                    }

                    commentator.stop ("second step", NULL, "det");
                    commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                        << "Iterations done " << iteration.iterations() << "(" 
                        << iteration.iterations2() << " )\n";
                    d = k*beta;

                }

                return d ;
	
            }
    
} // end of LinBox namespace
#endif // __HYBRID_DET_H




    
