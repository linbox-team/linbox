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

#ifndef __LINBOX_hybrid_det_H
#define __LINBOX_hybrid_det_H

//#include "linbox/blackbox/diagonal.h"
//#include "linbox/blackbox/compose.h"
//#include "linbox/solutions/methods.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/sparse.h"
//#include "linbox/blackbox/blas-blackbox.h"
//#include "linbox/matrix/blas-matrix.h"
//#include "linbox/algorithms/blackbox-container.h"
//#include "linbox/algorithms/blackbox-container-symmetric.h"
//#include "linbox/algorithms/massey-domain.h"
//#include "linbox/algorithms/blas-domain.h"
//#include "linbox/vector/vector-traits.h"
//#include "linbox/util/prime-stream.h"
//#include "linbox/util/debug.h"

//#include "linbox/solutions/solve.h"
//#include "linbox/field/gmp-rational.h"
//#include "linbox/field/gmp-integers.h"
#include "linbox/field/PID-integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/algorithms/last-invariant-factor.h"
//#include "linbox/field/PIR-modular-int32.h"
//#include "linbox/field/PIR-ntl-ZZ_p.h"
//#include "linbox/field/ntl-ZZ_p.h"

// Namespace in which all LinBox library code resides

//#include "linbox/algorithms/cra.h"
//#include "linbox/field/modular.h"
#include "linbox/field/modular-double.h"
//#include "linbox/field/givaro-zpz.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

#include "linbox/solutions/det.h"

namespace LinBox 
{
   
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



                //          int iter_count;
                //          int iter_count2;
                //          Vector <PID_integer>:: Dense moduli;
                //          Vector <PID_integer>:: Dense primes;

            IntegerModularDetReduced(const Blackbox& b, const MyMethod& n, const integer& divisor, const size_t& f)
                    : A(b), M(n), beta(divisor), factor(f) {
                moduli.resize(factor);
                primes.resize(factor);
                iter_count = 0;
                iter_count2 = 0;

            }

            template<typename Field>
            typename Field::Element& operator()(typename Field::Element& d, const Field& F) {

                if (beta > 1) {
                    if (iter_count2 < factor) {
                        Field D(primes[iter_count2]);
                        typename Field::Element pbeta;
                        typename Field::Element kbeta;
                        typename Field::Element current_moduli;
                        D.init(pbeta, beta);
                        D.init(current_moduli, moduli[iter_count2]);
                        D.div(kbeta, current_moduli,pbeta);
                        ++this->iter_count2;
                        return d=kbeta;
                    }
                }

                typedef typename Blackbox::template rebind<Field>::other FBlackbox;
                FBlackbox Ap(A,F);
                detin( d, Ap, M);

                if (beta > 1) {
                    typename Field::Element y;
                    F.init(y,beta);
                    F.div(d,d,y);
                }

                if (iter_count < factor) {
                    moduli[iter_count] = d;
                }
                ++this->iter_count;

                return d;
            }

            void Beta(Integer& b) { beta = b; iter_count2=0;}

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
                //commentator.setReportStream(std::cout);
                typedef Modular<double> myModular;
                typedef typename Blackbox::Field Integers;
                typedef typename Integers::Element Integer;

                commentator.start ("Integer Determinant - hybrid version ", "det");
                size_t myfactor=5;
                size_t early_counter=0;

            	//double a = 0.0;//0.28;//0.0013//if time(lif)/time(10lu) < a * log(lif) then calculate bonus

                Integer lif = 1;
                Integer bonus = 1;
                Integer beta = 1;
                d=1;
		
		double p_size = 26-(int)ceil(log((double)A.rowdim())*0.7213475205);

                RandomPrimeIterator genprime( (Integer)p_size );
		//cout << "prime size: " << p_size << "\n";
                EarlySingleCRA<myModular> cra(4UL);
                IntegerModularDetReduced<Blackbox,MyMethod> iteration(A, M, beta,myfactor);

                //if (A.rowdim() < 200 ) {
                //    cra(d,iteration,genprime);
                //    commentator.stop ( "first step", NULL, "det");
                //    commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //        << "Iterations done " << iteration.iterations() << "\n";
                //} else { 
                Integer res;

		Timer BT;
		double time1, time2;
		BT.start();
		
                while ( early_counter < myfactor && !cra.terminated() ) {
			++genprime;
			while (cra.noncoprime(*genprime)) ++genprime;
                        myModular D(*genprime);
                        iteration.primes[early_counter] = *genprime;
                            //		prime(p, early_counter);
                        myModular::Element r;
                        D.init(r,0);
                        cra.progress( D, iteration(r, D));
                        ++early_counter; 
                }
		
		BT.stop();
		time1 = BT.usertime()/myfactor;
		if (time1 < 0) time1 = 0;
                cra.result(res);

                if (early_counter < myfactor) {
                    /* determinant found */
			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << myfactor << "\n";
			commentator.stop ( "first step", NULL, "det");
                        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                            << "Iterations done " << iteration.iterations()<< "\n";
                        return d=res;
                }
		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << "no very early termination \n";
		/* turn to LU when matrix size small - not to be used at the moment */ 
		//if (A.rowdim() < 50 ) {
                //while ( !cra.terminated() ) {
		//	genprime.randomPrime(p);
                //        myModular D(p);
                //        iteration.primes[early_counter] = p;
                            //		prime(p, early_counter);
                //        myModular::Element r;
                //        D.init(r,0);
                //        cra.progress( D, iteration(r, D));
                //        ++early_counter; 
                //}
		//commentator.stop ( "zero step", NULL, "det");
		//commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //        << "Iterations done " << iteration.iterations() << "\n";
		//return d;
		//}
		
                //RandomPrime genprime1( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));
                //Integers ZZ;
                //RationalSolver < Integers , myModular, RandomPrime, DixonTraits > RSolver(A. field(), genprime); 
		RationalSolver < Integers , myModular, RandomPrimeIterator, DixonTraits > RSolver;
		
		typename Vector<Integers>:: Dense r_num1 (A. coldim());
		
                LastInvariantFactor < Integers ,RationalSolver < Integers, myModular, RandomPrimeIterator, DixonTraits > >  LIF(RSolver);

		BT.start();
		if (LIF.lastInvariantFactor1(lif, r_num1, A)==0) {
			//if (lif==0) {
			d = 0;
			commentator.stop ("is 0", NULL, "det");
			return d;
		}
		BT.stop();
		time2 = BT.usertime();
		if (time2 < 0) time2 =0;

		//commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //            << "5 LU time: " << time1 << " LIF time: " << time2 << ")\n";
		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << "lif calculated\n";
                //LIF.lastInvariantFactor_Bonus(lif, bonus, A);	
	
                //if (lif==0) {
		//	d = 0;
		//      commentator.stop ("done", NULL, "det");
                //        return d;
		//}
 
                //if (bonus == 1) {
		//	d = lif;
                //        commentator.stop ("done", NULL, "det");
                //        return d;
                //}
                        
                beta = lif*bonus;
                iteration.Beta(beta);

                //RandomPrime genprime2( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
                EarlySingleCRA< Modular<double> > cra2(4UL);
                Integer k = 1;

                early_counter = 0;
                while ( early_counter < myfactor && !cra2.terminated() ) {
                        myModular D(iteration.primes[early_counter]);
                        myModular::Element r;
                        D.init(r,0);
                        cra2.progress( D, iteration(r, D) );
                        ++early_counter; 
                }

                if (early_counter < myfactor) {
                            /* determinant found */
			k = cra2.result(res);
			commentator.stop ("second step ", NULL, "det");
			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                            << "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
			
                }
                else  if (0/* time2 < a*log2(lif)*time1/p_size*/) {
			typename Vector<Integers>:: Dense r_num2 (A. coldim());
			Integer lif2=1;
			LIF.lastInvariantFactor1(lif2,r_num2,A);
			LIF.bonus(bonus,lif, lif2, r_num1,r_num2);
			
			//lif2 = lif;
			if ((bonus > 1) || (lif2 != lif)) {
				//cout << "lif: " <<lif <<",\n     "<< lif2<< "przed\n";
				lif = lcm(lif, lif2);
				//cout << "lif: " <<lif << "po";
				beta=lif*bonus;
				iteration.Beta(beta);
				//iteration.Restart2() // included in Beta();
				//iteration.Moduli(moduli);
				//iteration.Primes(primes);
				k=1;
			        EarlySingleCRA< Modular<double> > cra3(4UL);

				early_counter = 0;
				while ( (early_counter < myfactor) && (!cra3.terminated() )) {
					myModular D(iteration.primes[early_counter]);
					myModular::Element r;
					D.init(r,0);
					cra3.progress( D, iteration(r, D));
					++early_counter;
					//iteration.Inc();
				}

				if (early_counter < myfactor) {
					/* determinant found based on the initial LU */
					cra3.result(k);
					commentator.stop ("third step - recalc", NULL, "det");
					commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
						<< "Iterations done " << iteration.iterations()<<"(" << iteration.iterations2() << ")\n";
						//<< "bonus size " << log2(bonus) << "\n";
					
				} else {
					/* enter the cra loop */
					//cra3(k,iteration, genprime);
					while (!cra3.terminated()) {
						++genprime;
						while (cra3.noncoprime(*genprime)) ++genprime;
						myModular D(*genprime);
						myModular::Element r;
						D.init(r,0);
						cra3.progress( D, iteration(r, D));			
					}
					cra3.result(k);
				       	commentator.stop ("third step, bonus > 1", NULL, "det");
					commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
						<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
                                                //<< "bonus size " << log2(bonus) << "\n";
				}
			} else {
				//cra2(k,iteration, genprime);
				while (!cra2.terminated()) {
					++genprime;
					while (cra2.noncoprime(*genprime)) ++genprime;
					myModular D(*genprime);
					myModular::Element r;
					D.init(r,0);
					cra2.progress( D, iteration(r, D));
				}
				cra2.result(k);
				commentator.stop ("third step, bonus = 1", NULL, "det");
				commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
					<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";		
			}
		} else {
				while (!cra2.terminated()) {
					++genprime;
					while (cra2.noncoprime(*genprime)) ++genprime;
					myModular D(*genprime);
					myModular::Element r;
					D.init(r,0);
					cra2.progress( D, iteration(r, D));
				}
				cra2.result(k);
				commentator.stop ("second step+", NULL, "det");
				commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
					<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";

                           /* enter the cra loop */
                        //cra2(k,iteration, genprime);
                }

                //commentator.stop ("second step", NULL, "det");
                //commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //        << "Iterations done " << iteration.iterations() << "(" 
                //        << iteration.iterations2() << " )\n";
		d = k*beta;
		
		Integer tmp;
		
                commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                            << "1 LU time: " << time1 << " LIF time: " << time2 << ")\n";
		commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                        << "det/lif " << k<< "\n";
		//commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		//	<< "determinant size " << Integers::log2(tmp,abs(d))<< " lif size "<< Integers::log2(tmp,beta) << "\n";

		return d ;
	
	    }
#if 0
        template <class Integers, class MyMethod>
        typename Integers::Element & lif_cra_det (typename Integers::Element         &d,
                                                        const SparseMatrix<Integers >                            &A,
                                                        const RingCategories::IntegerTag          &tag,
                                                        const MyMethod                            &M)
            {

                //commentator.setReportStream(std::cout);
                typedef Modular<double> myModular;
                //typedef PID_integer Integers;
                typedef typename Integers::Element Integer;

                commentator.start ("Integer Determinant - hybrid version for sparse matrices", "det");
                size_t myfactor=5;
                size_t early_counter=0;

                //double a = 0.0;//0.28;//0.0013i			//if time(lif)/time(10lu) < a * log(lif) then calculate bonus

                Integer lif = 1;
                Integer bonus = 1;
                Integer beta = 1;
                d=1;

                double p_size = 26-(int)ceil(log((double)A.rowdim())*0.7213475205);

                RandomPrime genprime( (Integer)p_size );
                commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << "prime size: " << p_size << "\n"; 
                ChineseRemainder< myModular > cra(3UL);
                IntegerModularDetReduced<SparseMatrix<Integers >,MyMethod> iteration(A, M, beta,myfactor);

                //if (A.rowdim() < 200 ) {
                //    cra(d,iteration,genprime);
                //    commentator.stop ( "first step", NULL, "det");
                //    commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //        << "Iterations done " << iteration.iterations() << "\n";
                //} else {
                Integer p;
                Integer res;

                Timer BT;
                double time1, time2;
                BT.start();

                while ( early_counter < myfactor && !cra.terminated() ) {
                        genprime.randomPrime(p);
                        while (cra.noncoprime(p)) genprime.randomPrime(p);
                        myModular D(p);
                        iteration.primes[early_counter] = p;
                            //          prime(p, early_counter);
                        myModular::Element r;
                        D.init(r,0);
                        cra.progress( D, iteration(r, D));
                        ++early_counter;
                }

                BT.stop();
                time1 = BT.usertime()/myfactor;
                if (time1 < 0) time1 = 0;
                cra.result(res);

                if (early_counter < myfactor) {
                     //determinant found 
                        commentator.stop ( "first step", NULL, "det");
                        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                            << "Iterations done " << iteration.iterations()<< "\n";
                        return d=res;
                }
                //cout << "no very early termination \n";
                // turn to LU when matrix size small - not to be used at the moment 
                //if (A.rowdim() < 50 ) {
                //while ( !cra.terminated() ) {
                //      genprime.randomPrime(p);
                //        myModular D(p);
                //        iteration.primes[early_counter] = p;
                            //          prime(p, early_counter);
                //        myModular::Element r;
                //        D.init(r,0);
                //        cra.progress( D, iteration(r, D));
                //        ++early_counter;
                //}
                //commentator.stop ( "zero step", NULL, "det");
                //commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //        << "Iterations done " << iteration.iterations() << "\n";
                //return d;
                //}

                //RandomPrime genprime1( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));
                //Integers ZZ;
                typedef RationalSolver < Integers , myModular, RandomPrime, BlockHankelTraits > Solver;
		Solver RSolver(A. field(), genprime);

                typename Vector<Integers>:: Dense r_num1 (A. coldim());

                LastInvariantFactor < Integers ,Solver >  LIF(RSolver);

                BT.start();
                if (LIF.lastInvariantFactor1(lif, r_num1, A)==0) {
                        //if (lif==0) {
                        d = 0;
                        commentator.stop ("is 0", NULL, "det");
                        return d;
                }
                BT.stop();
                time2 = BT.usertime();
                if (time2 < 0) time2 =0;

                //commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //            << "5 LU time: " << time1 << " LIF time: " << time2 << ")\n";
                //cout << "lif calculated\n";
                //LIF.lastInvariantFactor_Bonus(lif, bonus, A);

                //if (lif==0) {
                //      d = 0;
                //      commentator.stop ("done", NULL, "det");
                //        return d;
                //}

                //if (bonus == 1) {
                //      d = lif;
                //        commentator.stop ("done", NULL, "det");
                //        return d;
                //}

                beta = lif*bonus;
                iteration.Beta(beta);

                //RandomPrime genprime2( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));
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
                            // determinant found 
                        k = cra2.result(res);
                        commentator.stop ("second step ", NULL, "det");
                        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                            << "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";

                }
                else  if (0) { //time2 < a*log2(lif)*time1/p_size) {
                        typename Vector<Integers>:: Dense r_num2 (A. coldim());
                        Integer lif2=1;
                        LIF.lastInvariantFactor1(lif2,r_num2,A);
                        LIF.bonus(bonus,lif, lif2, r_num1,r_num2);

                        //lif2 = lif;
                        if ((bonus > 1) || (lif2 != lif)) {
                                //cout << "lif: " <<lif <<",\n     "<< lif2<< "przed\n";
                                lif = lcm(lif, lif2);
                                //cout << "lif: " <<lif << "po";
                                beta=lif*bonus;
                                iteration.Beta(beta);
                                //iteration.Restart2() // included in Beta();
                                //iteration.Moduli(moduli);
                                //iteration.Primes(primes);
                                k=1;
                                ChineseRemainder< Modular<double> > cra3(3UL);

                                early_counter = 0;
                                while ( (early_counter < myfactor) && (!cra3.terminated() )) {
                                        myModular D(iteration.primes[early_counter]);
                                        myModular::Element r;
                                        D.init(r,0);
                                        cra3.progress( D, iteration(r, D));
                                        ++early_counter;
                                        //iteration.Inc();
                                }

                                if (early_counter < myfactor) {
                                        // determinant found based on the initial LU 
                                        cra3.result(k);
                                        commentator.stop ("third step - recalc", NULL, "det");
                                        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                                                << "Iterations done " << iteration.iterations()<<"(" << iteration.iterations2() << ")\n";
                                                //<< "bonus size " << log2(bonus) << "\n";

                                } else {
                                        // enter the cra loop 
                                        //cra3(k,iteration, genprime);
                                        while (!cra3.terminated()) {
                                                genprime.randomPrime(p);
                                                while (cra3.noncoprime(p)) genprime.randomPrime(p);
                                                myModular D(p);
                                                myModular::Element r;
                                                D.init(r,0);
                                                cra3.progress( D, iteration(r, D));
                                        }
                                        cra3.result(k);
                                        commentator.stop ("third step, bonus > 1", NULL, "det");
                                        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                                                << "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
                                                //<< "bonus size " << log2(bonus) << "\n";
                                }
                        } else {
                                //cra2(k,iteration, genprime);
                                while (!cra2.terminated()) {
                                        genprime.randomPrime(p);
                                        while (cra2.noncoprime(p)) genprime.randomPrime(p);
                                        myModular D(p);
                                        myModular::Element r;
                                        D.init(r,0);
                                        cra2.progress( D, iteration(r, D));
                                }
                                cra2.result(k);
                                commentator.stop ("third step, bonus = 1", NULL, "det");
                                commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                                        << "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
                        }
                } else {
                                while (!cra2.terminated()) {
                                        genprime.randomPrime(p);
                                        while (cra2.noncoprime(p)) genprime.randomPrime(p);
                                        myModular D(p);
                                        myModular::Element r;
                                        D.init(r,0);
                                        cra2.progress( D, iteration(r, D));
                                }
                                cra2.result(k);
                                commentator.stop ("second step+", NULL, "det");
                                commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                                        << "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";

                           // enter the cra loop 
                        //cra2(k,iteration, genprime);
                }

                //commentator.stop ("second step", NULL, "det");
                //commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //        << "Iterations done " << iteration.iterations() << "("
                //        << iteration.iterations2() << " )\n";
                d = k*beta;

                Integer tmp;

                commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                            << "1 LU time: " << time1 << " LIF time: " << time2 << ")\n";
                commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                        << "det/lif " << k<< "\n";
                //commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                //        << "determinant size " << Integers::log2(tmp,abs(d))<< " lif size "<< Integers::log2(tmp,beta) << "\n";

                return d ;

            }
#endif

} // end of LinBox namespace
#endif // __LINBOX_hybrid_det_H




    
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
