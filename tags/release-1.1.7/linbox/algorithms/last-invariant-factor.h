/* Copyright (C) LinBox
 *
 *  Author: Zhendong Wan
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

#ifndef __LINBOX_last_invariant_factor_H
#define __LINBOX_last_invariant_factor_H

#include <linbox/util/debug.h>
#include <linbox/algorithms/default.h>
#include <linbox/algorithms/rational-solver.h>
#include <utility>

namespace LinBox {
	
/** \brief This is used in a Smith Form algorithm.

This computes the last invariant factor of an integer matrix,
whether zero or not, by rational solving.
*/
	template<class _Ring,
		class _Solver>
		
		class LastInvariantFactor {
			
		public:
		
			typedef _Ring Ring;	
			typedef _Solver Solver;
			typedef typename Ring::Element Integer;
			
		protected:
			
			Ring r;
			Solver solver;
			int threshold;

		public:	

	/** _Ring, an integer ring,
	 *  _Solver, a function which solves Ax = b over the quotient field of _Ring.
         */
			LastInvariantFactor(const Solver& _solver = Solver(),
					    const Ring& _r =Ring(),
					    int _threshold =DEFAULTLIFTHRESHOLD)
				:r(_r),solver(_solver), threshold(_threshold) {
				
				if ( _threshold <= 1) threshold = DEFAULTLIFTHRESHOLD;
			}

			void setThreshold (int _threshold) {
				if (_threshold > 1) {
					threshold = _threshold; 
				}
			}
			
			int getThreshold () const {
				return threshold;
			}
			
			const Solver& getSolver() const {
				return solver;
			}

			void setSolver(const Solver& s) {
				solver = s;
			}
			
			/** \brief Compute the last invariant factor of an integer matrix,
			 * by solving linear system,
			 * ignoring these factors of primes in list PrimeL
			 */
			template<class IMatrix, class Vector>
				Integer& lastInvariantFactor(Integer& lif, const IMatrix& A, 
							     const Vector& PrimeL) const{
			
				r.init(lif, 1);
				int count = 0;
				SolverReturnStatus tmp;
				// Storage of rational solution
				std::vector<Integer> r_num (A. coldim()); Integer r_den;
				//std::vector<std::pair<Integer, Integer> > result (A.coldim());
				//typename std::vector<std::pair<Integer, Integer> >::iterator result_p;
				// vector b, RHS, 32-bit int is good enough
				std::vector<int> b(A.rowdim());
				typename std::vector<int>::iterator b_p;
				typename Vector::const_iterator Prime_p;

				Integer pri, quo, rem;
				
				for (; count < threshold; ++ count) {
					// assign b to be a random vector
					for (b_p = b.begin(); b_p != b.end(); ++ b_p) {
						* b_p = rand() % 268435456 - 134217728; // may need to change to use ring's random gen.
						                // dpritcha, 2004-07-26
					}
					
					// try to solve Ax = b over Ring
					tmp = solver.solveNonsingular(r_num, r_den, A, b);
					// If no solution found
					if (tmp != SS_OK) {
						r.init (lif, 0);
						break;
					}
					
					r. lcmin (lif, r_den);
				}
				// filter out primes in PRIMEL from lif.					
				if (!r. isZero (lif))
					for ( Prime_p = PrimeL.begin(); 
					      Prime_p != PrimeL.end();
					      ++ Prime_p) {
						r.init (pri, *Prime_p);
						do {
							r.quoRem(quo,rem,lif,pri);
							if (r.isZero(rem)) r.assign(lif,quo);
							else break;
						}
						while (true);
					}
				
				return lif;
			}
			
			/** \brief Compute the last invariant factor of an integer matrix,
			 * by solving linear system,
			 * ignoring these factors of primes in list PrimeL
			 * Implement the bonus in ref{....}
			 */
			template<class IMatrix, class Vector>
				Integer& lastInvariantFactor_Bonus(Integer& lif, Integer& bonus, const IMatrix& A, 
							     const Vector& PrimeL) const{
			
				r. init(lif, 1);
				r. init (bonus, 1);
				int count = 0;
				SolverReturnStatus tmp1, tmp2;
				// Storage of rational solution
				std::vector<Integer> r1_num (A. coldim()), r2_num (A. coldim()); Integer r1_den, r2_den;
				//std::vector<std::pair<Integer, Integer> > result (A.coldim());
				//typename std::vector<std::pair<Integer, Integer> >::iterator result_p;
				// vector b, RHS, 32-bit int is good enough
				std::vector<int> b1(A. rowdim()), b2(A. rowdim());
				typename std::vector<int>::iterator b_p;
				typename Vector::const_iterator Prime_p;
				Integer pri, quo, rem;
				
				for (; count < (threshold + 1) / 2; ++ count) {
					// assign b to be a random vector
					for (b_p = b1. begin(); b_p != b1. end(); ++ b_p) {
						* b_p = rand();
					}
					for (b_p = b2. begin(); b_p != b2. end(); ++ b_p) {
						* b_p = rand();
					}
					// try to solve Ax = b1, b2 over Ring
					tmp1 = solver. solveNonsingular(r1_num, r1_den, A, b1);
					tmp2 = solver. solveNonsingular(r2_num, r2_den, A, b2);
					// If no solution found
					if ((tmp1 != SS_OK) || (tmp2 != SS_OK)){
						r.init (lif, 0);
						break;
					}

					r. lcm (lif, lif, r1_den);
					r. lcm (lif, lif, r2_den);

					// compute the bonus
					Integer g, d, a11, a12, a21, a22, l, c_bonus, c_l;
					typename std::vector<Integer>::iterator num1_p, num2_p;
					std::vector<Integer> r1 (A. rowdim());
					std::vector<Integer> r2 (A. rowdim());
					typename std::vector<Integer>::iterator r1_p, r2_p;
					r. init (l, 0);
					int i;
					for (i = 0; i < 20; ++ i) {
						for (r1_p = r1. begin(), r2_p = r2. begin(); r1_p != r1. end(); ++ r1_p, ++ r2_p) {
							r. init (*r1_p, rand());
							r. init (*r2_p, rand());
						}
						r. init (a11, 0); r. init (a12, 0); r. init (a21, 0); r. init (a22, 0);
						for (r1_p = r1. begin(), num1_p = r1_num. begin(); r1_p != r1. end(); ++ r1_p, ++ num1_p)
							r. axpyin (a11, *r1_p, *num1_p);
						for (r1_p = r1. begin(), num2_p = r2_num. begin(); r1_p != r1. end(); ++ r1_p, ++ num2_p)
							r. axpyin (a12, *r1_p, *num2_p);
						for (r2_p = r2. begin(), num1_p = r1_num. begin(); r2_p != r2. end(); ++ r2_p, ++ num1_p)
							r. axpyin (a21, *r2_p, *num1_p);
						for (r2_p = r2. begin(), num2_p = r2_num. begin(); r2_p != r2. end(); ++ r2_p, ++ num2_p)
							r. axpyin (a22, *r2_p, *num2_p);
						// g = gcd (a11, a12, a21, a22)
						r. gcd (g, a11, a12); r. gcdin (g, a21); r. gcdin (g, a22);
						// d = a11 a22 - a12 a21
						r. mul (d, a12, a21); r. negin (d); r. axpyin (d, a11, a22);
						if (! r. isZero (g)) r. div (c_l, d, g);
						r. gcdin (l, c_l);
					}

					if (!r. isZero (l) ) {
						r. gcd (c_bonus, r1_den, r2_den);
						r. gcdin (l, c_bonus);
						r. divin (c_bonus, l);
					}

					r. lcmin (bonus, c_bonus);
				}

				// filter out primes in PRIMEL from lif.					
				if (!r. isZero (lif)) 
					for ( Prime_p = PrimeL.begin(); Prime_p != PrimeL.end(); ++ Prime_p) {
						r.init (pri, *Prime_p);
						do {
							r.quoRem(quo,rem,lif,pri);
							if (r.isZero(rem)) r.assign(lif,quo);
							else break;
						} while (true);
					} 
				r. gcdin (bonus, lif);
				if (!r. isZero (bonus))
					for ( Prime_p = PrimeL.begin(); Prime_p != PrimeL.end(); ++ Prime_p) {
						r.init (pri, *Prime_p);
						do {
							r.quoRem(quo,rem,bonus,pri);
							if (r.isZero(rem)) r.assign(lif,quo);
							else break;
						} while (true);
					} 

				
				return lif;
			}

                        template<class IMatrix, class Vector>
                        Integer& lastInvariantFactor1(Integer& lif, Vector& r_num, const IMatrix& A, const bool oldMatrix=false) const {
                                //cout << "enetering lif\n";
                                SolverReturnStatus tmp;

                                if (r_num.size()!=A. coldim()) return lif=0;
                                Integer r_den;
                                std::vector<Integer> b(A.rowdim());
                                typename std::vector<Integer>::iterator b_p;
                                //typename Vector::const_iterator Prime_p;

                                Integer pri, quo, rem;

                                // assign b to be a random vector
                                for (b_p = b.begin(); b_p != b.end(); ++ b_p) {
                                        * b_p = rand() % 268435456 - 134217728; // may need to change to use ring's random gen.
                                        // dpritcha, 2004-07-26
                                }
                                //report <<"try to solve Ax = b over Ring";
                                // try to solve Ax = b over Ring
                                //cout << "trying to solve\n";
                                tmp = solver.solveNonsingular(r_num, r_den, A, b,oldMatrix);
                                // If no solution found
                                if (tmp != SS_OK) {
                                        //r.init (lif, 0);
                                        //break;
                                        return lif=0;
                                }
                                //report << "r_den: "<< r_den;

                                r. lcmin (lif, r_den);
				if (r_den != lif) {
                                	Integer den,t;
                                	r. lcm(den,r_den,lif);
                                	r. div(t, den, r_den);
                                	typename std::vector<Integer>::iterator num_p = r_num.begin();
                                	for (; num_p != r_num. end(); ++num_p) {
	                                	r. mulin(*num_p, t);
					}	
				}
                                return lif;
                        }

                        template<class Vector>
                        Integer& bonus(Integer& bonus, const Integer r1_den,const Integer r2_den, Vector& r1_num, Vector& r2_num) const {
                                if (bonus==0) bonus=1;
                                if (r1_num.size() != r2_num.size()) return bonus=0;
                                Integer g, d, a11, a12, a21, a22, c_bonus, l, c_l;
                                typename std::vector<Integer>::iterator num1_p, num2_p;
                                std::vector<Integer> r1 (r1_num. size());
                                std::vector<Integer> r2 (r2_num. size());
                                typename std::vector<Integer>::iterator r1_p, r2_p;
                                r. init (l, 0);
                                int i;
                                for (i = 0; i < 20; ++ i) {
                                        for (r1_p = r1. begin(), r2_p = r2. begin(); r1_p != r1. end(); ++ r1_p, ++ r2_p) {
                                                r. init (*r1_p, rand());
                                                r. init (*r2_p, rand());
                                        }
                                        r. init (a11, 0); r. init (a12, 0); r. init (a21, 0); r. init (a22, 0);
                                        for (r1_p = r1. begin(), num1_p = r1_num. begin(); r1_p != r1. end(); ++ r1_p, ++ num1_p)
                                                r. axpyin (a11, *r1_p, *num1_p);
                                        for (r1_p = r1. begin(), num2_p = r2_num. begin(); r1_p != r1. end(); ++ r1_p, ++ num2_p)
                                                r. axpyin (a12, *r1_p, *num2_p);
                                        for (r2_p = r2. begin(), num1_p = r1_num. begin(); r2_p != r2. end(); ++ r2_p, ++ num1_p)
                                                r. axpyin (a21, *r2_p, *num1_p);
                                        for (r2_p = r2. begin(), num2_p = r2_num. begin(); r2_p != r2. end(); ++ r2_p, ++ num2_p)
                                                r. axpyin (a22, *r2_p, *num2_p);
                                        // g = gcd (a11, a12, a21, a22)
                                        r. gcd (g, a11, a12); r. gcdin (g, a21); r. gcdin (g, a22);
                                        // d = a11 a22 - a12 a21
                                        r. mul (d, a12, a21); r. negin (d); r. axpyin (d, a11, a22);
                                        if (! r. isZero (g)) r. div (c_l, d, g);
                                        r. gcdin (l, c_l);
                                }
                                if (!r. isZero (l) ) {
                                        r. gcd (c_bonus, r1_den, r2_den);
                                        r. gcdin (l, c_bonus);
                                        r. divin (c_bonus, l);
                                }

                                r. lcmin (bonus, c_bonus);
                                return bonus;
                        }

			
			/** \brief Compute the last invariant factor.
			 */
			template<class IMatrix>
			  Integer& lastInvariantFactor(Integer& lif, const IMatrix& A)  const {

				std::vector<Integer> empty_v;
				lastInvariantFactor (lif, A, empty_v);
				return lif;
			}

			/** \brief Compute the last invariant factor with bonus
			 */
			template<class IMatrix>
			  Integer& lastInvariantFactor_Bonus(Integer& lif, Integer& bonus, const IMatrix& A)  const {

				std::vector<Integer> empty_v;
				lastInvariantFactor_Bonus (lif, bonus, A, empty_v);
				return lif;
			}
	
	};
}


#endif //__LINBOX_last_invariant_factor_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
