/** File: last-invariant-factor.h
 *  Author: Zhendong Wan
 */

#ifndef __LINBOX__LAST_INVARIANT_FACTOR_H__
#define __LINBOX__LAST_INVARIANT_FACTOR_H__

#include <linbox/util/debug.h>
#include <linbox/algorithms/default.h>
#include <utility>

namespace LinBox {
	
	/** _Ring, an integer ring,
	 *  _Solver, soler Ax = b over the quotient field of _Ring.
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
			
			
			/** @memo Compute the last invariant factor of an integer matrix,
			 * by solving linear system,
			 * ignoring these factors of primes in list PrimeL
			 */
			template<class IMatrix, class Vector>
				Integer& lastInvariantFactor(Integer& lif, const IMatrix& A, 
							     const Vector& PrimeL)  const {
				
			
				r.init(lif, 1);
				
				int count = 0;

				int tmp;
				
				std::vector<std::pair<Integer, Integer> > result (A.coldim());
				
				typename std::vector<std::pair<Integer, Integer> >::iterator result_p;
						       
				// vector b
				std::vector<int> b(A.rowdim());
				
				std::vector<int>::iterator b_p;

				typename Vector::const_iterator Prime_p;

				Integer pri, quo, rem;
				
				for (; count < threshold; ++ count) {
					
					// assign b to be a random vector
					for (b_p = b.begin(); b_p != b.end(); ++ b_p) {
						* b_p = rand();
					}
					
					// try to solve Ax = b over Ring
					tmp = solver.solve(result, A, b, false);
					
					// If no solution found
					if (tmp) {
						r.init (lif, 0);
						break;
					}
					
					//prev = lif;

					// compute lcm (den[0], den[1], ..., den[n-1])
					for (result_p = result.begin(); 
					     result_p != result.end();
					     ++ result_p)
						r. lcm (lif,lif,result_p -> second);
					
					// filter out primes in PRIMEL from lif.					
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

					/*
					if ( prev == lif ) ++ count;
					else count = 1;
					*/
				}
				
				return lif;
			}
			
			/** memo Compute the last invariant factor.
			 */
			template<class IMatrix>
				Integer& lastInvariantFactor(Integer& lif, const IMatrix& A)  const {

				std::vector<Integer> empty_v;

				lastInvariantFactor (lif, A, empty_v);

				return lif;
			}
	
	};
}


#endif
