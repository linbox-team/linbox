/** File: ith-invaraint-factor.h
 *  Author: Zhendong Wan
 */
#ifndef __LINBOX__ITH_INVARIANT_FACTOR_H__
#define __LINBOX__ITH_INVARIANT_FACTOR_H__

#include <linbox/util/debug.h>
#include <linbox/algorithms/default.h>
#include <linbox/blackbox/compose-traits.h>
#include <linbox/blackbox/random-matrix-traits.h>

namespace LinBox {
	
	template<class _Ring,
		class _LastInvariantFactor,
		class _Compose,
		class _RandomMatrix>
		
		class IthInvariantFactor {

		public:
		
		typedef _LastInvariantFactor LastInvariantFactor;
		
		typedef _Ring Ring;

		typedef _Compose Compose;

		typedef _RandomMatrix RandomMatrix;
		
		typedef typename Ring::Element Integer;
		
		protected:
	
			Ring r;
			
			LastInvariantFactor lif;

			Compose compose;

			RandomMatrix randomMatrix;
			
			int threshold;

			double crossover;
			
			public:			
			
			IthInvariantFactor(const Ring& _r = Ring(),
					   const LastInvariantFactor& _lif =LastInvariantFactor(),
					   const Compose& _compose =Compose(),
					   const RandomMatrix& _randomMatrix = RandomMatrix(),
					   int _threshold =DEFAULTIIFTHRESHOLD)
				:r(_r), lif(_lif), compose(_compose), randomMatrix(_randomMatrix), 
				threshold(_threshold), crossover(CROSSOVER) {
				
				if (_threshold <= 0) threshold = DEFAULTIIFTHRESHOLD;
			}
			
			void setThreshold (int _threshold) {
				if (_threshold > 0) {
					threshold = _threshold; 
				}
				
			}
			
			int getThreshold () const {
				return threshold;
			}
			
			void setCrossover(double t) {
				if(0 <= t <= 1)
					crossover = t;

			}

			double getCrossover() const {
				return crossover;
			}

			LastInvariantFactor& getLastInvariantFactor () {
				return lif;
			}
			
			const LastInvariantFactor& getLastInvariantFactor () const {
				return lif;
			}

			/** @memo Compute the i-th invariant factor of A, 
			 *  missing those factors of prime in PrimeL list.
			 *  It implements the adaptive algorithm of EGV and EGV+.
			 */
			template<class IMatrix, class Vector>
				Integer& ithInvariantFactor(Integer& iif, const IMatrix& A, 
							    int i, Vector& PrimeL) const {	
				
				// some check
				linbox_check(0 < i);
				linbox_check((unsigned int)i <= A.rowdim());
				linbox_check((unsigned int)i <= A.coldim());
				
				
			
				// if iif is the last invariant factor of A
				if ( ((unsigned int)i == A.rowdim()) && (A.rowdim() == A.coldim())) {

					lif.lastInvariantFactor(iif, A, PrimeL);
					return iif;
				}				
					
				r.init (iif, 0);

				int count;

				Integer prev, tmp_i;

				typename RandomMatrixTraits<IMatrix>::value_type *L, *U;

				typename RandomMatrixTraits<IMatrix>::value_type *R, *V;
				
				typename ComposeTraits<IMatrix>::value_type* LAR, *AUV; 


				// repeat threshold times
				for (count =0; count < threshold; ++ count) {

					r.assign (prev, iif);

					// Use A + UV
					if ((A.rowdim() == A.coldim()) && (i > crossover * A.rowdim())) {
						
						randomMatrix.randomMatrix(U, r, A.rowdim(), A.coldim() - i);
						
						randomMatrix.randomMatrix(V, r, A.rowdim() - i, A.coldim());
						
						compose.composeBig(AUV, A, *U, *V);

						// compute the last invariant factor of RAL
						lif.lastInvariantFactor(tmp_i, *AUV, PrimeL);
					
						// free memory.
						delete U;
						delete V;
						delete AUV;
						
					}
					// Use LAR
					else {
						randomMatrix.randomMatrix(L, r, i, A.rowdim());
						
						randomMatrix.randomMatrix(R, r, A.coldim(), i);
						
						compose.composeSmall(LAR, *L, A, *R);
						
						lif.lastInvariantFactor(tmp_i, *LAR, PrimeL);

						//free memory
						delete L;
						delete R;
						delete LAR;
					}
				
						
					r.gcd(iif, tmp_i, prev);
			
					// if iif reaches one
					if ( r.isOne(iif) ) break;		
					
				}
					
				return iif;
			}
			
		 	/** @memo Compute the i-th invariant factor of A.
			 *  It implements the adaptive algorithm of EGV and EGV+.
			 */
			template<class IMatrix>
				Integer& ithInvariantFactor(Integer& iif, const IMatrix& A, int i) const {   
				
				std::vector<Integer> empty_v;

				ithInvariantFactor (iif, A, i, empty_v);

				return iif;
			}


	};
}


#endif				

			

			
