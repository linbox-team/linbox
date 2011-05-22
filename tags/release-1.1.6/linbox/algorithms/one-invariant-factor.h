/* File: ith-invaraint-factor.h
 *  Author: Zhendong Wan
 */
#ifndef __LINBOX__ITH_INVARIANT_FACTOR_H__
#define __LINBOX__ITH_INVARIANT_FACTOR_H__

#include <linbox/util/debug.h>
#include <linbox/algorithms/default.h>
#include <linbox/blackbox/compose.h>
#include <linbox/blackbox/random-matrix-traits.h>

namespace LinBox {
	
/// \brief Limited doc so far.
	template<class _Ring,
		class _LastInvariantFactor,
		class _Compose,
		class _RandomMatrix>
		
		class OneInvariantFactor {

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
			
			OneInvariantFactor(const Ring& _r = Ring(),
					   const LastInvariantFactor& _lif =LastInvariantFactor(),
					   const Compose& _compose =Compose(),
					   const RandomMatrix& _randomMatrix = RandomMatrix(),
					   int _threshold =DEFAULTOIFTHRESHOLD)
				:r(_r), lif(_lif), compose(_compose), randomMatrix(_randomMatrix), 
				threshold(_threshold), crossover(CROSSOVER) {
				
				if (_threshold <= 0) threshold = DEFAULTOIFTHRESHOLD;
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

			/** \brief Compute the i-th invariant factor of A, 
			 *  ignoring those factors of prime in PrimeL list.
			 *  It implements EGV++ (by bds), the adaptive algorithm of EGV and EGV+.
			 */
			template<class IMatrix, class Vector>
				Integer& oneInvariantFactor(Integer& oif, const IMatrix& A, 
							    int i, Vector& PrimeL) const {	
				
				// some check
				linbox_check(0 < i);
				linbox_check((unsigned int)i <= A.rowdim());
				linbox_check((unsigned int)i <= A.coldim());
				
				
			
				// if oif is the last invariant factor of A
				if ( ((unsigned int)i == A.rowdim()) && (A.rowdim() == A.coldim())) {

					lif.lastInvariantFactor(oif, A, PrimeL);
					return oif;
				}				
					
				r.init (oif, 0);

				int count;

				Integer prev, tmp_i;

				//typename RandomMatrixTraits<IMatrix>::value_type *L, *U;

				typename RandomMatrixTraits<IMatrix>::value_type *R, *L;
				
				typename ComposeTraits<IMatrix>::value_type* LAR;//*AUV; 


				// repeat threshold times
				for (count =0; count < threshold; ++ count) {

					r.assign (prev, oif);

					/*
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
					else {
					}
					*/
					// Always use LAR please refer ISSAC'04 paper by BDS and ZW
					randomMatrix.randomMatrix(L, r, i, A.rowdim());
						
					randomMatrix.randomMatrix(R, r, A.coldim(), i);
						
					compose.compose(LAR, *L, A, *R);
						
					lif.lastInvariantFactor(tmp_i, *LAR, PrimeL);

					//free memory
					delete L;
					delete R;
					delete LAR;
				
					r.gcd(oif, tmp_i, prev);
			
					// if oif reaches one
					if ( r.isOne(oif) ) break;		
					
				}
					
				return oif;
			}
			
		 	/** \brief Compute the i-th invariant factor of A.
			 *  It implements the adaptive algorithm of EGV and EGV+.
			 */
			template<class IMatrix>
				Integer& oneInvariantFactor(Integer& oif, const IMatrix& A, int i) const {   
				
				std::vector<Integer> empty_v;

				oneInvariantFactor (oif, A, i, empty_v);

				return oif;
			}

			/** \brief Compute the i-th invariant factor of A with bonus, 
			 *  ignoring those factors of prime in PrimeL list.
			 *  It implements EGV++ (by bds), the adaptive algorithm of EGV and EGV+.
			 */
			template<class IMatrix, class Vector>
				Integer& oneInvariantFactor_Bonus(Integer& oif, Integer& bonus, const IMatrix& A, 
							    int i, Vector& PrimeL) const {	
				// some check
				linbox_check(0 < i);
				linbox_check((unsigned int)i <= A.rowdim());
				linbox_check((unsigned int)i <= A.coldim());
				
				// if oif is the last invariant factor of A
				if ( ((unsigned int)i == A.rowdim()) && (A.rowdim() == A.coldim())) {
					lif.lastInvariantFactor_Bonus(oif, bonus, A, PrimeL);
					return oif;
				}				
					
				r.init (oif, 0); r. init (bonus, 0);
				int count;
				Integer prev, tmp_i, p_bonus;
				//typename RandomMatrixTraits<IMatrix>::value_type *L, *U;
				typename RandomMatrixTraits<IMatrix>::value_type *R, *L;
				typename ComposeTraits<IMatrix>::value_type* LAR;//*AUV; 
				// repeat threshold times
				for (count =0; count < threshold; ++ count) {
					r.assign (prev, oif); r. assign (p_bonus, bonus);
					// Always use LAR please refer ISSAC'04 papre by BDS and ZW
					randomMatrix.randomMatrix(L, r, i, A.rowdim());
					randomMatrix.randomMatrix(R, r, A.coldim(), i);
					compose.compose(LAR, *L, A, *R);
					lif.lastInvariantFactor_Bonus(tmp_i, bonus, *LAR, PrimeL);

					//free memory
					delete L;
					delete R;
					delete LAR;
					r. gcd(oif, tmp_i, prev);
					r. gcdin (bonus, p_bonus);
					// if oif reaches one
					if ( r.isOne(oif) ) break;		
					
				}
					
				return oif;
			}
			
		 	/** \brief Compute the i-th invariant factor of A.
			 *  It implements the adaptive algorithm of EGV and EGV+.
			 */
			template<class IMatrix>
				Integer& oneInvariantFactor_Bonus(Integer& oif, Integer& bonus, const IMatrix& A, int i) const {   
				
				std::vector<Integer> empty_v;

				oneInvariantFactor_Bonus (oif, bonus, A, i, empty_v);

				return oif;
			}


	};
}


#endif				
