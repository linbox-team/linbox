/** -*- mode:C++ -*- */
/** File: rational-solver.h
 *  Author: Zhendong Wan
 */

#ifndef __LINBOX_RATIONAL_SOLVER_H__
#define __LINBOX_RATIONAL_SOLVER_H__

#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/lifting-container.h>
#include <linbox/algorithms/rational-reconstruction.h>
#include <linbox/algorithms/matrix-inverse.h>
#include <linbox/algorithms/matrix-mod.h>


namespace LinBox {
	
	/** _Ring integer ring
	 *  _Field, finite field for lifting
	 */

	template<class _Ring,
		 class _Field,
		 class _RandomPrime>
		
	class RationalSolver {
		
	public:
			
		typedef _Ring Ring;
		typedef _Field Field;
		typedef _RandomPrime RandomPrime;

	protected:
		
		Ring r;

		RandomPrime rp;

	public:

		RationalSolver(const Ring& _r = Ring(), const RandomPrime& _rp = RandomPrime()) : 
			r (_r), rp (_rp) {}

		/** @memo solve Ax = b with coefficiencis in Ring 
		 *  over the quotient field of the Ring 
		 *  by lifting over Field
		 *  If A is not full rank, return 1,
		 *  otherwise return 0, and compute the solution 
		 *  and put in answer
		 *  answer is a vector of pair (num, den).
		 */
		template<class IMatrix, class Vector1, class Vector2>
		long solve(Vector1& answer, const IMatrix& A,
			   const Vector2& b) const {
			
			return solve (answer, A, b, false);
		}

		/** General case, risky algorithm
		 *  Mark oldMatrix is for optimal reason,
		 *  tell if A is the previous one
		 *  Suggest not to use it.
		 */
		template<class IMatrix, class Vector1, class Vector2>
		long solve(Vector1& answer, const IMatrix& A, 
			   const Vector2& b, bool oldMatrix) const {
			
			// history sensitive data for optimal reason
			static const IMatrix* IMP = 0;

			static DenseMatrix<Field>* FMP;
			
			static typename RandomPrime::Prime_Type prime;

			static long notfr;

			// if input matrix A is different one.
			if (!oldMatrix) {

				//delete IMP;
				
				delete FMP;

				IMP = &A;

				prime = rp. randomPrime();

				Field F(prime);
				
				MatrixMod::mod (FMP, A, F);
				
				notfr = MatrixInverse::matrixInverseIn(*FMP);

			}

			// not full rank
			if(notfr) 
				return 1;

			LiftingContainer<IMatrix, DenseMatrix<Field> >
				lc(*IMP, *FMP, b, prime);

			RationalReconstruction<LiftingContainer<IMatrix,
				DenseMatrix<Field> > > 
				re (lc);

			re.getRational(answer);

			return 0;
		}

	};
}

#endif
