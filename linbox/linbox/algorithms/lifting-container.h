/* -*- mode: C++ -*- */
/* File: lifting-container.h
 *  Author: Zhendong Wan
 */
#ifndef __LIFTING_CONTAINER_H__
#define __LIFTING_CONTAINER_H__

#include <linbox/util/debug.h>
#include <vector>
#include <algorithm>
#include <linbox/blackbox/apply.h>

namespace LinBox {

	/** @memo Use to solve Ax = b in the rational field by p-adic lifting. 
	 * @doc LiftingContainer takes two template parameters.
	 *  @param _IMatrix, matrix type over the integer ring
	 *  @param _FMatrix, matrix type over the finite field 
	 *  The member function next will output the p-adic digit
	 */

	template<class _IMatrix, class _FMatrix>
		class LiftingContainer {
		
		public:
			typedef _IMatrix IMatrix;
			typedef _FMatrix FMatrix;
			
			typedef typename IMatrix::Field Ring;
			typedef typename FMatrix::Field Field;
			
			typedef typename Ring::Element Integer;
			typedef typename Field::Element Element;	       
			
		protected:
			
			// an integer matrix, orginal data
			const IMatrix& A;
			
			// Should be the inverse of A mod p
			const FMatrix& Ap;
			
			// an object of integer ring
			const Ring r;
			
			// finite field 
			const Field f;
			
			// prime
			Integer p;
		
			// res = (b - Ax) / m, m = current modulus.
			mutable std::vector<Integer> res;

			const VectorDomain<Ring> vd;
			
		public:
			
			// return size of solution
			inline int size() const { return res.size(); }
			
			// return the field
			const Field& field() const { return f; }

			// return the ring
			const Ring& ring() const { return r; }

			// return the prime
			const Integer& prime () const { return p; }

			/** @memo constructed from an integer matrix, 
			 *  its inverse over finite field,  a vector, 
			 *  a prime.
			 */		       
			
			template<class Vector, class Prime_Type>
			LiftingContainer (const IMatrix& _A, const FMatrix& _Ap, const Vector& b, const Prime_Type& _p)
				: A(_A), Ap(_Ap), r (_A.field()), f (_Ap.field()), vd (_A.field()) {
				
				// linbox check if A is a square matrix
				//linbox_check(A.rowdim() == A.coldim());
 			
				linbox_check(_A.coldim() == b.size());
				
				r. init (p, _p);
			
				typename Vector::const_iterator b_p;
				typename std::vector<Integer>::iterator p0;
				
				// initialize res = b
				res.resize (A. coldim());

				for (p0 = res.begin(), b_p = b.begin();
				     p0 != res.end(); ++ p0, ++ b_p)
					
					r.init(*p0, *b_p);
			}



			/** @memo next (digit) 
			 *  Generate next p-adic digit of solution fo x Ax = b.
			 */
			template<class Vector>
			inline Vector& next(Vector& digit) const {
				
				linbox_check (digit.size() == A. rowdim());
				
				std::vector<Element> v1 (A. coldim());

				typename std::vector<Element>::iterator p1;

				typename std::vector<Integer>::iterator p0;
				
				// v1 =  res mod p
				for ( p1 = v1. begin(), p0 = res. begin();
				      p0 != res. end(); ++ p0, ++ p1) 
					
					f. init (*p1, *p0);

				// compute next p-adic digit
				apply (digit, Ap, v1);
			
				// prepare for computing next digit				
				// update tmp_integerv = A * digit
				
				std::vector <Integer> v2 (A. coldim());
				
				// v2 = A digit
				apply (v2, A, digit);
				
				// update res -= v2 
				vd. subin (res, v2);

				// update res = res / p
				for ( p0 = res.begin(); p0 != res.end(); ++ p0) {
					
					r.divin(*p0, p);

				}

				return digit;

			}
		};

	
}       

#endif
