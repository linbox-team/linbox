/** File: smith.h
 *  Author: Zhendong Wan
 *  Implementation of EGV and EGV+ algorithm
 */

#ifndef __LINBOX__SMITH_FORM_H__
#define __LINBOX__SMITH_FORM_H__

#include <linbox/util/debug.h>
#include <linbox/algorithms/default.h>

namespace LinBox {
	
	template <class _Ring,
		  class _IthInvariantFactor,
		  class _Rank>
		
		class SmithForm {
			
		public:
		
		typedef _Ring Ring;
		
		typedef _IthInvariantFactor IthInvariantFactor;
		
		typedef _Rank Rank;		
		
		typedef typename Ring::Element Integer;

		protected:
		
		IthInvariantFactor iif;
		Rank rank;
		Ring r;

		public:
		
		// constructor
		SmithForm(const IthInvariantFactor& _iif =IthInvariantFactor(),
			  const Rank& _rank =Rank(), const Ring& _r = Ring(),
			  int _iifthreshold =DEFAULTIIFTHRESHOLD, int _lifthreshold =DEFAULTLIFTHRESHOLD)
			
			: iif(_iif),rank(_rank),r(_r) { 
			
			iif.setThreshold(_iifthreshold);
			iif.getLastInvariantFactor().setThreshold( _lifthreshold);
		}
			
		void setIIFThreshold (int _iifthreshold =DefaultIIFThreshold) {
			iif.setThreshold(_iifthreshold);
		}
		
		void setLIFThreshold (int _lifthreshold =DefaultLIFThreshold) {
			iif.getLastInvariantFactor().setThreshold(_lifthreshold);
		}
		
		int getIIFThreshold() const {
			return iif.getThreshold();
		}
		
		int getLIFThreshold() const {
			return iif.getLastInvariantFactor().getlIFThreshold();
		}
		
		
		/** @memo compute the Smith Form of an integer matrix,
		 *  missing these factors of primes in PrimeL
		 */
		template<class IMatrix, class Vector, class VectorP>
			Vector&  smithForm(Vector& sf, const IMatrix& A, const VectorP& PrimeL) const{
				
			// check if there are enough spaces in sf to store all invariant factors of A
			linbox_check(sf.size() >= (A.rowdim() <= A.coldim() ? A.rowdim() : A.coldim()));
						
			typename Vector::iterator p;

			Integer zero;
				
			long Ar = rank.rank(A);

			r.init (zero,0);
				
			// set k-th invariant factor to zero for all k > Ar
			for (p = sf.begin() + Ar; p!= sf.end(); ++p)
				r.assign(*p,zero);

			// A is a zero matrix
			if (Ar == 0)
				return sf;
				
				
			// compute first invariant factor of A
			firstInvariantFactor(sf[0], A, PrimeL);
				
			// if rank(A) == 1
			if (Ar == 1) 
				return sf;

				
			iif.ithInvariantFactor(sf[Ar - 1], A, Ar, PrimeL);
				
			// binary search smith form
			smithFormBinarySearch (sf, A, 1, Ar, PrimeL);

			return sf;
		}

			
		/** @memo compute the Smith Form of an integer matrix
		 */
		template<class IMatrix, class Vector>
			Vector&  smithForm(Vector& sf, const IMatrix& A) const{

			std::vector<Integer> empty_v;
				
			smithForm (sf, A, empty_v);

			return sf;

		}

		protected:			

		/** @memo compute the 1st invariant factor, = GCD (all element in A),
		 *  missing these factors of primes in PrimeL
		 */
		template<class IMatrix, class Vector>
			Integer& firstInvariantFactor(Integer& fif, const IMatrix& A, const Vector& PrimeL) const {

			r.init(fif,0);			
				
			typename IMatrix::ConstRawIterator A_p;

			for (A_p = A.rawBegin(); A_p != A.rawEnd(); ++ A_p) {
				
				if (!r.isZero(*A_p)) {
					r.gcd(fif, fif, *A_p);
					
					// if tmp == 1, break
					if (r.isOne(fif)) return fif;
				}
			}


			if (r.isZero(fif)) return fif;

			Integer p, quo, rem;

			typename Vector::const_iterator Prime_p;
			       
			// filter out primes in PRIME from lif
			for ( Prime_p = PrimeL.begin(); Prime_p != PrimeL.end(); ++ Prime_p) {
					
				r.init (p, *Prime_p);
					
				do {
					r.quoRem(quo,rem,fif,p);
							
					if (r.isZero( rem )) r.assign(fif,quo);
					else break;
				}
				while (true);

			}												
				
			return fif;
				
		}
		

		/** @memo Binary search invariant factors between i and j, missing those factors in PrimeL
		 *  suppose sf[i - 1], sf [j - 1] are ith and jth invariant factor of A
		 *  i <= j
		 */

		template<class IMatrix, class Vector, class VectorP>
			Vector& smithFormBinarySearch (Vector& sf, const IMatrix& A, int i, int j, const VectorP& PrimeL) const {

#ifdef WANDEBUG
			std::cout << "Binary Search invariant factors [" << i << ", "<< j << "]\n " << std::flush;
#endif
				
			typename Vector::iterator p;

			// if no invariant factor between i and j
			if (j <= i + 1) return sf;
			
			// if i-th invariant factor == j-th invariant factor
			if (r.areEqual(sf[i - 1], sf[j - 1])) {
				for (p = sf.begin() + i; p != sf.begin() + (j -1); ++ p)
					r.assign (*p, sf[i-1]);
				return sf;
			}

			int mid = (i + j) / 2;
#ifdef WANDEBUG
			std::cout << "Start to compute " << mid << "-th invariant factor:\n" << std::flush;
#endif 
			
			iif.ithInvariantFactor (sf[mid - 1], A, mid, PrimeL);


#ifdef WANDEBUG
			std::cout << mid <<"-th invairant factor of A = " ;
			r.write (std::cout, sf[mid -1]);
			std::cout << "\n" << std::flush;
#endif

			// recursively binary search all k-invariant factors, where i <= k <= mid
			smithFormBinarySearch (sf, A, i, mid, PrimeL);

			// recurseively binary search all k-invariant factors, where mid <= k < j
			smithFormBinarySearch (sf, A, mid, j, PrimeL);
			
			return sf;
		}

	};

}

#endif 

			
			

			
			
				
			


			
		

		
