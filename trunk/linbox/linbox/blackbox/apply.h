/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* File: apply.h
 *  Author: Zhendong Wan
 */

/* Reserve for possible optimal.
 */

#ifndef __LINBOX_APPLY_H__
#define __LINBOX_APPLY_H__

#include <linbox-config.h>
#include <linbox/integer.h>
#include <linbox/util/debug.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/field/ntl-ZZ.h>
#include <vector>

#ifdef __LINBOX_BLAS_AVAILABLE
#include <linbox/fflas/fflas.h>
#endif

namespace LinBox {
	
	// general case, y = A x
	template<class OutV, class Matrix, class InV>
	inline OutV& apply (OutV& y, const Matrix& A, const InV& x) {
		
		return A. apply (y, x);

	}

	template<class OutV, class Matrix, class InV>
	inline OutV& applyTranspose (OutV& y, const Matrix& A, const InV& x) {

		return A. applyTranspose (y, x);
		
	}

	
	template<class Domain>
	class BlasApply {
	  	 
	public:
		typedef typename Domain::Element    Element;
		typedef std::vector<Element>         Vector;

		BlasApply(const Domain& D) : _D(D), _MD(D) { 
			_D.characteristic(_prime); 
			_D.init(_one,1UL); 
			_D.init(_zero,0UL);
		}
	  	  		

		//#ifdef __LINBOX_BLAS_AVAILABLE
		inline Vector& applyV(Vector                        &y,
				     const BlasMatrix<Element>     &A, 
				     const Vector                  &x) const {
	    
			if (( _prime > 0) && ( _prime <  67108863)) {	    	   				

				FFLAS::fgemv( _D, FFLAS::FflasNoTrans, 
					      A.rowdim(), A.coldim(),
					      _one,
					      A.getPointer(), A.getStride(),
					      &x[0],1,
					      _zero,
					      &y[0],1);  	      
			}
			else {
				_MD.vectorMul (y, A, x);	      
			}
			return y;
		}		
		//#endif
	  
		template<class OutV, class Matrix, class InV>
		inline OutV& applyV(OutV           &y, 
				   const Matrix   &A, 
				   const InV      &x) const {
			return apply(y,A,x);
		}


	private:
		Domain         _D;
		integer    _prime;
		Element _one,_zero;
		MatrixDomain<Domain> _MD;
	  

	};

} // end of namespace LinBox		
#endif
