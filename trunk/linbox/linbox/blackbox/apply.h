/** File: apply.h
 *  Author: Zhendong Wan
 */

/** Reserve for possible optimal.
 */

#ifndef __LINBOX_APPLY_H__
#define __LINBOX_APPLY_H__

#include <linbox/util/debug.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/field/ntl-ZZ.h>
#include <vector>

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


}		
#endif
