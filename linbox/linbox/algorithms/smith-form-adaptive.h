/* Implement the adaptive algorithm for Smith form compuation
 * authors: bds and zw
 */
#ifndef __LINBOX_SMITH_FORM_ADAPTIVE_H__
#define __LINBOX_SMITH_FORM_ADAPTIVE_H__

#include <vector>
#include <linbox/integer.h>
#include <linbox/blackbox/dense.h>

namespace LinBox {

	long prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

	const int NPrime = 25;

	/* Compute the local smith form at prime p, when modular (p^e) fits in long
	*/
	template <class Ring>
	void compute_local_long (std::vector<integer>& s, const DenseMatrix<Ring>& A, long p, long e);

	/* Compute the local smith form at prime p, when modular (p^e) doesnot fit in long
	*/
	template <class Ring>
	void compute_local_big (std::vector<integer>& s, const DenseMatrix<Ring>& A, long p, long e);

	/* Compute the local smith form at prime p
	*/
	template <class Ring>
	void compute_local (std::vector<integer>& s, const DenseMatrix<Ring>& A, long p, long e);

	/* Compute the k-smooth part of the invariant factor, where k = 100.
	 * @param sev is the exponent part ...
	 * By local smith form and rank computation
	 */
	template <class Ring>
	void SmithFormSmooth (std::vector<integer>& s, const DenseMatrix<Ring>& A, long r, const std::vector<long>& sev);
			
	/* Compute the k-rough part of the invariant factor, where k = 100.
	 * By EGV+ algorithm or Iliopoulos' algorithm for Smith form.
	*/
	template <class Ring>
	void SmithFormRough  (std::vector<integer>& s, const DenseMatrix<Ring>& A, integer m );

	/** \brief Smith form of a dense matrix by adaptive algorithm.
	 *
	 * Compute the largest invariant factor, then, based on that, 
	 * compute the rough and smooth part, separately.
	 */
	template <class Ring>
	void SmithForm (std::vector<integer>& s, const DenseMatrix<Ring>& A);
}

#include <linbox/algorithms/smith-form-adaptive.inl>
#endif
