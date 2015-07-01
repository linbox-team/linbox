/* Implement the adaptive algorithm for Smith form compuation
 * authors: bds and zw
 */
#ifndef __LINBOX_SMITH_FORM_ADAPTIVE_H__
#define __LINBOX_SMITH_FORM_ADAPTIVE_H__

#include <vector>
#include <linbox/integer.h>
#include <linbox/blackbox/dense.h>

namespace LinBox {

class SmithFormAdaptive {
	public:

	static const long prime[];// = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

	static const int NPrime;// = 25;

	/* Compute the local smith form at prime p, when modular (p^e) fits in long
	 * Should work with SparseMatrix and DenseMatrix
	*/
	template <class Matrix>
	static void compute_local_long (std::vector<integer>& s, const Matrix& A, long p, long e);

	/* Compute the local smith form at prime p, when modular (p^e) doesnot fit in long
	 * Should work with SparseMatrix and DenseMatrix
	*/
	template <class Matrix>
	static void compute_local_big (std::vector<integer>& s, const Matrix& A, long p, long e);

	/* Compute the local smith form at prime p
	*/
	template <class Matrix>
	static void compute_local (std::vector<integer>& s, const Matrix& A, long p, long e);

	/* Compute the k-smooth part of the invariant factor, where k = 100.
	 * @param sev is the exponent part ...
	 * By local smith form and rank computation
	 * Should work with SparseMatrix and DenseMatrix
	 */
	template <class Matrix>
	static void smithFormSmooth (std::vector<integer>& s, const Matrix& A, long r, const std::vector<long>& sev);
			
	/* Compute the k-rough part of the invariant factor, where k = 100.
	 * By EGV+ algorithm or Iliopoulos' algorithm for Smith form.
	 * Should work with DenseMatrix
	*/
	template <class Matrix>
	static void smithFormRough  (std::vector<integer>& s, const Matrix& A, integer m );

	/* Compute the Smith form via valence algorithms
	 * Compute the local Smith form at each possible prime
	 * r >= 2;
	 * Should work with SparseMatrix and DenseMatrix
	 */
	template <class Matrix>
	static void smithFormVal (std::vector<integer>&s, const Matrix& A, long r, const std::vector<long>& sev);

	/** \brief Smith form of a dense matrix by adaptive algorithm.
	 *
	 * Compute the largest invariant factor, then, based on that, 
	 * compute the rough and smooth part, separately.
	 * Should work with SparseMatrix and DenseMatrix
	 */
	template <class Matrix>
	static void smithForm (std::vector<integer>& s, const Matrix& A);
	/** Specialization for dense case*/
	template <class IRing>
	static void smithForm (std::vector<integer>& s, const DenseMatrix<IRing>& A);
};
	const long SmithFormAdaptive::prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
	const int SmithFormAdaptive::NPrime = 25;
}

#include <linbox/algorithms/smith-form-adaptive.inl>
#endif
