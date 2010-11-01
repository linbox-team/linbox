/* elim.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.waterloo.ca>
 *
 * --------------------------------------------
 *
 * Licensed under the GNU Lesser General Public License. See COPYING for
 * details.
 *
 * Elimination code for lookahead block Lanczos
 */

#ifndef __LINBOX_eliminator_H
#define __LINBOX_eliminator_H

#include "linbox/linbox-config.h"

#include <vector>

#include "linbox/field/archetype.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/dense.h"
#include "linbox/matrix/dense-submatrix.h"
#include "linbox/solutions/methods.h"

// Fix for Solaris wierdness
#undef _S
#undef _M
#undef _N
#undef _P
#undef _U
#undef _T

namespace LinBox 
{

/** Elimination system
 *
 * This is the supporting elimination system for a lookahead-based
 * variant of block Lanczos.
 */
template <class Field, class Matrix = DenseMatrixBase<typename Field::Element> >
class Eliminator
{
    public:

	typedef typename Field::Element Element;

	/** Permutation
	 *
	 * A permutation is represented as a vector of pairs, each
	 * pair representing a transposition. Thus a permutation
	 * requires O(n log n) storage and O(n log n) application
	 * time, as opposed to the lower bound of O(n) for
	 * both. However, this allows us to decompose a permutation
	 * easily into its factors, thus eliminating the need for
	 * additional auxillary storage in each level of the
	 * Gauss-Jordan transform recursion. Additionally, we expect
	 * to use this with dense matrices that are "close to
	 * generic", meaning that the rank should be high and there
	 * should be relatively little need for transpositions. In
	 * practice, we therefore expect this to beat the vector
	 * representation. The use of this representation does not
	 * affect the analysis of the Gauss-Jordan transform, since
	 * each step where a permutation is applied also requires
	 * matrix multiplication, which is strictly more expensive.
	 */
	typedef std::pair<unsigned int, unsigned int> Transposition;
	typedef std::vector<Transposition> Permutation;

	/** Constructor
	 * @param F Field over which to operate
	 * @param traits @ref{SolverTraits} structure describing user
	 *               options for the solver 
	 */
	Eliminator (const Field &F, unsigned int N);

	/** Destructor
	 */
	~Eliminator ();

	/** Two-sided Gauss-Jordan transform
	 *
	 * @param Ainv Inverse of nonsingular part of A
	 * @param Tu Row dependencies
	 * @param Tv Column dependencies
	 * @param P Row permutation
	 * @param Q Column permutation
	 * @param A Input matrix
	 * @param rank Rank of A
	 */
	template <class Matrix1, class Matrix2, class Matrix3, class Matrix4>
	void twoSidedGaussJordan (Matrix1       &Ainv,
				  Permutation   &P,
				  Matrix2       &Tu,
				  Permutation   &Q,
				  Matrix3       &Tv,
				  const Matrix4 &A,
				  unsigned int  &rank);

	/** Permute the input and invert it
	 *
	 * Compute the pseudoinverse of the input matrix A and return
	 * it. First apply the permutation given by the lists leftPriorityIdx
	 * and rightPriorityIdx to the input matrix so that independent
	 * columns and rows are more likely to be found on the first indices
	 * in those lists. Zero out the rows and columns of the inverse
	 * corresponding to dependent rows and columns of the input. Set S and
	 * T to boolean vectors such that S^T A T is invertible and of maximal
	 * size.
	 *
	 * @param W Output inverse
	 * @param S Output vector S
	 * @param T Output vector T
	 * @param leftPriorityIdx Priority indices on the left
	 * @param Rightpriorityidx Priority indices on the right
	 * @param A Input matrix A
	 * @return Reference to inverse matrix
	 */
	Matrix &permuteAndInvert (Matrix                  &W,
				  std::vector<bool>       &S,
				  std::vector<bool>       &T,
				  std::list<unsigned int> &rightPriorityIdx,
				  Permutation             &Qp,
				  unsigned int            &rank,
				  const Matrix            &A);

	/** Perform a Gauss-Jordan transform using a recursive algorithm
	 *
	 * Upon completion, we have UPA = R, where R is of reduced row
	 * echelon form
	 *
	 * @param U Output matrix U
	 * @param P Output permutation P
	 * @param A Input matrix A
	 * @return Reference to U
	 */
	template <class Matrix1, class Matrix2, class Matrix3, class Matrix4>
	Matrix1 &gaussJordan (Matrix1                   &U,
			      std::vector<unsigned int> &profile,
			      Permutation               &P,
			      Matrix2                   &Tu,
			      Permutation               &Q,
			      Matrix3                   &Tv,
			      unsigned int              &rank,
			      typename Field::Element   &det,
			      const Matrix4             &A);

	/** Retrieve the total user time spent permuting and inverting
	 */
	double getTotalTime () const { return _total_time; }

	/** Retrieve the total user time spent inverting only
	 */
	double getInvertTime () const { return _invert_time; }

	/** Write the filter vector to the given output stream
	 */
	std::ostream &writeFilter (std::ostream &out, const std::vector<bool> &v) const;

	/** Write the given permutation to the output stream
	 */
	std::ostream &writePermutation (std::ostream &out, const Permutation &P) const;

    private:
	// Compute the kth indexed Gauss-Jordan transform of the input
	Matrix &kthGaussJordan (unsigned int                  &r,
				typename Field::Element       &d,
				unsigned int                   k,
				unsigned int                   s,
				unsigned int                   m,
				const typename Field::Element &d0);

	// Set the given matrix to the identity
	template <class Matrix1>
	Matrix1 &setIN (Matrix1 &A) const;

	// Add d * I_N to A
	template <class Matrix1>
	Matrix1 &adddIN (Matrix1                       &A,
			 const typename Field::Element &d) const;

	// Clean out the given priority index list and add new elements as needed
	void cleanPriorityIndexList (std::list<unsigned int> &list,
				     std::vector<bool>       &S,
				     std::vector<bool>       &old_S) const;

	// Permute the given bit vector
	template <class Iterator>
	std::vector<bool> &permute (std::vector<bool>  &v,
				    Iterator            P_start,
				    Iterator            P_end) const;

	// Construct a permutation from the given priority list
	Permutation &buildPermutation (Permutation &P, const std::list<unsigned int> &pidx) const;

	// Prepare a minimal permutation based on the given permutation
	Permutation &buildMinimalPermutation (Permutation &P, unsigned int rank,
					      unsigned int dim, const Permutation &Pold);

	Permutation &buildMinimalPermutationFromProfile (Permutation &P, unsigned int rank,
							 unsigned int dim, const std::vector<unsigned int> &profile);

	// Private variables

	const Field                      &_F;
	VectorDomain<Field>               _VD;
	MatrixDomain<Field>               _MD;
	unsigned int                      _N;

	typename Field::Element           _one;

	// Temporaries used in the computation

	mutable Permutation               _P;

	mutable DenseMatrixBase<Element>  _A;         // Variable
	mutable DenseMatrixBase<Element>  _U;         // Variable
	mutable DenseMatrixBase<Element>  _tmp;

	// These record the independent rows and columns found during the
	// elimination process

	mutable std::vector<bool>         _S;         // Independent rows
	mutable std::vector<bool>         _T;         // Independent columns

	std::vector<unsigned int>         _profile;
	unsigned int                      _profile_idx;

	// Timer information

	double                            _total_time;
	double                            _invert_time;

	// Priority indices for rows
	std::vector<unsigned int>         _indices;
};

} // namespace LinBox

#include "eliminator.inl"

#endif // __LINBOX_eliminator_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
