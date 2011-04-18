/* eliminator.inl
 * Copyright (C) 2002, 2003 LinBox, Bradford Hovinen
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

#ifndef __LINBOX_eliminator_INL 
#define __LINBOX_eliminator_INL

#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/util/debug.h"
#include "linbox/solutions/methods.h"
#include "linbox/matrix/dense-submatrix.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"

#include "eliminator.h"

#undef _F
#undef _N
#undef _S
#undef _T

namespace LinBox 
{

std::ostream &reportPermutation
	(std::ostream &out,
	 const std::vector<std::pair<unsigned int, unsigned int> > &P) 
{
	std::vector<std::pair<unsigned int, unsigned int> >::const_iterator i;

	out << "  ";

	for (i = P.begin (); i != P.end (); ++i)
		out << "(" << i->first << " " << i->second << ")";

	out << std::endl;

	return out;
}

template <class Field, class Matrix>
Eliminator<Field, Matrix>::Eliminator (const Field &F, unsigned int N) 
	: _F (F), _VD (F), _MD (F), _N (N), _S (N), _T (N)
{
	_F.init (_one, 1);
}

template <class Field, class Matrix>
Eliminator<Field, Matrix>::~Eliminator () 
{
}

template <class Field, class Matrix>
template <class Matrix1, class Matrix2, class Matrix3, class Matrix4>
void Eliminator<Field, Matrix>::twoSidedGaussJordan
	(Matrix1       &Ainv,
	 Permutation   &P,
	 Matrix2       &Tu,
	 Permutation   &Q,
	 Matrix3       &Tv,
	 const Matrix4 &A,
	 unsigned int  &rank)
{
	typename Field::Element d, dinv;
	std::vector<unsigned int>::iterator i;
	unsigned int idx;

	_A.resize (A.rowdim (), A.coldim ());
	_U.resize (A.rowdim (), A.rowdim ());
	_tmp.resize (A.rowdim (), A.coldim ());
	_profile.resize (A.coldim ());
	_indices.resize (A.rowdim ());

	for (i = _indices.begin (), idx = 0; i != _indices.end (); ++i, ++idx)
		*i = idx;

	_MD.subin (_A, _A);

	DenseSubmatrix<Element> A1 (_A, 0, 0, A.rowdim (), A.coldim ());
	_MD.copy (A1, A);

	setIN (_U);
	_P.clear ();
	_profile_idx = 0;

	kthGaussJordan (rank, d, 0, 0, A.coldim (), _one);

	buildMinimalPermutation (P, rank, A.rowdim (), _P);
	buildMinimalPermutationFromProfile (Q, rank, A.coldim (), _profile);

	_MD.permuteColumns (_U, _P.rbegin (), _P.rend ());
	_MD.permuteColumns (_U, P.begin (), P.end ());

	DenseSubmatrix<Element> Tu1 (Tu, rank, 0, A.rowdim () - rank, rank);
	DenseSubmatrix<Element> U2 (_U, rank, 0, A.rowdim () - rank, rank);
	_MD.copy (Tu1, U2);

	DenseSubmatrix<Element> Ainv1 (Ainv, 0, 0, rank, rank);
	DenseSubmatrix<Element> U1 (_U, 0, 0, rank, rank);
	_F.inv (dinv, d);
	_MD.mul (Ainv1, U1, dinv);

	DenseSubmatrix<Element> Tv1 (Tv, 0, rank, rank, A.coldim () - rank);
	DenseSubmatrix<Element> U3 (_U, 0, 0, rank, A.rowdim ());
	DenseSubmatrix<Element> A2 (_A, 0, rank, A.rowdim (), A.coldim () - rank);
	_MD.copy (A1, A);
	_MD.permuteColumns (A1, Q.begin (), Q.end ());
	_MD.permuteRows (A2, P.begin (), P.end ());
	_MD.mul (Tv1, U3, A2);
	_MD.negin (Tv1);
}

/* permuteAndInvert
 *
 * Compute the pseudoinverse of the input matrix A and return
 * it. First apply the permutation given by the lists leftPriorityIdx
 * and rightPriorityIdx to the input matrix so that independent
 * columns and rows are more likely to be found on the first indices
 * in those lists. Zero out the rows and columns of the inverse
 * corresponding to dependent rows and columns of the input. Set S and
 * T to boolean vectors such that S^T A T is invertible and of maximal
 * size.
 */

template <class Field, class Matrix>
Matrix &Eliminator<Field, Matrix>::permuteAndInvert
	(Matrix                    &W,
	 std::vector<bool>         &S,
	 std::vector<bool>         &T,
	 std::list<unsigned int>   &rightPriorityIdx,
	 Permutation               &Qp,
	 unsigned int              &rank,
	 const Matrix              &A)
{
	typename Field::Element d;     // Determinant of input A, up to sign

	typename Matrix::ConstRowIterator ai;
	typename Matrix::RowIterator _ai, wi, ui;

	std::vector<unsigned int>::iterator i;
	unsigned int idx;

 	Timer timer;

 	timer.start ();

	linbox_check (_N == S.size ());
	linbox_check (_N == T.size ());
	linbox_check (_N == W.rowdim ());
	linbox_check (_N == W.coldim ());
	linbox_check (_N == A.rowdim ());
	linbox_check (_N == A.coldim ());

#ifdef ELIM_DETAILED_TRACE
	commentator.start ("Computing W, S, T", "Eliminator::permuteAndInvert", _N);

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix:" << std::endl;
	_MD.write (report, A);
#endif

	/* Apply initial permutations to A, copying to _A */

	buildPermutation (Qp, rightPriorityIdx); // Column permutation

#ifdef ELIM_DETAILED_TRACE
	report << "Column permutation: ";
	reportPermutation (report, Qp) << endl;
#endif

	_A.resize (A.rowdim (), A.coldim ());
	_U.resize (A.rowdim (), A.rowdim ());
	_tmp.resize (A.rowdim (), A.coldim ());
	_indices.resize (A.rowdim ());

	for (i = _indices.begin (), idx = 0; i != _indices.end (); ++i, ++idx)
		*i = idx;

	_MD.copy (_A, A);

	_MD.permuteColumns (_A, Qp.begin (), Qp.end ());

	/* Initialise temporaries for the computation */

	setIN (_U);
	_P.clear ();
	_profile.resize (A.coldim ());
	_profile_idx = 0;
	std::fill (_T.begin (), _T.end (), false);

	/* Run the computation */

	kthGaussJordan (rank, d, 0, 0, _A.coldim (), _one);
	
	/* Set _S based on the permutation */
	std::fill (_S.begin (), _S.begin () + rank, true);
	std::fill (_S.begin () + rank, _S.end (), false);
	permute (_S, _P.rbegin (), _P.rend ());

	permute (_T, Qp.rbegin (), Qp.rend ());

	/* Apply final permutations to _U, copying to W */
	DenseSubmatrix<Element> U1 (_U, rank, 0, _U.rowdim () - rank, _U.coldim ());
	_MD.subin (U1, U1);

	_MD.permuteColumns (_U, _P.rbegin (), _P.rend ());

	/* Divide _U by the determinant and copy to W */
	_F.invin (d);
	_MD.mulin (_U, d);

	_MD.subin (W, W);

	typename std::vector<unsigned int>::iterator pi;

	for (pi = _profile.begin (), ui = _U.rowBegin (); pi != _profile.begin () + rank; ++ui, ++pi)
		_VD.copy (*(W.rowBegin () + *pi), *ui);

//	_MD.permuteRows (W, Qp.rbegin (), Qp.rend ());

	/* Rebuild leftPriorityIdx and rightPriorityIdx */
	cleanPriorityIndexList (rightPriorityIdx, _T, T);

	/* Reverse the row priority index list */
	std::reverse (_indices.begin (), _indices.end ());

	S = _S;
	T = _T;

#ifdef ELIM_DETAILED_TRACE
	report << "Computed W:" << std::endl;
	_MD.write (report, W);

	commentator.stop ("done", NULL, "Eliminator::permuteAndInvert");
#endif

 	timer.stop ();
 	_total_time += timer.usertime ();

	return W;
}

template <class Field, class Matrix>
template <class Matrix1, class Matrix2, class Matrix3, class Matrix4>
Matrix1 &Eliminator<Field, Matrix>::gaussJordan
	(Matrix1                   &U,
	 std::vector<unsigned int> &profile,
	 Permutation               &P,
	 Matrix2                   &Tu,
	 Permutation               &Q,
	 Matrix3                   &Tv,
	 unsigned int              &rank,
	 typename Field::Element   &det,
	 const Matrix4             &A)
{
	typename Field::Element dinv;
	std::vector<unsigned int>::iterator i;
	unsigned int idx;

	_A.resize (A.rowdim (), A.coldim ());
	_U.resize (A.rowdim (), A.rowdim ());
	_tmp.resize (A.rowdim (), A.coldim ());
	_profile.resize (A.coldim ());
	_indices.resize (A.rowdim ());

	for (i = _indices.begin (), idx = 0; i != _indices.end (); ++i, ++idx)
		*i = idx;

	setIN (_U);
	_P.clear ();
	_MD.copy (_A, A);
	_profile_idx = 0;
	kthGaussJordan (rank, det, 0, 0, A.coldim (), _one);

	buildMinimalPermutation (P, rank, A.rowdim (), _P);
	buildMinimalPermutationFromProfile (Q, rank, A.coldim (), _profile);

	_MD.permuteColumns (_U, _P.rbegin (), _P.rend ());
	_MD.permuteColumns (_U, P.begin (), P.end ());

	DenseSubmatrix<Element> Tu1 (Tu, rank, 0, A.rowdim () - rank, rank);
	DenseSubmatrix<Element> U2 (_U, rank, 0, A.rowdim () - rank, rank);
	_MD.copy (Tu1, U2);

	DenseSubmatrix<Element> Ainv1 (U, 0, 0, rank, rank);
	DenseSubmatrix<Element> U1 (_U, 0, 0, rank, rank);
	_F.inv (dinv, det);
	_MD.mul (Ainv1, U1, dinv);

	DenseSubmatrix<Element> Tv1 (Tv, 0, rank, rank, A.coldim () - rank);
	DenseSubmatrix<Element> U3 (_U, 0, 0, rank, A.rowdim ());
	DenseSubmatrix<Element> A2 (_A, 0, rank, A.rowdim (), A.coldim () - rank);
	_MD.copy (_A, A);
	_MD.permuteColumns (_A, Q.begin (), Q.end ());
	_MD.permuteRows (A2, P.begin (), P.end ());
	_MD.mul (Tv1, U3, A2);
	_MD.negin (Tv1);

	profile.resize (rank);
	std::copy (_profile.begin (), _profile.begin () + rank, profile.begin ());

	return U;
}

template <class Field, class Matrix>
std::ostream &Eliminator<Field, Matrix>::writeFilter (std::ostream &out, const std::vector<bool> &v) const
{
	std::vector<bool>::const_iterator i;

	for (i = v.begin (); i != v.end (); ++i) {
		if (*i)
			out << "1 ";
		else
			out << "0 ";
	}

	return out;
}

template <class Field, class Matrix>
std::ostream &Eliminator<Field, Matrix>::writePermutation (std::ostream &out, const Permutation &P) const
{
	Permutation::const_iterator i;

	for (i = P.begin (); i != P.end (); ++i)
		out << '(' << i->first << ' ' << i->second << ')';

	return out;
}

/* Perform the kth indexed Gauss-Jordan transform on _A, storing the
 * transformation matrix in _U and the permutation in _P. The caller is
 * responsible for ensuring that _U and _P are the identity and that _A is set
 * to a copy of the input on the initial call.
 */

template <class Field, class Matrix>
Matrix &Eliminator<Field, Matrix>::kthGaussJordan
	(unsigned int                  &r,
	 typename Field::Element       &d,
	 unsigned int                   k,
	 unsigned int                   s,
	 unsigned int                   m,
	 const typename Field::Element &d0)
{
	unsigned int i;

#ifdef ELIM_DETAILED_TRACE
	commentator.start ("kth indexed Gauss-Jordan transform", "Eliminator::kthGaussJordan");

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "k = " << k << std::endl;
	report << "Starting column: " << s << std::endl;
	report << "Column dimension: " << m << std::endl;

	DenseMatrixBase<Element> Acopy (_A);

	unsigned int P_start = _P.size ();
#endif

	DenseSubmatrix<Element> Ap (_A, k, s, _A.rowdim () - k, m);

	if (_MD.isZero (Ap)) {
		r = 0;
		_F.assign (d, d0);
	}
	else if (m == 1) {
		// Find minimal index i > k with _A[i, 1] != 0
		for (i = 0; i < _A.rowdim (); ++i)
			if (_indices[i] >= k && !_F.isZero (_A.getEntry (_indices[i], s)))
				break;

		linbox_check (i < _A.rowdim ());

		_T[s] = true;  // This column is independent

		if (_indices[i] != k) _P.push_back (Transposition (_indices[i], k));

		r = 1;
		_A.getEntry (d, _indices[i], s);

		typename Matrix::ColIterator Uk = _U.colBegin () + k;
		typename Matrix::ColIterator A1 = _A.colBegin () + s;

		_VD.neg (*Uk, *A1);
		_F.assign ((*Uk)[_indices[i]], (*Uk)[k]);
		_F.assign ((*Uk)[k], d0);

		std::swap (_indices[i], _indices[k]);

		for (i = k + 1; i < _U.rowdim (); ++i)
			_U.setEntry (i, i, d);

		_profile[_profile_idx++] = s;
	}
	else {
		unsigned int m1 = m / 2;
		unsigned int m2 = m - m1;
		unsigned int r1, r2;
		typename Field::Element d1, d0inv, d1inv, d1neg;

		DenseSubmatrix<Element> B (_A, 0, s + m1, _A.rowdim (), m2);

		unsigned int P_start = _P.size ();

		kthGaussJordan (r1, d1, k, s, m1, d0);

		unsigned int P_end = _P.size ();

		_MD.permuteRows (B, _P.begin () + P_start, _P.end ());

		unsigned int l1 = _U.rowdim () - (k + r1);

		DenseSubmatrix<Element> a (_U,    0,      k,      k,  r1);
		DenseSubmatrix<Element> u (_U,    k,      k,      r1, r1);
		DenseSubmatrix<Element> c (_U,    k + r1, k,      l1, r1);

		DenseSubmatrix<Element> et (_tmp, 0,      0,      k,  m2);
		DenseSubmatrix<Element> gt (_tmp, 0,      0,      l1, m2);

		DenseSubmatrix<Element> e (_A,    0,      s + m1, k,  m2);
		DenseSubmatrix<Element> f (_A,    k,      s + m1, r1, m2);
		DenseSubmatrix<Element> g (_A,    k + r1, s + m1, l1, m2);

		_F.inv (d0inv, d0);

		_MD.mul (et, a, f);
		_MD.mulin (e, d1);
		_MD.addin (e, et);
		_MD.mulin (e, d0inv);

		_MD.mul (gt, c, f);
		_MD.mulin (g, d1);
		_MD.addin (g, gt);
		_MD.mulin (g, d0inv);

		_MD.leftMulin (u, f);
		_MD.mulin (f, d0inv);

#ifdef ELIM_DETAILED_TRACE
		report << "(" << k << ") Matrix A prepared for second recursive call: " << std::endl;
		_MD.write (report, _A);
#endif

		kthGaussJordan (r2, d, k + r1, s + m1, m2, d1);

#ifdef ELIM_DETAILED_TRACE
		report << "(" << k << ") Transform U after recursive calls: " << std::endl;
		_MD.write (report, _U);
#endif

		DenseSubmatrix<Element> U1 (_U, 0, k, _U.rowdim (), r1);

		_F.neg (d1neg, d1);
		adddIN (_U, d1neg);
		_MD.permuteRows (U1, _P.begin () + P_end, _P.end ());
		adddIN (_U, d1);

#ifdef ELIM_DETAILED_TRACE
		report << "(" << k << ") P2 U P2^-1: " << std::endl;
		_MD.write (report, _U);
#endif

		r = r1 + r2;

		unsigned int l2 = _U.rowdim () - (k + r);

		DenseSubmatrix<Element> a1    (_U, 0,      k,      k,  r1);
		DenseSubmatrix<Element> u1    (_U, k,      k,      r1, r1);
		DenseSubmatrix<Element> c1    (_U, k + r1, k,      r2, r1);
		DenseSubmatrix<Element> c1bar (_U, k + r,  k,      l2, r1);

		DenseSubmatrix<Element> &a11 = a1;
		DenseSubmatrix<Element> &u11 = u1;
		DenseSubmatrix<Element> &u21 = c1;
		DenseSubmatrix<Element> &c11 = c1bar;

		DenseSubmatrix<Element> a11t  (_tmp, 0,    0,      k,  r1);
		DenseSubmatrix<Element> u11t  (_tmp, 0,    0,      r1, r1);
		DenseSubmatrix<Element> c11t  (_tmp, 0,    0,      l2, r1);

		DenseSubmatrix<Element> a2    (_U, 0,      k + r1, k,  r2);
		DenseSubmatrix<Element> a2bar (_U, k,      k + r1, r1, r2);
		DenseSubmatrix<Element> u2    (_U, k + r1, k + r1, r2, r2);
		DenseSubmatrix<Element> c2    (_U, k + r,  k + r1, l2, r2);

		_F.inv (d1inv, d1);

		_MD.mul (a11t, a2, c1);
		_MD.mulin (a1, d);
		_MD.addin (a11, a11t);   // a11 <- d * a1 + a2 * c1
		_MD.mulin (a11, d1inv);

		_MD.mul (u11t, a2bar, c1);
		_MD.mulin (u1, d);
		_MD.addin (u11, u11t);   // u11 <- d * u1 + a2bar * c1
		_MD.mulin (u11, d1inv);

		_MD.mul (c11t, c2, c1);
		_MD.mulin (c1bar, d);
		_MD.addin (c11, c11t);   // c11 <- d * c1bar + c2 * c1
		_MD.mulin (c11, d1inv);

		_MD.leftMulin (u2, c1);  // u21 <- u2 * c1
		_MD.mulin (u21, d1inv);
	}

#ifdef ELIM_DETAILED_TRACE
	report << "(" << k << ") Finished U: " << std::endl;
	_MD.write (report, _U);

	report << "(" << k << ") Finished P: " << std::endl;
	reportPermutation (report, _P);

	typename Field::Element dinv, d0inv;

	_F.inv (dinv, d);
	_F.inv (d0inv, d0);

	DenseMatrixBase<Element> R (_A.rowdim () - k, _A.coldim () - s);
	DenseSubmatrix<Element> Atest (Acopy, k, s, _A.rowdim () - k, _A.coldim () - s);
	DenseSubmatrix<Element> Utest (_U, k, k, _U.rowdim () - k, _U.coldim () - k);
	_MD.permuteRows (Acopy, _P.begin () + P_start, _P.end ());

	report << "(" << k << ") PA: " << std::endl;
	_MD.write (report, Acopy);

	_MD.mul (R, Utest, Atest);
	_MD.mulin (R, dinv);
	_MD.mulin (R, d0inv);

	report << "(" << k << ") R:=1/d U 1/d0 PA: " << std::endl;
	_MD.write (report, R);

	commentator.stop ("done", NULL, "Eliminator::kthGaussJordan");
#endif

	return _U;
}

template <class Field, class Matrix>
template <class Matrix1>
Matrix1 &Eliminator<Field, Matrix>::adddIN
	(Matrix1                       &A,
	 const typename Field::Element &d) const
{
	typename Matrix1::RowIterator i;
	unsigned int idx;

	for (i = A.rowBegin (), idx = 0; i != A.rowEnd (); ++i, ++idx)
		_F.addin ((*i)[idx], d);

	return A;
}

template <class Field, class Matrix>
template <class Matrix1>
Matrix1 &Eliminator<Field, Matrix>::setIN (Matrix1 &A) const 
{
	linbox_check (A.coldim () == A.rowdim ());

	typename Matrix1::RowIterator i;
	size_t i_idx;

	for (i = A.rowBegin (), i_idx = 0; i != A.rowEnd (); ++i, ++i_idx) {
		_VD.subin (*i, *i);
		_F.assign ((*i)[i_idx], _one);
	}

	return A;
}

/* Clean out the given priority index list and add new elements as needed */

template <class Field, class Matrix>
void Eliminator<Field, Matrix>::cleanPriorityIndexList
	(std::list<unsigned int> &list,
	 std::vector<bool>       &S,
	 std::vector<bool>       &old_S) const
{
	std::list<unsigned int>::iterator li;
	std::vector<bool>::iterator si, old_si;
	unsigned int idx;

	for (li = list.begin (); li != list.end ();) {
		if (S[*li])
			li = list.erase (li);
		else
			++li;
	}

	for (si = S.begin (), old_si = old_S.begin (), idx = 0; si != S.end (); ++si, ++old_si, ++idx) {
		if (!*si && *old_si)
			list.push_back (idx);
	}
}

/* Permute the entries given bit vector using the given permutation */

template <class Field, class Matrix>
template <class Iterator>
std::vector<bool> &Eliminator<Field, Matrix>::permute (std::vector<bool> &v, Iterator P_start, Iterator P_end) const
{
	Iterator i;

	for (i = P_start; i != P_end; ++i) {
		bool tmp = v[i->first];
		v[i->first] = v[i->second];
		v[i->second] = tmp;
	}

	return v;
}

template <class Field, class Matrix>
typename Eliminator<Field, Matrix>::Permutation &
Eliminator<Field, Matrix>::buildPermutation (Permutation &P, const std::list<unsigned int> &pidx) const
{
	unsigned int offset, current;
	Permutation::iterator i;
	std::list<unsigned int>::const_iterator li;

	P.clear ();

	for (li = pidx.begin (), offset = 0; li != pidx.end (); ++li, ++offset) {
		if (*li > offset)
			P.push_back (Transposition (*li, offset));
		else if (*li < offset) {
			// We need to figure out to what place the original
			// bubbled. We'll use a naive algorithm here, since I
			// don't anticipate this being a problem too much, and
			// it's O(n) in any case.

			current = *li;

			for (i = P.begin (); i != P.end (); ++i)
				if (i->second == current)
					current = i->first;

			if (current != offset)
				P.push_back (Transposition (current, offset));
		}
	}

	return P;
}

template <class Field, class Matrix>
typename Eliminator<Field, Matrix>::Permutation &
Eliminator<Field, Matrix>::buildMinimalPermutation (Permutation &P, unsigned int rank,
						    unsigned int dim, const Permutation &Pold) 
{
	Permutation::const_reverse_iterator j;
	unsigned int idx, idx2;

	P.clear ();

	std::fill (_S.begin (), _S.begin () + rank, true);
	std::fill (_S.begin () + rank, _S.begin () + dim, false);

	for (j = Pold.rbegin (); j != Pold.rend (); ++j) {
		bool tmp = _S[j->first];
		_S[j->first] = _S[j->second];
		_S[j->second] = tmp;
	}

	idx = 0;
	idx2 = dim - 1;

	while (idx < rank && idx2 >= rank) {
		while (_S[idx] && idx < rank) ++idx;
		while (!_S[idx2] && idx2 >= rank) --idx2;

		if (idx < rank && idx2 >= rank)
			P.push_back (Transposition (idx, idx2));

		++idx;
		--idx2;
	}

	return P;
}

template <class Field, class Matrix>
typename Eliminator<Field, Matrix>::Permutation &
Eliminator<Field, Matrix>::buildMinimalPermutationFromProfile (Permutation &P, unsigned int rank,
							       unsigned int dim, const std::vector<unsigned int> &profile) 
{
	typename std::vector<unsigned int>::const_iterator j;
	unsigned int idx = 0;

	P.clear ();

	for (j = profile.begin (); j != profile.begin () + rank; ++j, ++idx)
		if (*j != idx)
			P.push_back (Transposition (*j, idx));

	return P;
}

} // namespace LinBox

#endif // __LINBOX_eliminator_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
