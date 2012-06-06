/* eliminator.inl
 * Copyright (C) 2002, 2003 LinBox, Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.waterloo.ca>
 *
 * --------------------------------------------
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========

 * Elimination code for lookahead block Lanczos
 */

#ifndef __LINBOX_eliminator_INL
#define __LINBOX_eliminator_INL

#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/util/debug.h"
#include "linbox/solutions/methods.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"

#include "eliminator.h"


namespace LinBox
{

	std::ostream &reportPermutation (std::ostream &out,
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
	Eliminator<Field, Matrix>::Eliminator (const Field &F, unsigned int N) :
		_field (F), _VD (F), _MD (F), _number (N), _indepRows (N), _indepCols (N)
	{
		_field.init (_one, 1);
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

		_matA.resize (A.rowdim (), A.coldim ());
		_matU.resize (A.rowdim (), A.rowdim ());
		_tmp.resize (A.rowdim (), A.coldim ());
		_profile.resize (A.coldim ());
		_indices.resize (A.rowdim ());

		for (i = _indices.begin (), idx = 0; i != _indices.end (); ++i, ++idx)
			*i = idx;

		_MD.subin (_matA, _matA);

		BlasMatrix<Field> A1 (_matA, 0, 0, A.rowdim (), A.coldim ());
		_MD.copy (A1, A);

		setIN (_matU);
		_perm.clear ();
		_profile_idx = 0;

		kthGaussJordan (rank, d, 0, 0, A.coldim (), _one);

		buildMinimalPermutation (P, rank, A.rowdim (), _perm);
		buildMinimalPermutationFromProfile (Q, rank, A.coldim (), _profile);

		_MD.permuteColumns (_matU, _perm.rbegin (), _perm.rend ());
		_MD.permuteColumns (_matU, P.begin (), P.end ());

		BlasMatrix<Field> Tu1 (Tu, rank, 0, A.rowdim () - rank, rank);
		BlasMatrix<Field> U2 (_matU, rank, 0, A.rowdim () - rank, rank);
		_MD.copy (Tu1, U2);

		BlasMatrix<Field> Ainv1 (Ainv, 0, 0, rank, rank);
		BlasMatrix<Field> U1 (_matU, 0, 0, rank, rank);
		_field.inv (dinv, d);
		_MD.mul (Ainv1, U1, dinv);

		BlasMatrix<Field> Tv1 (Tv, 0, rank, rank, A.coldim () - rank);
		BlasMatrix<Field> U3 (_matU, 0, 0, rank, A.rowdim ());
		BlasMatrix<Field> A2 (_matA, 0, rank, A.rowdim (), A.coldim () - rank);
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

		linbox_check (_number == S.size ());
		linbox_check (_number == T.size ());
		linbox_check (_number == W.rowdim ());
		linbox_check (_number == W.coldim ());
		linbox_check (_number == A.rowdim ());
		linbox_check (_number == A.coldim ());

#ifdef ELIM_DETAILED_TRACE
		commentator().start ("Computing W, S, T", "Eliminator::permuteAndInvert", _number);

		std::ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input matrix:" << std::endl;
		_MD.write (report, A);
#endif

		/* Apply initial permutations to A, copying to _matA */

		buildPermutation (Qp, rightPriorityIdx); // Column permutation

#ifdef ELIM_DETAILED_TRACE
		report << "Column permutation: ";
		reportPermutation (report, Qp) << endl;
#endif

		_matA.resize (A.rowdim (), A.coldim ());
		_matU.resize (A.rowdim (), A.rowdim ());
		_tmp.resize (A.rowdim (), A.coldim ());
		_indices.resize (A.rowdim ());

		for (i = _indices.begin (), idx = 0; i != _indices.end (); ++i, ++idx)
			*i = idx;

		_MD.copy (_matA, A);

		_MD.permuteColumns (_matA, Qp.begin (), Qp.end ());

		/* Initialise temporaries for the computation */

		setIN (_matU);
		_perm.clear ();
		_profile.resize (A.coldim ());
		_profile_idx = 0;
		std::fill (_indepCols.begin (), _indepCols.end (), false);

		/* Run the computation */

		kthGaussJordan (rank, d, 0, 0, _matA.coldim (), _one);

		/* Set _indepRows based on the permutation */
		std::fill (_indepRows.begin (), _indepRows.begin () + rank, true);
		std::fill (_indepRows.begin () + rank, _indepRows.end (), false);
		permute (_indepRows, _perm.rbegin (), _perm.rend ());

		permute (_indepCols, Qp.rbegin (), Qp.rend ());

		/* Apply final permutations to _matU, copying to W */
		BlasMatrix<Field> U1 (_matU, rank, 0, _matU.rowdim () - rank, _matU.coldim ());
		_MD.subin (U1, U1);

		_MD.permuteColumns (_matU, _perm.rbegin (), _perm.rend ());

		/* Divide _matU by the determinant and copy to W */
		_field.invin (d);
		_MD.mulin (_matU, d);

		_MD.subin (W, W);

		typename std::vector<unsigned int>::iterator pi;

		for (pi = _profile.begin (), ui = _matU.rowBegin (); pi != _profile.begin () + rank; ++ui, ++pi)
			_VD.copy (*(W.rowBegin () + *pi), *ui);

		//	_MD.permuteRows (W, Qp.rbegin (), Qp.rend ());

		/* Rebuild leftPriorityIdx and rightPriorityIdx */
		cleanPriorityIndexList (rightPriorityIdx, _indepCols, T);

		/* Reverse the row priority index list */
		std::reverse (_indices.begin (), _indices.end ());

		S = _indepRows;
		T = _indepCols;

#ifdef ELIM_DETAILED_TRACE
		report << "Computed W:" << std::endl;
		_MD.write (report, W);

		commentator().stop ("done", NULL, "Eliminator::permuteAndInvert");
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

		_matA.resize (A.rowdim (), A.coldim ());
		_matU.resize (A.rowdim (), A.rowdim ());
		_tmp.resize (A.rowdim (), A.coldim ());
		_profile.resize (A.coldim ());
		_indices.resize (A.rowdim ());

		for (i = _indices.begin (), idx = 0; i != _indices.end (); ++i, ++idx)
			*i = idx;

		setIN (_matU);
		_perm.clear ();
		_MD.copy (_matA, A);
		_profile_idx = 0;
		kthGaussJordan (rank, det, 0, 0,  (unsigned int) A.coldim (), _one);

		buildMinimalPermutation (P, rank,  (unsigned int) A.rowdim (), _perm);
		buildMinimalPermutationFromProfile (Q, rank,  (unsigned int) A.coldim (), _profile);

		_MD.permuteColumns (_matU, _perm.rbegin (), _perm.rend ());
		_MD.permuteColumns (_matU, P.begin (), P.end ());

		BlasMatrix<Field> Tu1 (Tu, rank, 0, A.rowdim () - rank, rank);
		BlasMatrix<Field> U2 (_matU, rank, 0, A.rowdim () - rank, rank);
		_MD.copy (Tu1, U2);

		BlasMatrix<Field> Ainv1 (U, 0, 0, rank, rank);
		BlasMatrix<Field> U1 (_matU, 0, 0, rank, rank);
		_field.inv (dinv, det);
		_MD.mul (Ainv1, U1, dinv);

		BlasMatrix<Field> Tv1 (Tv, 0, rank, rank, A.coldim () - rank);
		BlasMatrix<Field> U3 (_matU, 0, 0, rank, A.rowdim ());
		BlasMatrix<Field> A2 (_matA, 0, rank, A.rowdim (), A.coldim () - rank);
		_MD.copy (_matA, A);
		_MD.permuteColumns (_matA, Q.begin (), Q.end ());
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

	/* Perform the kth indexed Gauss-Jordan transform on _matA, storing the
	 * transformation matrix in _matU and the permutation in _perm. The caller is
	 * responsible for ensuring that _matU and _perm are the identity and that _matA is set
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
		commentator().start ("kth indexed Gauss-Jordan transform", "Eliminator::kthGaussJordan");

		std::ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "k = " << k << std::endl;
		report << "Starting column: " << s << std::endl;
		report << "Column dimension: " << m << std::endl;

		BlasMatrix<Field> Acopy (_matA);

		unsigned int P_start = _perm.size ();
#endif

		BlasMatrix<Field> Ap (_matA, k, s, _matA.rowdim () - k, m);

		if (_MD.isZero (Ap)) {
			r = 0;
			_field.assign (d, d0);
		}
		else if (m == 1) {
			// Find minimal index i > k with _matA[i, 1] != 0
			for (i = 0; i < _matA.rowdim (); ++i)
				if (_indices[i] >= k && !_field.isZero (_matA.getEntry (_indices[i], s)))
					break;

			linbox_check (i < _matA.rowdim ());

			_indepCols[s] = true;  // This column is independent

			if (_indices[i] != k) _perm.push_back (Transposition (_indices[i], k));

			r = 1;
			_matA.getEntry (d, _indices[i], s);

			typename Matrix::ColIterator Uk = _matU.colBegin () + k;
			typename Matrix::ColIterator A1 = _matA.colBegin () + s;

			_VD.neg (*Uk, *A1);
			_field.assign ((*Uk)[_indices[i]], (*Uk)[k]);
			_field.assign ((*Uk)[k], d0);

			std::swap (_indices[i], _indices[k]);

			for (i = k + 1; i < _matU.rowdim (); ++i)
				_matU.setEntry (i, i, d);

			_profile[_profile_idx++] = s;
		}
		else {
			unsigned int m1 = m / 2;
			unsigned int m2 = m - m1;
			unsigned int r1, r2;
			typename Field::Element d1, d0inv, d1inv, d1neg;

			BlasMatrix<Field> B (_matA, 0, s + m1, _matA.rowdim (), m2);

			unsigned int P_start = (unsigned int) _perm.size ();

			kthGaussJordan (r1, d1, k, s, m1, d0);

			unsigned int P_end =  (unsigned int) _perm.size ();

			_MD.permuteRows (B, _perm.begin () + P_start, _perm.end ());

			unsigned int l1 = (unsigned int) _matU.rowdim () - (k + r1);

			BlasMatrix<Field> a (_matU,    0,      k,      k,  r1);
			BlasMatrix<Field> u (_matU,    k,      k,      r1, r1);
			BlasMatrix<Field> c (_matU,    k + r1, k,      l1, r1);

			BlasMatrix<Field> et (_tmp, 0,      0,      k,  m2);
			BlasMatrix<Field> gt (_tmp, 0,      0,      l1, m2);

			BlasMatrix<Field> e (_matA,    0,      s + m1, k,  m2);
			BlasMatrix<Field> f (_matA,    k,      s + m1, r1, m2);
			BlasMatrix<Field> g (_matA,    k + r1, s + m1, l1, m2);

			_field.inv (d0inv, d0);

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
			_MD.write (report, _matA);
#endif

			kthGaussJordan (r2, d, k + r1, s + m1, m2, d1);

#ifdef ELIM_DETAILED_TRACE
			report << "(" << k << ") Transform U after recursive calls: " << std::endl;
			_MD.write (report, _matU);
#endif

			BlasMatrix<Field> U1 (_matU, 0, k, _matU.rowdim (), r1);

			_field.neg (d1neg, d1);
			adddIN (_matU, d1neg);
			_MD.permuteRows (U1, _perm.begin () + P_end, _perm.end ());
			adddIN (_matU, d1);

#ifdef ELIM_DETAILED_TRACE
			report << "(" << k << ") P2 U P2^-1: " << std::endl;
			_MD.write (report, _matU);
#endif

			r = r1 + r2;

			unsigned int l2 = (unsigned int) _matU.rowdim () - (k + r);

			BlasMatrix<Field> a1    (_matU, 0,      k,      k,  r1);
			BlasMatrix<Field> u1    (_matU, k,      k,      r1, r1);
			BlasMatrix<Field> c1    (_matU, k + r1, k,      r2, r1);
			BlasMatrix<Field> c1bar (_matU, k + r,  k,      l2, r1);

			BlasMatrix<Field> &a11 = a1;
			BlasMatrix<Field> &u11 = u1;
			BlasMatrix<Field> &u21 = c1;
			BlasMatrix<Field> &c11 = c1bar;

			BlasMatrix<Field> a11t  (_tmp, 0,    0,      k,  r1);
			BlasMatrix<Field> u11t  (_tmp, 0,    0,      r1, r1);
			BlasMatrix<Field> c11t  (_tmp, 0,    0,      l2, r1);

			BlasMatrix<Field> a2    (_matU, 0,      k + r1, k,  r2);
			BlasMatrix<Field> a2bar (_matU, k,      k + r1, r1, r2);
			BlasMatrix<Field> u2    (_matU, k + r1, k + r1, r2, r2);
			BlasMatrix<Field> c2    (_matU, k + r,  k + r1, l2, r2);

			_field.inv (d1inv, d1);

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
		_MD.write (report, _matU);

		report << "(" << k << ") Finished P: " << std::endl;
		reportPermutation (report, _perm);

		typename Field::Element dinv, d0inv;

		_field.inv (dinv, d);
		_field.inv (d0inv, d0);

		BlasMatrix<Field> R (_matA.rowdim () - k, _matA.coldim () - s);
		BlasMatrix<Field> Atest (Acopy, k, s, _matA.rowdim () - k, _matA.coldim () - s);
		BlasMatrix<Field> Utest (_matU, k, k, _matU.rowdim () - k, _matU.coldim () - k);
		_MD.permuteRows (Acopy, _perm.begin () + P_start, _perm.end ());

		report << "(" << k << ") PA: " << std::endl;
		_MD.write (report, Acopy);

		_MD.mul (R, Utest, Atest);
		_MD.mulin (R, dinv);
		_MD.mulin (R, d0inv);

		report << "(" << k << ") R:=1/d U 1/d0 PA: " << std::endl;
		_MD.write (report, R);

		commentator().stop ("done", NULL, "Eliminator::kthGaussJordan");
#endif

		return _matU;
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
			_field.addin ((*i)[idx], d);

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
			_field.assign ((*i)[i_idx], _one);
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

		std::fill (_indepRows.begin (), _indepRows.begin () + rank, true);
		std::fill (_indepRows.begin () + rank, _indepRows.begin () + dim, false);

		for (j = Pold.rbegin (); j != Pold.rend (); ++j) {
			bool tmp = _indepRows[j->first];
			_indepRows[j->first] = _indepRows[j->second];
			_indepRows[j->second] = tmp;
		}

		idx = 0;
		idx2 = dim - 1;

		while (idx < rank && idx2 >= rank) {
			while (_indepRows[idx] && idx < rank) ++idx;
			while (!_indepRows[idx2] && idx2 >= rank) --idx2;

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


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

