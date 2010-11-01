
/* linbox/matrix/matrix-domain.inl
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_domain_INL
#define __LINBOX_matrix_domain_INL

#include "linbox/matrix/transpose-matrix.h"

namespace LinBox
{

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::copyRow (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.copy (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::copyCol (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.copy (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
bool MatrixDomain<Field>::areEqualRow (const Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		if (!_VD.areEqual (*i, *j))
			return false;

	return true;
}

template <class Field>
template <class Matrix1, class Matrix2>
bool MatrixDomain<Field>::areEqualCol (const Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		if (!_VD.areEqual (*i, *j))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
bool MatrixDomain<Field>::isZeroRow (const Matrix &A) const
{
	typename Matrix::ConstRowIterator i;

	i = A.rowBegin ();

	for (; i != A.rowEnd (); ++i)
		if (!_VD.isZero (*i))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
bool MatrixDomain<Field>::isZeroCol (const Matrix &A) const
{
	typename Matrix::ConstColIterator i;

	i = A.colBegin ();

	for (; i != A.colEnd (); ++i)
		if (!_VD.isZero (*i))
			return false;

	return true;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::addRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (A.coldim () == C.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;
	typename Matrix3::ConstRowIterator k;

	i = C.rowBegin ();
	j = A.rowBegin ();
	k = B.rowBegin ();

	for (; i != C.rowEnd (); ++i, ++j, ++k)
		_VD.add (*i, *j, *k);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::addCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (A.coldim () == C.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::ConstColIterator k;

	i = C.colBegin ();
	j = A.colBegin ();
	k = B.colBegin ();

	for (; i != C.colEnd (); ++i, ++j, ++k)
		_VD.add (*i, *j, *k);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::addinRow (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.addin (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::addinCol (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.addin (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::subRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (A.coldim () == C.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;
	typename Matrix3::ConstRowIterator k;

	i = C.rowBegin ();
	j = A.rowBegin ();
	k = B.rowBegin ();

	for (; i != C.rowEnd (); ++i, ++j, ++k)
		_VD.sub (*i, *j, *k);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::subCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (A.coldim () == C.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::ConstColIterator k;

	i = C.colBegin ();
	j = A.colBegin ();
	k = B.colBegin ();

	for (; i != C.colEnd (); ++i, ++j, ++k)
		_VD.sub (*i, *j, *k);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::subinRow (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.subin (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::subinCol (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.subin (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::negRow (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.neg (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::negCol (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.neg (*i, *j);

	return A;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomain<Field>::neginRow (Matrix &A) const
{
	typename Matrix::RowIterator i;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i)
		_VD.negin (*i);

	return A;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomain<Field>::neginCol (Matrix &A) const
{
	typename Matrix::ColIterator i;

	for (i = A.colBegin (); i != A.colEnd (); ++i)
		_VD.negin (*i);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::mulRowRowCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstRowIterator i;
	typename Matrix3::ConstColIterator j;
	typename Matrix1::RowIterator l1;
	typename Matrix1::Row::iterator l2;

	for (i = A.rowBegin (), l1 = C.rowBegin (); i != A.rowEnd (); ++i, ++l1)
		for (j = B.colBegin (), l2 = l1->begin (); j != B.colEnd (); ++j, ++l2)
			_VD.dot (*l2, *i, *j);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::mulColRowCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstRowIterator i;
	typename Matrix3::ConstColIterator j;
	typename Matrix1::ColIterator l1;
	typename Matrix1::Col::iterator l2;

	for (j = B.colBegin (), l1 = C.colBegin (); j != B.colEnd (); ++j, ++l1)
		for (i = A.rowBegin (), l2 = l1->begin (); i != A.rowEnd (); ++i, ++l2)
			_VD.dot (*l2, *i, *j);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::mulRowRowRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstRowIterator i = A.rowBegin ();
	typename Matrix1::RowIterator j = C.rowBegin ();

	TransposeMatrix<const Matrix3> BT (B);

	for (; i != A.rowEnd (); ++i, ++j)
		vectorMul (*j, BT, *i);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::mulColColCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix3::ConstColIterator i = B.colBegin ();
	typename Matrix1::ColIterator j = C.colBegin ();

	for (; i != B.colEnd (); ++i, ++j)
		vectorMul (*j, A, *i);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomain<Field>::leftMulin (const Matrix1 &A, Matrix2 &B) const
{
	linbox_check (A.rowdim () == A.coldim ());
	linbox_check (A.coldim () == B.rowdim ());

	typename LinBox::Vector<Field>::Dense t (A.rowdim ());

	typename Matrix2::ColIterator i;
	typename Matrix1::ConstRowIterator j;

	typename LinBox::Vector<Field>::Dense::iterator k;

	for (i = B.colBegin (); i != B.colEnd (); ++i) {
		for (j = A.rowBegin (), k = t.begin (); j != A.rowEnd (); ++j, ++k)
			_VD.dot (*k, *j, *i);

		_VD.copy (*i, t);
	}

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::rightMulin (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (B.rowdim () == B.coldim ());

	typename LinBox::Vector<Field>::Dense t (B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstColIterator j;

	typename LinBox::Vector<Field>::Dense::iterator k;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i) {
		for (j = B.colBegin (), k = t.begin (); j != B.colEnd (); ++j, ++k)
			_VD.dot (*k, *i, *j);

		_VD.copy (*i, t);
	}

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::mulRow (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a) const
{
	linbox_check (C.rowdim () == B.rowdim ());
	linbox_check (C.coldim () == B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = C.rowBegin ();
	j = B.rowBegin ();

	for (; i != C.rowEnd (); ++i, ++j)
		_VD.mul (*i, *j, a);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::mulCol (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a) const
{
	linbox_check (C.rowdim () == B.rowdim ());
	linbox_check (C.coldim () == B.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;

	i = C.colBegin ();
	j = B.colBegin ();

	for (; i != C.colEnd (); ++i, ++j)
		_VD.mul (*i, *j, a);

	return C;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomain<Field>::mulinRow (Matrix &B, const typename Field::Element &a) const
{
	typename Matrix::RowIterator i;

	for (i = B.rowBegin (); i != B.rowEnd (); ++i)
		_VD.mulin (*i, a);

	return B;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomain<Field>::mulinCol (Matrix &B, const typename Field::Element &a) const
{
	typename Matrix::ColIterator i;

	for (i = B.colBegin (); i != B.colEnd (); ++i)
		_VD.mulin (*i, a);

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::axpyinRowRowCol (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const
{
	linbox_check (A.coldim () == X.rowdim ());
	linbox_check (A.rowdim () == Y.rowdim ());
	linbox_check (X.coldim () == Y.coldim ());

	typename Matrix2::ConstRowIterator i;
	typename Matrix3::ConstColIterator j;
	typename Matrix1::RowIterator l1;
	typename Matrix1::Row::iterator l2;

	typename Field::Element t;

	for (i = A.rowBegin (), l1 = Y.rowBegin (); i != A.rowEnd (); ++i, ++l1) {
		for (j = X.colBegin (), l2 = l1->begin (); j != X.colEnd (); ++j, ++l2) {
			_VD.dot (t, *i, *j);
			_F.addin (*l2, t);
		}
	}

	return Y;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::axpyinColRowCol (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const
{
	linbox_check (A.coldim () == X.rowdim ());
	linbox_check (A.rowdim () == Y.rowdim ());
	linbox_check (X.coldim () == Y.coldim ());

	typename Matrix2::ConstRowIterator i;
	typename Matrix3::ConstColIterator j;
	typename Matrix1::ColIterator l1;
	typename Matrix1::Col::iterator l2;

	typename Field::Element t;

	for (j = X.colBegin (), l1 = Y.colBegin (); j != X.colEnd (); ++j, ++l1) {
		for (i = A.rowBegin (), l2 = l1->begin (); i != A.rowEnd (); ++i, ++l2) {
			_VD.dot (t, *i, *j);
			_F.addin (*l2, t);
		}
	}

	return Y;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::axpyinRowRowRow (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const
{
	linbox_check (A.coldim () == X.rowdim ());
	linbox_check (A.rowdim () == Y.rowdim ());
	linbox_check (X.coldim () == Y.coldim ());

	typename LinBox::Vector<Field>::Dense t (X.coldim ());

	typename Matrix2::ConstRowIterator i = A.rowBegin ();
	typename Matrix1::RowIterator j = Y.rowBegin ();

	TransposeMatrix<const Matrix3> XT (X);

	for (; i != A.rowEnd (); ++i, ++j) {
		vectorMul (t, XT, *i);
		_VD.addin (*j, t);
	}

	return Y;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &MatrixDomain<Field>::axpyinColColCol (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const
{
	linbox_check (A.coldim () == X.rowdim ());
	linbox_check (A.rowdim () == Y.rowdim ());
	linbox_check (X.coldim () == Y.coldim ());

	typename LinBox::Vector<Field>::Dense t (A.rowdim ());

	typename Matrix3::ConstColIterator i = X.colBegin ();
	typename Matrix1::ColIterator j = Y.colBegin ();

	for (; i != X.colEnd (); ++i, ++j) {
		vectorMul (t, A, *i);
		_VD.addin (*j, t);
	}

	return Y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector1::iterator j = w.begin ();

            // JGD 02.09.2008 : when sizes differ
            // A must decide if dot is possible, not w 
// 	for (; j != w.end (); ++j, ++i)
// 		_VD.dot (*j, v, *i);
	for (; i != A.rowEnd (); ++j, ++i)
		_VD.dot (*j, v, *i);

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						 VectorCategories::SparseSequenceVectorTag) const
{
	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	w.clear ();

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, v, *i);

		if (!_F.isZero (t))
			w.push_back (std::pair<size_t, typename Field::Element> (idx, t));
	}

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						 VectorCategories::SparseAssociativeVectorTag) const
{
	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	w.clear ();

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, v, *i);

		if (!_F.isZero (t))
			w[idx] = t;
	}

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						 VectorCategories::SparseParallelVectorTag) const
{
	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	w.first.clear ();
	w.second.clear ();

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, v, *i);

		if (!_F.isZero (t)) {
			w.first.push_back (idx);
			w.second.push_back (t);
		}
	}

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Field>::mulColDense
	(const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j = v.begin ();

	VD.subin (w, w);

	for (; j != v.end (); ++j, ++i)
		VD.axpyin (w, *j, *i);

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						 VectorCategories::DenseVectorTag,
						 VectorCategories::SparseSequenceVectorTag) const
{
	linbox_check (A.rowdim () == w.size ());

	typename Vector2::const_iterator j = v.begin ();

	_VD.subin (w, w);

	for (; j != v.end (); ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		_VD.axpyin (w, j->second, *i);
	}

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						 VectorCategories::DenseVectorTag,
						 VectorCategories::SparseAssociativeVectorTag) const
{
	linbox_check (A.rowdim () == w.size ());

	typename Vector2::const_iterator j = v.begin ();

	_VD.subin (w, w);

	for (; j != v.end (); ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		_VD.axpyin (w, j->second, *i);
	}

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						 VectorCategories::DenseVectorTag,
						 VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.rowdim () == w.size ());

	typename Vector2::first_type::const_iterator j_idx = v.first.begin ();
	typename Vector2::second_type::const_iterator j_elt = v.second.begin ();

	_VD.subin (w, w);

	for (; j_idx != v.first.end (); ++j_idx, ++j_elt) {
		typename Matrix::ConstColIterator i = A.colBegin () + *j_idx;
		_VD.axpyin (w, *j_elt, *i);
	}

	return w;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector1::iterator j = y.begin ();

	typename Field::Element t;

	for (; j != y.end (); ++j, ++i) {
		_VD.dot (t, x, *i);
		_F.addin (*j, t);
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::SparseSequenceVectorTag) const
{
	typename LinBox::Vector<Field>::Dense t1 (A.coldim ()), t2 (A.rowdim ());

	_VD.copy (t1, x);
	_VD.copy (t2, y);
	axpyin (t2, A, t1);
	_VD.copy (y, t2);

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::SparseAssociativeVectorTag) const
{
	typename LinBox::Vector<Field>::Dense t1 (A.coldim ()), t2 (A.rowdim ());

	_VD.copy (t1, x);
	_VD.copy (t2, y);
	axpyin (t2, A, t1);
	_VD.copy (y, t2);

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::SparseParallelVectorTag) const
{
	typename LinBox::Vector<Field>::Dense t1 (A.coldim ()), t2 (A.rowdim ());

	_VD.copy (t1, x);
	_VD.copy (t2, y);
	axpyin (t2, A, t1);
	_VD.copy (y, t2);

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j = x.begin ();

	for (; j != x.end (); ++j, ++i)
		_VD.axpyin (y, *j, *i);

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::SparseSequenceVectorTag) const
{
	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::iterator j = x.begin ();

	int diff;

	typename LinBox::Vector<Field>::Dense t (A.rowdim ());

	_VD.copy (t, y);

	while (j != x.end ()) {
		_VD.axpyin (t, j->second, *i);
		diff = j->first; ++j;
		diff -= j->first;
		i -= diff;
	}

	_VD.copy (y, t);

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::SparseAssociativeVectorTag) const
{
	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::iterator j = x.begin ();

	int diff;

	typename LinBox::Vector<Field>::Dense t (A.rowdim ());

	_VD.copy (t, y);

	while (j != x.end ()) {
		_VD.axpyin (t, j->second, *i);
		diff = j->first; ++j;
		diff -= j->first;
		i -= diff;
	}

	_VD.copy (y, t);

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<Field>::axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
						    VectorCategories::SparseParallelVectorTag) const
{
	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::iterator j_idx = x.first.begin ();
	typename Vector2::iterator j_elt = x.second.begin ();

	int diff;

	typename LinBox::Vector<Field>::Dense t (A.rowdim ());

	_VD.copy (t, y);

	for (; j_idx != x.first.end (); ++j_elt) {
		_VD.axpyin (t, *j_elt, *i);
		diff = *j_idx; ++j_idx;
		diff -= *j_idx;
		i -= diff;
	}

	_VD.copy (y, t);

	return y;
}

template <class Field>
template <class Matrix1, class Blackbox, class Matrix2>
Matrix1 &MatrixDomain<Field>::blackboxMulLeft (Matrix1 &C, const Blackbox &A, const Matrix2 &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ColIterator i = C.colBegin ();
	typename Matrix2::ConstColIterator j = B.colBegin ();

	for (; i != C.colEnd (); ++i, ++j)
		A.apply (*i, *j);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Blackbox>
Matrix1 &MatrixDomain<Field>::blackboxMulRight (Matrix1 &C, const Matrix2 &A, const Blackbox &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::RowIterator i = C.rowBegin ();
	typename Matrix2::ConstRowIterator j = A.rowBegin ();

	for (; i != C.rowEnd (); ++i, ++j)
		B.applyTranspose (*i, *j);

	return C;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomain<Field>::permuteRowsByRow (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end) const
{
	Iterator i;
	typename Matrix::RowIterator j, k;

	for (i = P_start; i != P_end; ++i) {
		j = A.rowBegin () + i->first;
		k = A.rowBegin () + i->second;

		_VD.swap (*j, *k);
	}

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomain<Field>::permuteRowsByCol (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end) const
{
	typename Matrix::ColIterator j;

	for (j = A.colBegin (); j != A.colEnd (); ++j)
		_VD.permute (*j, P_start, P_end);

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomain<Field>::permuteColsByRow (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end) const
{
	typename Matrix::RowIterator j;

	for (j = A.rowBegin (); j != A.rowEnd (); ++j)
		_VD.permute (*j, P_start, P_end);

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomain<Field>::permuteColsByCol (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end) const
{
	Iterator i;
	typename Matrix::ColIterator j, k;

	for (i = P_start; i != P_end; ++i) {
		j = A.colBegin () + i->first;
		k = A.colBegin () + i->second;

		_VD.swap (*j, *k);
	}

	return A;
}

/* FIXME: These methods are undocumented, and I'm unclear what they are supposed
 * to do 
 */

#if 0

/*M1<-M2**k;
 */
template<class Matrix1, class Matrix2>
Matrix1& MatrixDomain::pow_apply(Matrix1& M1, const Matrix2& M2, unsigned long int k) const
{
	linbox_check((M1.rowdim()==M1.coldim())&&
		     (M2.rowdim()==M2.coldim())&&
		     (M1.rowdim()==M2.rowdim()));

  
	typename Matrix1::RawIterator p=M1.rawBegin();
	for(;p!=M1.rawEnd();++p)
		M1.field().init(*p,0);
	for(p=M1.rawBegin();p<M1.rawEnd();)
	{
		M1.field().init(*p,1);
		p=p+M1.rowdim()+1;
	}
    
    
	for(int i=0;i<k;++i)
		mulin_R(M1,M2);
    
	return M1;
}
    
  
template<class Matrix1, class Matrix2>
Matrix1& MatrixDomain::pow_horn(Matrix1& M1, const Matrix2& M2, unsigned long int k) const
{
	linbox_check((M1.rowdim()==M1.coldim())&&
		     (M2.rowdim()==M2.coldim())&&
		     (M1.rowdim()==M2.rowdim()));
    
	if(k==0)
	{
		typename Matrix1::RawIterator p=M1.rawBegin();
		for(;p!=M1.rawEnd();++p)
			M1.field().init(*p,0);
		for(p=M1.rawBegin();p<M1.rawEnd();)
		{
			M1.field().init(*p,1);
			p+=M1.rowdim()+1;
		}
		return M1;
	}
    
	typename Matrix1::RawIterator p1;
	typename Matrix2::ConstRawIterator p2;
	for(p1=M1.rawBegin(),p2=M2.rawBegin();p1!=M1.rawEnd();++p1,++p2)
		M1.field().assign(*p1,*p2);
  
	std::vector<bool> bit;
	bit.reserve(sizeof(unsigned long)*4);
	while(k>0)
	{
		bit.push_back(k%2);
		k/=2;
	};

    
	std::vector<bool>::reverse_iterator p=bit.rbegin();
	++p;
	Matrix1 temp(M1);
	for(;p!=bit.rend();++p)
	{
		temp=M1;
		mulin_L(M1,temp);
		if(*p)
			mulin_L(M1,M2);

	}
      
	return M1;     
}

#endif

} // namespace LinBox

#endif // __LINBOX_matrix_domain_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
