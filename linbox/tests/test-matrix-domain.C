/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-matrix-domain.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Test suite for MatrixDomain. This also effectively tests DenseMatrixBase,
 * SparseMatrix0Base, and TransposeMatrix
 */

/* ERRORS:
 *
 * X- Ambiguous specializations with RolColMatrixTag
 *  - Can't use leftMulin or rightMulin on some tests
 *  - VectorDomain needs subin, addin, etc. with multiple vector representations
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/vector-domain.h"
#include "linbox/field/matrix-domain.h"
#include "linbox/vector/stream.h"
#include "linbox/blackbox/dense-base.h"
#include "linbox/blackbox/sparse-base.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/dense-submatrix.h"

#include "test-common.h"

using namespace std;
using namespace LinBox;

class SingularMatrix {};

template <class Field, class Matrix1, class Matrix2>
void eliminate (MatrixDomain<Field> &MD, Matrix1 &M, Matrix2 &pivotRow,
		size_t row, size_t col, size_t rowdim, size_t coldim) 
{
	DenseMatrixBase<typename Matrix1::Element> pivotCol (rowdim, 1);
	DenseSubmatrix<typename Matrix1::Element> realPivotCol (M, row, col, rowdim, 1);
	DenseSubmatrix<typename Matrix1::Element> block (M, row, col, rowdim, coldim);

	MD.neg (pivotCol, realPivotCol);
	MD.axpyin (block, pivotCol, pivotRow);
}

/* Dumb elimination code
 *
 * This computes the inverse of a nonsingular matrix. It throws SingularMatrix
 * in the case that no inverse exists
 */

template <class Field, class Matrix1, class Matrix2>
Matrix1 &inv (MatrixDomain<Field> &MD, Matrix1 &res, const Matrix2 &A) 
{
	linbox_check (res.coldim () == A.coldim ());
	linbox_check (res.rowdim () == A.rowdim ());

	DenseMatrixBase<typename Matrix1::Element> M (res.rowdim (), res.coldim () * 2);
	DenseSubmatrix<typename Matrix1::Element> M1 (M, 0, 0, res.rowdim (), res.coldim ());
	DenseSubmatrix<typename Matrix1::Element> M2 (M, 0, res.coldim (), res.rowdim (), res.coldim ());

	StandardBasisStream<Field, typename DenseSubmatrix<typename Matrix1::Element>::Row> stream (MD.field (), res.coldim ());
	typename DenseSubmatrix<typename Matrix1::Element>::RowIterator ip = M2.rowBegin ();

	for (; ip != M2.rowEnd (); ++ip)
		stream >> *ip;

	MD.copy (M1, A);

	unsigned int idx;
	typename Field::Element Mjj_inv;

	for (idx = 0; idx < M.rowdim (); ++idx) {
		if (MD.field ().isZero (M.getEntry (idx, idx))) {
			typename DenseMatrixBase<typename Matrix1::Element>::ColIterator col;
			typename DenseMatrixBase<typename Matrix1::Element>::Col::iterator i;
			unsigned int c_idx = idx + 1;

			col = M.colBegin () + idx;
			i = col->begin () + idx + 1;

			while (MD.field ().isZero (*i) && i != col->end ()) ++i, ++c_idx;

			if (i == col->end ())
				throw SingularMatrix ();
			else {
				typename DenseMatrixBase<typename Matrix1::Element>::RowIterator row1, row2;

				row1 = M.rowBegin () + idx;
				row2 = M.rowBegin () + c_idx;

				std::swap_ranges (row1->begin () + idx, row1->end (), row2->begin () + idx);
			}
		}

		MD.field ().inv (Mjj_inv, M.getEntry (idx, idx));
		DenseSubmatrix<typename Matrix1::Element> realPivotRow (M, idx, idx, 1, M.coldim () - idx);
		MD.mulin (realPivotRow, Mjj_inv);

		if (idx > 0)
			eliminate (MD, M, realPivotRow, 0, idx, idx, M.coldim () - idx);

		if (idx < M.rowdim () - 1)
			eliminate (MD, M, realPivotRow, idx + 1, idx,
				   M.rowdim () - idx - 1, M.coldim () - idx);
	}

	MD.copy (res, M2);
	return res;
}

/* Test 1: Copy and areEqual
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testCopyEqual (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix copy, areEqual" << ends;
	commentator.start (str.str ().c_str (), "testCopyEqual");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	MD.copy (M1, M);

	report << "Output matrix M1:" << endl;
	MD.write (report, M1);

	if (!MD.areEqual (M1, M)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M and M1 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testCopyEqual");

	return ret;
}

/* Test 2: subin and isZero
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testSubinIsZero (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix subin, isZero" << ends;
	commentator.start (str.str ().c_str (), "testSubinIsZero");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	MD.copy (M1, M);
	MD.subin (M1, M);

	report << "Output matrix M1:" << endl;
	MD.write (report, M1);

	if (!MD.isZero (M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M1 is not zero" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSubinIsZero");

	return ret;
}

/* Test 3: add-neg and sub
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testAddNegSub (Field &F, const char *text, const Matrix &M1, const Matrix &M2) 
{
	ostringstream str;

	str << "Testing " << text << " matrix add-neg, sub" << ends;
	commentator.start (str.str ().c_str (), "testAddNegSub");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M3 (M1.rowdim (), M1.coldim ());
	DenseMatrixBase<typename Field::Element> M4 (M1.rowdim (), M1.coldim ());
	DenseMatrixBase<typename Field::Element> M5 (M1.rowdim (), M1.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M1:" << endl;
	MD.write (report, M1);

	report << "Input matrix M2:" << endl;
	MD.write (report, M2);

	MD.neg (M3, M2);
	MD.add (M4, M1, M3);
	MD.sub (M5, M1, M2);

	report << "Output matrix M3 = -M2:" << endl;
	MD.write (report, M3);

	report << "Output matrix M4 = M1 + -M2:" << endl;
	MD.write (report, M4);

	report << "Output matrix M5 = M1 - M2:" << endl;
	MD.write (report, M5);

	if (!MD.areEqual (M4, M5)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M4 and M5 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAddNegSub");

	return ret;
}

/* Test 4: addin-negin and sub
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testAddinNeginSub (Field &F, const char *text, const Matrix &M1, const Matrix &M2) 
{
	ostringstream str;

	str << "Testing " << text << " matrix addin-negin, sub" << ends;
	commentator.start (str.str ().c_str (), "testAddinNeginSub");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M3 (M1.rowdim (), M1.coldim ());
	DenseMatrixBase<typename Field::Element> M4 (M1.rowdim (), M1.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M1:" << endl;
	MD.write (report, M1);

	report << "Input matrix M2:" << endl;
	MD.write (report, M2);

	MD.copy (M3, M2);
	MD.negin (M3);

	report << "Output matrix M3 = -M2:" << endl;
	MD.write (report, M3);

	MD.addin (M3, M1);

	report << "Output matrix M3 = M1 + -M2:" << endl;
	MD.write (report, M3);

	MD.sub (M4, M1, M2);

	report << "Output matrix M4 = M1 - M2:" << endl;
	MD.write (report, M4);

	if (!MD.areEqual (M3, M4)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M3 and M4 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAddinNeginSub");

	return ret;
}

/* Test 5: inv and mul
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testInvMul (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix multiplication" << ends;
	commentator.start (str.str ().c_str (), "testInvMul");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M1 (M.rowdim (), M.coldim ());
	DenseMatrixBase<typename Field::Element> M2 (M.rowdim (), M.coldim ());
	DenseMatrixBase<typename Field::Element> M3 (M.rowdim (), M.coldim ());

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.coldim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	try {
		inv (MD, M1, M);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvMul");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse M1:" << endl;
	MD.write (report, M1);

	MD.mul (M2, M, M1);

	report << "Computed product M2:" << endl;
	MD.write (report, M2);

	MD.subin (M2, I);

	if (!MD.isZero (M2)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M2 is not the identity" << endl;
		ret = false;
	}

	MD.mul (M3, M1, M);

	report << "Computed product M3:" << endl;
	MD.write (report, M3);

	MD.subin (M3, I);

	if (!MD.isZero (M3)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M3 is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvMul");

	return ret;
}

/* Test 6: inv and leftMulin
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testInvLeftMulin (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " left in-place matrix multiplication" << ends;
	commentator.start (str.str ().c_str (), "testInvLeftMulin");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	try {
		inv (MD, M1, M);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvLeftMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return false;
	}

	report << "Computed inverse M1:" << endl;
	MD.write (report, M1);

	MD.leftMulin (M, M1);

	report << "Computed product M1:" << endl;
	MD.write (report, M1);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.coldim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	MD.subin (M1, I);

	if (!MD.isZero (M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M1 is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvLeftMulin");

	return ret;
}

/* Test 7: inv and rightMulin
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testInvRightMulin (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " right in-place matrix multiplication" << ends;
	commentator.start (str.str ().c_str (), "testInvRightMulin");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	try {
		inv (MD, M1, M);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvRightMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse M1:" << endl;
	MD.write (report, M1);

	MD.rightMulin (M1, M);

	report << "Computed product M1:" << endl;
	MD.write (report, M1);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.coldim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	MD.subin (M1, I);

	if (!MD.isZero (M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M1 is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvRightMulin");

	return ret;
}

/* Test 8: add-mul and axpyin
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testAddMulAxpyin (Field &F, const char *text, const Matrix &M1, const Matrix &M2, const Matrix &M3) 
{
	ostringstream str;

	str << "Testing " << text << " matrix add-mul, axpyin" << ends;
	commentator.start (str.str ().c_str (), "testAddMulAxpyin");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M4 (M1.rowdim (), M1.coldim ());
	DenseMatrixBase<typename Field::Element> M5 (M1.rowdim (), M1.coldim ());
	DenseMatrixBase<typename Field::Element> M6 (M1.rowdim (), M1.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M1:" << endl;
	MD.write (report, M1);

	report << "Input matrix M2:" << endl;
	MD.write (report, M2);

	report << "Input matrix M3:" << endl;
	MD.write (report, M3);

	MD.mul (M4, M2, M3);
	MD.add (M5, M1, M4);

	report << "Computed matrix M5 = M1 + M2 * M3 (add-mul):" << endl;
	MD.write (report, M5);

	MD.copy (M6, M1);
	MD.axpyin (M6, M2, M3);

	report << "Computed matrix M6 = M1 + M2 * M3 (axpyin):" << endl;
	MD.write (report, M6);

	if (!MD.areEqual (M5, M6)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M5 and M6 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAddMulAxpyin");

	return ret;
}

/* Test 9: m-v mul by e_i and sub
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testMVMulSub (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix-vector mul" << ends;
	commentator.start (str.str ().c_str (), "testMVMulSub");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, M.rowdim ());
	typename LinBox::Vector<Field>::Dense v (M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::ColIterator i = M1.colBegin ();

	for (; i != M1.colEnd (); ++i) {
		stream >> v;
		MD.vectorMul (*i, M, v);
	}

	report << "Output matrix M1:" << endl;
	MD.write (report, M1);

	MD.subin (M1, M);

	if (!MD.isZero (M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M and M1 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testMVMulSub");

	return ret;
}

/* Test 10: m-v axpy by e_i
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testMVAxpy (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix-vector axpy" << ends;
	commentator.start (str.str ().c_str (), "testMVAxpy");

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M1:" << endl;
	MD.write (report, M1);

	StandardBasisStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, M.rowdim ());
	typename LinBox::Vector<Field>::Dense v (M.coldim ()), w (M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = M1.rowBegin ();

	VD.subin (w, w);

	for (; i != M1.rowEnd (); ++i) {
		stream >> v;
		MD.vectorAxpyin (w, M, v);
	}

	report << "Output vector w:" << endl;
	VD.write (report, w);

	typename Field::Element one;
	F.init (one, 1);
	typename LinBox::Vector<Field>::Dense z (M.coldim (), one), w1 (M.rowdim ());

	MD.vectorMul (w1, M, z);

	report << "Output vector w1:" << endl;
	VD.write (report, w1);

	if (!VD.areEqual (w, w1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: VectorDomain reported vectors w and w1 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testMVAxpy");

	return ret;
}

/* Test 11: black box multiply by I on the left, test on random vectors
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testLeftBlackboxMul (Field &F, const char *text, const BlackboxArchetype<Vector> &A,
				 VectorStream<Vector> &stream) 
{
	ostringstream str;

	str << "Testing " << text << " matrix-black box left mul" << ends;
	commentator.start (str.str ().c_str (), "testLeftBlackboxMul");

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> I (A.coldim (), A.coldim ());
	DenseMatrixBase<typename Field::Element> AI (A.rowdim (), A.coldim ());

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> Istream (F, A.coldim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i;

	for (i = I.rowBegin (); i != I.rowEnd (); ++i)
		Istream >> *i;

	MD.blackboxMul (AI, A, I);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Output matrix AI:" << endl;
	MD.write (report, AI);

	typename LinBox::Vector<Field>::Dense v (A.rowdim ()), w1 (A.coldim ()), w2 (A.coldim ());

	while (stream) {
		stream >> v;

		A.apply (w1, v);
		MD.vectorMul (w2, AI, v);

		if (!VD.areEqual (w1, w2)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Av != AIv" << endl;
			ret = false;
		}
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testLeftBlackboxMul");

	return ret;
}

/* Test 12: black box multiply by I on the right, test on random vectors
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testRightBlackboxMul (Field &F, const char *text, const BlackboxArchetype<Vector> &A,
				  VectorStream<Vector> &stream) 
{
	ostringstream str;

	str << "Testing " << text << " matrix-black box mul" << ends;
	commentator.start (str.str ().c_str (), "testRightBlackboxMul");

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> I (A.rowdim (), A.rowdim ());
	DenseMatrixBase<typename Field::Element> IA (A.rowdim (), A.coldim ());

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> Istream (F, A.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i;

	for (i = I.rowBegin (); i != I.rowEnd (); ++i)
		Istream >> *i;

	MD.blackboxMul (IA, I, A);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Output matrix IA:" << endl;
	MD.write (report, IA);

	typename LinBox::Vector<Field>::Dense v (A.rowdim ()), w1 (A.coldim ()), w2 (A.coldim ());

	while (stream) {
		stream >> v;

		A.apply (w1, v);
		MD.vectorMul (w2, IA, v);

		if (!VD.areEqual (w1, w2)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Av != IAv" << endl;
			ret = false;
		}
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRightBlackboxMul");

	return ret;
}

template <class Field, class Vector, class Matrix, class Trait>
bool testMatrixDomain (const Field &F, const char *text,
		       const Matrix &M1, const Matrix &M2, const Matrix &M3,
		       const BlackboxArchetype<Vector> &A,
		       unsigned int iterations,
		       MatrixCategories::RowColMatrixTag<Trait>) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, stream.dim (), iterations);
	
	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testSubinIsZero (F, text, M1)) pass = false;
	if (!testAddNegSub (F, text, M1, M2)) pass = false;
	if (!testAddinNeginSub (F, text, M1, M2)) pass = false;
	if (!testInvMul (F, text, M1)) pass = false;
	if (!testInvLeftMulin (F, text, M1)) pass = false;
	if (!testInvRightMulin (F, text, M1)) pass = false;
	if (!testAddMulAxpyin (F, text, M1, M2, M3)) pass = false;
	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;
	if (!testLeftBlackboxMul (F, text, A, stream)) pass = false;
	if (!testRightBlackboxMul (F, text, A, stream)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Vector, class Matrix, class Trait>
bool testMatrixDomain (const Field &F, const char *text,
		       const Matrix &M1, const Matrix &M2, const Matrix &M3,
		       const BlackboxArchetype<Vector> &A,
		       unsigned int iterations,
		       MatrixCategories::RowMatrixTag<Trait>) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, stream.dim (), iterations);

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testSubinIsZero (F, text, M1)) pass = false;
	if (!testAddNegSub (F, text, M1, M2)) pass = false;
	if (!testAddinNeginSub (F, text, M1, M2)) pass = false;
	if (!testInvMul (F, text, M1)) pass = false;
	if (!testInvLeftMulin (F, text, M1)) pass = false;
	if (!testAddMulAxpyin (F, text, M1, M2, M3)) pass = false;
	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;
	if (!testLeftBlackboxMul (F, text, A, stream)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Vector, class Matrix, class Trait>
bool testMatrixDomain (const Field &F, const char *text,
		       const Matrix &M1, const Matrix &M2, const Matrix &M3,
		       const BlackboxArchetype<Vector> &A,
		       unsigned int iterations,
		       MatrixCategories::ColMatrixTag<Trait>) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, stream.dim (), iterations);

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testSubinIsZero (F, text, M1)) pass = false;
	if (!testAddNegSub (F, text, M1, M2)) pass = false;
	if (!testAddinNeginSub (F, text, M1, M2)) pass = false;
	if (!testInvMul (F, text, M1)) pass = false;
	if (!testInvRightMulin (F, text, M1)) pass = false;
	if (!testAddMulAxpyin (F, text, M1, M2, M3)) pass = false;
	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;
	if (!testRightBlackboxMul (F, text, A, stream)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long n = 100;
	static long m = 100;
	static long k = 20;
	static integer q = 2147483647U;
	static int iterations = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set row of test matrices to N (default 100)",                           TYPE_INT,     &n },
		{ 'm', "-m M", "Set column of test vectors to M (default 100)",                         TYPE_INT,     &m },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices (default 20)",     TYPE_INT,     &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",                      TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);

	typedef Modular<uint32> Field;
	typedef Field::Element Element;

	Field F (q);

	cout << "Matrix domain test suite" << endl << endl;
	cout.flush ();

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (2);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	DenseMatrixBase<Element> M1 (n, m);
	DenseMatrixBase<Element> M2 (n, m);
	DenseMatrixBase<Element> M3 (n, m);
	DenseMatrix<Field> A1 (F, n, m);

	RandomDenseStream<Field, DenseMatrixBase<Element>::Row> stream1 (F, n);

	DenseMatrixBase<Element>::RowIterator i;

	for (i = M1.rowBegin (); i != M1.rowEnd (); ++i)
		stream1 >> *i;

	for (i = M2.rowBegin (); i != M2.rowEnd (); ++i)
		stream1 >> *i;

	for (i = M3.rowBegin (); i != M3.rowEnd (); ++i)
		stream1 >> *i;

	for (i = A1.rowBegin (); i != A1.rowEnd (); ++i)
		stream1 >> *i;

	if (!testMatrixDomain (F, "dense", M1, M2, M3, A1, iterations,
			       MatrixTraits<DenseMatrixBase<Element> >::MatrixCategory ()))
		pass = false;

	SparseMatrix0Base<Element> M4 (n, m);
	SparseMatrix0Base<Element> M5 (n, m);
	SparseMatrix0Base<Element> M6 (n, m);
	SparseMatrix0<Field> A2 (F, n, m);

	RandomSparseStream<Field, SparseMatrix0Base<Element>::Row> stream2 (F, n, (double) k / (double) n);

	SparseMatrix0Base<Element>::RowIterator i2;

	for (i2 = M4.rowBegin (); i2 != M4.rowEnd (); ++i2)
		stream2 >> *i2;

	for (i2 = M5.rowBegin (); i2 != M5.rowEnd (); ++i2)
		stream2 >> *i2;

	for (i2 = M6.rowBegin (); i2 != M6.rowEnd (); ++i2)
		stream2 >> *i2;

	for (i2 = A2.rowBegin (); i2 != A2.rowEnd (); ++i2)
		stream2 >> *i2;

	if (!testMatrixDomain (F, "sparse row-wise", M4, M5, M6, A2, iterations,
			       MatrixTraits<SparseMatrix0Base<Element> >::MatrixCategory ()))
		pass = false;

	TransposeMatrix<SparseMatrix0Base<Element> > M7 (M4);
	TransposeMatrix<SparseMatrix0Base<Element> > M8 (M5);
	TransposeMatrix<SparseMatrix0Base<Element> > M9 (M6);

	if (!testMatrixDomain (F, "sparse column-wise", M7, M8, M9, A2, iterations,
			       MatrixTraits<TransposeMatrix<SparseMatrix0Base<Element> > >::MatrixCategory ()))
		pass = false;

	return pass ? 0 : -1;
}
