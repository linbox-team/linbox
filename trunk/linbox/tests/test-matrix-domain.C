
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
 * SparseMatrixBase, and TransposeMatrix
 */

/* ERRORS:
 *
 * X- Ambiguous specializations with RolColMatrixTag
 *  - Can't use leftMulin or rightMulin on some tests
 *  - VectorDomain needs subin, addin, etc. with multiple vector representations
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/sparse.h"
#include "linbox/blackbox/matrix-blackbox.h"
#include "linbox/matrix/dense-submatrix.h"

#include "test-common.h"

using namespace std;
using namespace LinBox;

class SingularMatrix {};

template <class Matrix>
TransposeMatrix<Matrix> transpose (Matrix &M)
{
	return TransposeMatrix<Matrix> (M);
}

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

// Version for square matrices

template <class Field, class Matrix>
static bool testInvMulSquare (Field &F, const char *text, const Matrix &M) 
{
	linbox_check (M.rowdim () == M.coldim ());

	ostringstream str;

	str << "Testing " << text << " matrix multiplication (square)" << ends;
	commentator.start (str.str ().c_str (), "testInvMulSquare");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.rowdim (), M.rowdim ());
	DenseMatrixBase<typename Field::Element> M2 (M.rowdim (), M.rowdim ());

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	try {
		inv (MD, Minv, M);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvMul");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.mul (M2, M, Minv);

	report << "Computed product M Minv:" << endl;
	MD.write (report, M2);

	if (!MD.areEqual (M2, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M Minv is not the identity" << endl;
		ret = false;
	}

	MD.mul (M2, Minv, M);

	report << "Computed product Minv M:" << endl;
	MD.write (report, M2);

	if (!MD.areEqual (M2, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvMulSquare");

	return ret;
}

// Version for over-determined matrices

template <class Field, class Matrix>
static bool testInvMulOver (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.coldim () <= M.rowdim ());

	ostringstream str;

	str << "Testing " << text << " matrix multiplication (over-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvMulOver");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.coldim (), M.coldim ());
	DenseMatrixBase<typename Field::Element> M2 (M.coldim (), M.coldim ());
	DenseMatrixBase<typename Field::Element> M3 (M.coldim (), M.coldim ());

	DenseMatrixBase<typename Field::Element> MTM (M.coldim (), M.coldim ());

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row>
		stream (F, M.coldim ());

	DenseMatrixBase<typename Field::Element> I (M.coldim (), M.coldim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	MD.mul (MTM, transpose (M), M);

	try {
		inv (MD, Minv, MTM);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvMul");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.mul (M2, transpose(M), M);
	MD.mul (M3, M2, Minv);

	report << "Computed product M^T M Minv:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M^T M Minv is not the identity" << endl;
		ret = false;
	}

	M2.resize (M.coldim (), M.rowdim ());

	MD.mul (M2, Minv, transpose (M));
	MD.mul (M3, M2, M);

	report << "Computed product Minv M^T M:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M^T M is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvMulOver");

	return ret;
}

// Version for under-determined matrices

template <class Field, class Matrix>
static bool testInvMulUnder (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.rowdim () <= M.coldim ());

	ostringstream str;

	str << "Testing " << text << " matrix multiplication (under-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvMulUnder");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.rowdim (), M.rowdim ());
	DenseMatrixBase<typename Field::Element> M2 (M.rowdim (), M.rowdim ());
	DenseMatrixBase<typename Field::Element> M3 (M.rowdim (), M.rowdim ());

	DenseMatrixBase<typename Field::Element> MMT (M.rowdim (), M.rowdim ());

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	MD.mul (MMT, M, transpose (M));

	try {
		inv (MD, Minv, MMT);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvMul");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.mul (M2, M, transpose (M));
	MD.mul (M3, M2, Minv);

	report << "Computed product M M^T Minv:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M M^T Minv is not the identity" << endl;
		ret = false;
	}

	M2.resize (M.rowdim (), M.coldim ());

	MD.mul (M2, Minv, M);
	MD.mul (M3, M2, transpose (M));

	report << "Computed product Minv M M^T:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M M^T is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvMulUnder");

	return ret;
}

/* Test 6: inv and leftMulin
 *
 * Return true on success and false on failure
 */

// Version for square matrices

template <class Field, class Matrix>
static bool testInvLeftMulinSquare (Field &F, const char *text, const Matrix &M) 
{
	linbox_check (M.rowdim () == M.coldim ());

	ostringstream str;

	str << "Testing " << text << " left in-place matrix multiplication (square)" << ends;
	commentator.start (str.str ().c_str (), "testInvLeftMulinSquare");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.rowdim (), M.rowdim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	try {
		inv (MD, Minv, M);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvLeftMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return false;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.leftMulin (M, Minv);

	report << "Computed product M Minv:" << endl;
	MD.write (report, Minv);

	if (!MD.areEqual (Minv, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvLeftMulinSquare");

	return ret;
}

// Version for over-determined matrices

template <class Field, class Matrix>
static bool testInvLeftMulinOver (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.coldim () <= M.rowdim ());

	ostringstream str;

	str << "Testing " << text << " left in-place matrix multiplication (over-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvLeftMulinOver");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.coldim (), M.coldim ());
	DenseMatrixBase<typename Field::Element> MTM (M.coldim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.coldim ());

	DenseMatrixBase<typename Field::Element> I (M.coldim (), M.coldim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	MD.mul (MTM, transpose (M), M);

	try {
		inv (MD, Minv, MTM);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvLeftMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return false;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.leftMulin (Minv, MTM);

	report << "Computed product Minv M^T M:" << endl;
	MD.write (report, MTM);

	if (!MD.areEqual (MTM, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M^T M is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvLeftMulinOver");

	return ret;
}

// Version for under-determined matrices

template <class Field, class Matrix>
static bool testInvLeftMulinUnder (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.rowdim () <= M.coldim ());

	ostringstream str;

	str << "Testing " << text << " left in-place matrix multiplication (under-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvLeftMulinUnder");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.rowdim (), M.rowdim ());
	DenseMatrixBase<typename Field::Element> MMT (M.rowdim (), M.rowdim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	MD.mul (MMT, M, transpose (M));

	try {
		inv (MD, Minv, MMT);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvLeftMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return false;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.leftMulin (Minv, MMT);

	report << "Computed product Minv M M^T:" << endl;
	MD.write (report, MMT);

	if (!MD.areEqual (MMT, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M M^T is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvLeftMulinUnder");

	return ret;
}

/* Test 7: inv and rightMulin
 *
 * Return true on success and false on failure
 */

// Version for square matrices

template <class Field, class Matrix>
static bool testInvRightMulinSquare (Field &F, const char *text, const Matrix &M) 
{
	linbox_check (M.rowdim () == M.coldim ());

	ostringstream str;

	str << "Testing " << text << " right in-place matrix multiplication (square)" << ends;
	commentator.start (str.str ().c_str (), "testInvRightMulinSquare");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.rowdim (), M.rowdim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	try {
		inv (MD, Minv, M);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvRightMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.rightMulin (Minv, M);

	report << "Computed product Minv M:" << endl;
	MD.write (report, Minv);

	if (!MD.areEqual (Minv, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvRightMulinSquare");

	return ret;
}

// Version for over-determined matrices

template <class Field, class Matrix>
static bool testInvRightMulinOver (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.coldim () <= M.rowdim ());

	ostringstream str;

	str << "Testing " << text << " right in-place matrix multiplication (over-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvRightMulinOver");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.coldim (), M.coldim ());
	DenseMatrixBase<typename Field::Element> MTM (M.coldim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.coldim ());

	DenseMatrixBase<typename Field::Element> I (M.coldim (), M.coldim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	MD.mul (MTM, transpose (M), M);

	try {
		inv (MD, Minv, MTM);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvRightMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.rightMulin (MTM, Minv);

	report << "Computed product M^T M Minv:" << endl;
	MD.write (report, MTM);

	if (!MD.areEqual (MTM, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M^T M Minv is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvRightMulinOver");

	return ret;
}

// Version for under-determined matrices

template <class Field, class Matrix>
static bool testInvRightMulinUnder (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.rowdim () <= M.coldim ());

	ostringstream str;

	str << "Testing " << text << " right in-place matrix multiplication (under-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvRightMulinUnder");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> Minv (M.rowdim (), M.rowdim ());
	DenseMatrixBase<typename Field::Element> MMT (M.rowdim (), M.rowdim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename DenseMatrixBase<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrixBase<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrixBase<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	MD.mul (MMT, M, transpose (M));

	try {
		inv (MD, Minv, MMT);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvRightMulin");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.rightMulin (MMT, Minv);

	report << "Computed product M M^T Minv:" << endl;
	MD.write (report, MMT);

	if (!MD.areEqual (MMT, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M M^T Minv is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvRightMulinUnder");

	return ret;
}

/* Test 8: add-mul and axpyin
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testAddMulAxpyin (Field &F, const char *text, Matrix &M1, Matrix &M2, Matrix &M3) 
{
	ostringstream str;

	str << "Testing " << text << " matrix add-mul, axpyin" << ends;
	commentator.start (str.str ().c_str (), "testAddMulAxpyin");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrixBase<typename Field::Element> M4 (M2.rowdim (), M3.coldim ());
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
	typename LinBox::Vector<Field>::Dense v (M.coldim ());
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

template <class Field, class Vector, class Blackbox>
static bool testLeftBlackboxMul (Field &F, const char *text, const Blackbox &A,
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

	MD.blackboxMulLeft (AI, A, I);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Output matrix AI:" << endl;
	MD.write (report, AI);

	typename LinBox::Vector<Field>::Dense v (A.coldim ()), w1 (A.rowdim ()), w2 (A.rowdim ());

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

template <class Field, class Vector, class Blackbox>
static bool testRightBlackboxMul (Field &F, const char *text, const Blackbox &A,
				  VectorStream<Vector> &stream) 
{
	ostringstream str;

	str << "Testing " << text << " matrix-black box right mul" << ends;
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

	MD.blackboxMulRight (IA, I, A);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Output matrix IA:" << endl;
	MD.write (report, IA);

	typename LinBox::Vector<Field>::Dense v (A.coldim ()), w1 (A.rowdim ()), w2 (A.rowdim ());

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

std::ostream &reportPermutation
	(std::ostream &out,
	 const std::vector<std::pair<unsigned int, unsigned int> > &P) 
{
	std::vector<std::pair<unsigned int, unsigned int> >::const_iterator i;

	for (i = P.begin (); i != P.end (); ++i)
		out << "(" << i->first << " " << i->second << ")";

	return out;
}

template <class Field, class Matrix>
bool testPermutation (const Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " permutations" << ends;
	commentator.start (str.str ().c_str (), "testPermutation");

	bool ret = true;

	MatrixDomain<Field> MD (F);
	MersenneTwister MT (time (NULL));

	typename MatrixDomain<Field>::Permutation P, Pinv;

	// Create a random row permutation
	for (unsigned int i = 0; i < M.rowdim (); ++i) {
		unsigned int row1, row2;

		do {
			row1 = MT.randomInt () % M.rowdim ();
			row2 = MT.randomInt () % M.rowdim ();
		} while (row1 == row2);

		P.push_back (typename MatrixDomain<Field>::Transposition (row1, row2));
	}

	// Construct the inverse of this transposition
	Pinv.resize (P.size ());
	std::copy (P.begin (), P.end (), Pinv.begin ());
	std::reverse (Pinv.begin (), Pinv.end ());

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Permutation P:    ";
	reportPermutation (report, P) << endl;
	report << "Permutation P^-1: ";
	reportPermutation (report, Pinv) << endl;

	// Apply the permutation and then its inverse to a copy of M
	Matrix M1 (M);

	MD.permuteRows (M1, P.begin (), P.end ());
 	MD.permuteRows (M1, Pinv.begin (), Pinv.end ());

	report << "Output matrix P^-1 PM:" << endl;
	MD.write (report, M1);

	// Compare M and M1
	if (!MD.areEqual (M, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: M != P^-1 PM" << endl;
		ret = false;
	}

	// Now we do exactly the same thing with columns
	P.clear ();
	Pinv.clear ();

	// Create a random column permutation
	for (unsigned int i = 0; i < M.coldim (); ++i) {
		unsigned int col1, col2;

		do {
			col1 = MT.randomInt () % M.coldim ();
			col2 = MT.randomInt () % M.coldim ();
		} while (col1 == col2);

		P.push_back (typename MatrixDomain<Field>::Transposition (col1, col2));
	}

	// Construct the inverse of this transposition
	Pinv.resize (P.size ());
	std::copy (P.begin (), P.end (), Pinv.begin ());
	std::reverse (Pinv.begin (), Pinv.end ());

	report << "Permutation P:    ";
	reportPermutation (report, P) << endl;
	report << "Permutation P^-1: ";
	reportPermutation (report, Pinv) << endl;

	// Apply the permutation and then its inverse to a copy of M
	MD.permuteColumns (M1, P.begin (), P.end ());
	MD.permuteColumns (M1, Pinv.begin (), Pinv.end ());

	report << "Output matrix MPP^-1:" << endl;
	MD.write (report, M1);

	// Compare M and M1
	if (!MD.areEqual (M, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: M != MPP^-1" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testPermutation");

	return ret;
}

template <class Field, class Blackbox, class Matrix>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       const Blackbox &A,
		       unsigned int iterations,
		       MatrixCategories::RowColMatrixTag)
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, A.coldim (), iterations);
	
	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testSubinIsZero (F, text, M1)) pass = false;
	if (!testAddNegSub (F, text, M1, M2)) pass = false;
	if (!testAddinNeginSub (F, text, M1, M2)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testInvMulSquare (F, text, M1)) pass = false;
		if (!testInvLeftMulinSquare (F, text, M1)) pass = false;
		if (!testInvRightMulinSquare (F, text, M1)) pass = false;
	}
	else if (M1.coldim () < M1.rowdim ()) {
		if (!testInvMulOver (F, text, M1)) pass = false;
		if (!testInvLeftMulinOver (F, text, M1)) pass = false;
		if (!testInvRightMulinOver (F, text, M1)) pass = false;
	}
	else if (M1.rowdim () < M1.coldim ()) {
		if (!testInvMulUnder (F, text, M1)) pass = false;
		if (!testInvLeftMulinUnder (F, text, M1)) pass = false;
		if (!testInvRightMulinUnder (F, text, M1)) pass = false;
	}

	if (!testAddMulAxpyin (F, text, M1, M2, M3)) pass = false;
	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;
	if (!testLeftBlackboxMul (F, text, A, stream)) pass = false;
	if (!testRightBlackboxMul (F, text, A, stream)) pass = false;
	if (!testPermutation (F, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Blackbox, class Matrix>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       const Blackbox &A,
		       unsigned int iterations,
		       MatrixCategories::RowMatrixTag) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, A.coldim (), iterations);

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testSubinIsZero (F, text, M1)) pass = false;
	if (!testAddNegSub (F, text, M1, M2)) pass = false;
	if (!testAddinNeginSub (F, text, M1, M2)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testInvMulSquare (F, text, M1)) pass = false;
		if (!testInvLeftMulinSquare (F, text, M1)) pass = false;
	}
	else if (M1.rowdim () < M1.coldim ()) {
		if (!testInvMulUnder (F, text, M1)) pass = false;
		if (!testInvLeftMulinUnder (F, text, M1)) pass = false;
	}

	if (!testAddMulAxpyin (F, text, M1, M2, M3)) pass = false;
	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;
	if (!testLeftBlackboxMul (F, text, A, stream)) pass = false;
	if (!testPermutation (F, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Blackbox, class Matrix>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       const Blackbox &A,
		       unsigned int iterations,
		       MatrixCategories::ColMatrixTag) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, A.coldim (), iterations);

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testSubinIsZero (F, text, M1)) pass = false;
	if (!testAddNegSub (F, text, M1, M2)) pass = false;
	if (!testAddinNeginSub (F, text, M1, M2)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testInvMulSquare (F, text, M1)) pass = false;
		if (!testInvRightMulinSquare (F, text, M1)) pass = false;
	}
	else if (M1.coldim () < M1.rowdim ()) {
		if (!testInvMulOver (F, text, M1)) pass = false;
		if (!testInvRightMulinOver (F, text, M1)) pass = false;
	}

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

	static long n = 50;
	static long m = 50;
	static long k = 10;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set row of test matrices to N.", TYPE_INT,     &n },
		{ 'm', "-m M", "Set column of test vectors to M.", TYPE_INT,     &m },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT,     &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	typedef Modular<uint32> Field;
	typedef Field::Element Element;

	Field F (q);

	commentator.start("Matrix domain test suite", "MatrixDomain");

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (2);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	DenseMatrixBase<Element> M1 (n, m);
	DenseMatrixBase<Element> M2 (n, m);
	DenseMatrixBase<Element> M3 (m, m);
	MatrixBlackbox<Field, DenseMatrixBase<Field::Element> > A1 (F, n, m);

	RandomDenseStream<Field, DenseMatrixBase<Element>::Row> stream1 (F, m);

	DenseMatrixBase<Element>::RowIterator i;

	for (i = M1.rowBegin (); i != M1.rowEnd (); ++i)
		stream1 >> *i;

	for (i = M2.rowBegin (); i != M2.rowEnd (); ++i)
		stream1 >> *i;

	for (i = M3.rowBegin (); i != M3.rowEnd (); ++i)
		stream1 >> *i;

	for (i = A1.rep ().rowBegin (); i != A1.rep ().rowEnd (); ++i)
		stream1 >> *i;

	if (!testMatrixDomain (F, "dense", M1, M2, M3, A1, iterations,
			       MatrixTraits<DenseMatrixBase<Element> >::MatrixCategory ()))
		pass = false;

	SparseMatrixBase<Element> M4 (n, m);
	SparseMatrixBase<Element> M5 (n, m);
	SparseMatrixBase<Element> M6 (m, m);
	MatrixBlackbox<Field, SparseMatrixBase<Field::Element> > A2 (F, n, m);

	RandomSparseStream<Field, SparseMatrixBase<Element>::Row> stream2 (F, (double) k / (double) n, m);

	SparseMatrixBase<Element>::RowIterator i2;

	for (i2 = M4.rowBegin (); i2 != M4.rowEnd (); ++i2)
		stream2 >> *i2;

	for (i2 = M5.rowBegin (); i2 != M5.rowEnd (); ++i2)
		stream2 >> *i2;

	for (i2 = M6.rowBegin (); i2 != M6.rowEnd (); ++i2)
		stream2 >> *i2;

	for (i2 = A2.rep ().rowBegin (); i2 != A2.rep ().rowEnd (); ++i2)
		stream2 >> *i2;

	if (!testMatrixDomain (F, "sparse row-wise", M4, M5, M6, A2, iterations,
			       MatrixTraits<SparseMatrixBase<Element> >::MatrixCategory ()))
		pass = false;

	TransposeMatrix<SparseMatrixBase<Element> > M7 (M4);
	TransposeMatrix<SparseMatrixBase<Element> > M8 (M5);
	TransposeMatrix<SparseMatrixBase<Element> > M9 (M6);

	if (!testMatrixDomain (F, "sparse column-wise", M7, M8, M9, A2, iterations,
			       MatrixTraits<TransposeMatrix<SparseMatrixBase<Element> > >::MatrixCategory ()))
		pass = false;

	commentator.stop("Matrix domain test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
