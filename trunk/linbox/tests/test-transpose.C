
/* tests/test-transpose.C
 *
 * ------------------------------------
 *
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
 *.
 */

/*! @file  tests/test-transpose.C
 * @ingroup tests
 *
 * @brief no doc
 *
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>


#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/givaro.h"
#ifdef __LINBOX_HAVE_NTL
#include "linbox/field/ntl.h"
#endif
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/matrix/blas-matrix.h"
// #include "linbox/blackbox/triplesbb.h"
#include "linbox/matrix/sparse-matrix.h"

#include "test-common.h"
#include "test-blackbox.h"

using namespace LinBox;



/*
template <class Field2, class Blackbox>
static bool testBBrebind (const Field2 &F2, const Blackbox& B)
{
    typedef typename Blackbox::template rebind<Field2>::other FBlackbox;

    FBlackbox A(B, F2);

    return testBlackbox(A);
}
*/

/* Test Transpose:
 * Check that Transpose<Blackbox> meets blackbox interface and behaves
 * as transpose.
 *
 * Return true on success and false on failure
 */
template <class Blackbox>
static bool testTransposeBlackbox(Blackbox & A)
{
	typedef typename Blackbox::Field::Element Element;
	commentator().start ("Testing Transpose", "testTranspose", 1);

	Transpose<Blackbox> B(A);

	bool ret = true, ret1;

	size_t m = A.rowdim(), n = A.coldim();
	VectorDomain<typename Blackbox::Field> VD (A.field());
	std::vector<Element> x(n), y(m), z(n), w(m);

	VD.random(x);
	A.apply(y, x);
	B.applyTranspose(w, x);
	ret1 = VD.areEqual(y, w);
	if (not ret1) commentator().report() << "A and B^T disagree, FAIL" << std::endl;
	ret = ret and ret1;

	VD.random(y);
	A.applyTranspose(x, y);
	B.apply(z, y);
	ret1 = VD.areEqual(x, z);
	if (not ret1) commentator().report() << "A^T and B disagree, FAIL" << std::endl;
	ret = ret and ret1;

    ret1 = testBlackbox(B);
	if (not ret1) commentator().report() << "testBlackbox A^T FAIL" << std::endl;
	ret = ret and ret1;

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testTranspose");

	return ret;
}


/* Test getEntry and setEntry of a Transpose<BB>.
 *
 * this test template can be instantiated only if the underlying class BB
 * has getEntry and setEntry (as dense and sparse matrix types should).
 *
 * For such matrices, A and A^T should share memory, so that results of
 * setEntry on one are reflected by getEntry on the other.
 */
template <class Matrix>
bool testTransposeMatrix(Matrix& A) {
	bool ret = true, ret1;
	typename Matrix::Field::Element x, y;
	A.field().init(x);
	A.field().init(y);
	Transpose<Matrix> B(A);
	size_t m = A.rowdim(), n = A.coldim();

	size_t i = (size_t)rand()%m, j = (size_t)rand()%n;
	A.getEntry(x, i, j);
	B.getEntry(y, j, i);
	ret = ret and (ret1 = A.field().areEqual(x, y));
	if (not ret1) commentator().report() << "A, A^T same getentry FAIL" << std::endl;

	i = (size_t)rand()%m, j = (size_t)rand()%n;
	A.setEntry(i, j, x);
	B.getEntry(y, j, i);
	ret = ret and (ret1 = A.field().areEqual(x, y));
	if (not ret1) commentator().report() << "A set, A^T getEntry FAIL" << std::endl;

	i = (size_t)rand()%m, j = (size_t)rand()%n;
	B.setEntry(i, j, x);
	A.getEntry(y, j, i);
	ret = ret and (ret1 = A.field().areEqual(x, y));
	if (not ret1) commentator().report() << "A^T set, A getEntry FAIL" << std::endl;

	return ret;
}
int main (int argc, char **argv)
{
	bool pass = true;

	static size_t m = 8;
	static size_t n = 10;
	static integer q = 101;

	static Argument args[] = {
		{ 'm', "-m M", "Set dimension of test matrices to NxN.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

#ifdef __LINBOX_HAVE_NTL_blah
//        typedef UnparametricField<NTL::zz_p> Field;
        typedef NTL_zz_p Field;
// 	NTL::zz_p::init(q1); // Done in the constructor
#else
	typedef Modular<int32_t> Field ;
#endif
    //typedef GivaroZpz< Givaro::Std32> Field2;
	Field F(q);

	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);

	commentator().start("transpose black box test suite", "transpose");

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	commentator().start("test on ScalarMatrix");
	Field::Element s; F.init(s, 5);
	ScalarMatrix<Field> A(F, n, n, s);
	pass = pass and testTransposeBlackbox(A);
	commentator().stop("test on ScalarMatrix");

	commentator().start("test on BlasMatrix");
	BlasMatrix<Field> B(F, m, n);
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j)
			B.setEntry(i, j, F.init(s, i*j));
	pass = pass and testTransposeBlackbox(B);
	pass = pass and testTransposeMatrix(B);
	commentator().stop("test on BlasMatrix");

	commentator().start("test on TriplesBB");
	SparseMatrix<Field,SparseMatrixFormat::TPL> C(F, m, n);
	for (size_t i = 0; i < min(m, n); ++i) C.setEntry(i, i, F.init(s, i+1));
	pass = pass and testTransposeBlackbox(C);
	pass = pass and testTransposeMatrix(C);
	commentator().stop("test on TriplesBB");

	commentator().stop("transpose black box test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

