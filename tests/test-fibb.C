/* tests/test-fibb.C
 * -bds
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
 *
 */


/*! @file  tests/test-fibb.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */

#include "linbox/linbox-config.h"
#include <iostream>
#include "linbox/util/commentator.h"
#include "givaro/modular-floating.h"
#include "test-blackbox.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/blackbox/fibb.h"
#include "linbox/blackbox/fibb-product.h"
#include "linbox/blackbox/diagonal.h"
//#include "linbox/matrix/permutation-matrix.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/blackbox/triangular-fibb.h"

using namespace LinBox;

/*
check rk = rank(A), dt = det(A)
Neg rk argument means we don't know rank, 
Zero dt with full rank  means we don't know det.
In these cases the test passes without checking rank,det.

Verify solve with apply (first provide consistent RHS).
check NSR in NS (no check of randomness) and NSB 
*/
template<class Field>
bool testFibb(FIBB<Field>& A, string title, const typename Field::Element& dt, int64_t rk = -1) 
{
	typedef typename Field::Element Element;
	commentator().start(title.c_str(), "fibb");
	ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	bool pass = true, trial;
	typedef DenseMatrix<Field> Matrix;
	const Field& F(A.field());
	BlasMatrixDomain<Field> MD(F);
	size_t m = A.rowdim();
	report << "rowdim " << m << std::endl;
	size_t n = A.coldim();
	report << "coldim " << n << std::endl;
	size_t r = A.rank(r);
	if (rk >= 0u) pass = (pass and (int64_t(r) == rk));
	if (m < 8u  and n < 8u) A.write(report << "A: " << std::endl) << std::endl;
	report << "rank " << r << ", expected " << rk << std::endl;

	Element d; F.init(d);
	A.det(d);
	if (rk == int64_t(n) and rk == int64_t(m) and F.isZero(dt))
		F.write(report << "det ", d) << ", unchecked" << std::endl;
	else 
		F.write(F.write(report << "det ", d) << ", expected ", dt) << std::endl;
	if (rk == int64_t(n) and not F.isZero(dt)) pass = pass and F.areEqual(d, dt);
	if (m == n)
	{	typename Field::Element d; F.init(d);
		if (r < m and not F.isZero(A.det(d))) 
		{	report << "problem in rank or det" << std::endl;
			pass = false;
		}
	}
	if (pass) report << "dims, rank, and det run fine (values partially checked)" << std::endl;

/* right side stuff */
{
#if 1
	// Solve consistent  
	size_t b = (n/2 > 5) ? 5 : n/2;
	Matrix B(F,m,b), X(F, n, b), Y(F, m, b), Z(F, n, b);
	X.random();
	A.applyRight(B, X); // B: B = AX
	if (not F.isZero(dt)) { // check for zero apply
		if (MD.isZero(X)) { 
			X.write(report << "unlucky X" << std::endl, Tag::FileFormat::Plain) << std::endl;
			pass = false;
		}
		if (MD.isZero(B)) { 
			B.write(report << "bad apply, B" << std::endl, Tag::FileFormat::Plain) << std::endl;
			pass = false;
		}
	}
	A.solveRight(Z, B); // Z: AZ = B = AX
	A.applyRight(Y, Z); // Y: Y = AZ = B
	trial = MD.areEqual(Y, B);
	if (m == n and m == r) trial = trial and MD.areEqual(Z,X); // nonsing case
	pass = pass and trial;
	if (not trial) report << "problem in solve Right" << std::endl;
	if (trial) report << " no problem in solve Right" << std::endl;
	// NSR
	A.nullspaceRandomRight(Y);
	A.applyRight(Z, Y);
	trial = MD.isZero(Z);
	// no check of randomness of cols of Y
	pass = pass and trial;
	if (not trial) report << "problem in NSR Right" << std::endl;
	if (trial) report << " no problem in NSR Right" << std::endl;
	// NSB
	Matrix N(F);
	A.nullspaceBasisRight(N);
	trial = N.rowdim() == n and N.coldim() == n-r; // proper num of cols.
	if (n != r)
	{	trial = trial and n-r == MD.rank(N); // indep cols.
		Matrix Nim(F, N.rowdim(), N.coldim()); 
		A.applyRight(Nim, N);
		trial = trial and MD.isZero(Nim); // in nullsp.
		if (n < 8) Nim.write(report << "Nim " << std::endl, Tag::FileFormat::Plain) << std::endl;
	}
	if (n < 8) N.write(report << "N " << std::endl, Tag::FileFormat::Plain) << std::endl;
	pass = pass and trial;
	if (not trial) report << "problem in NSB Right " << m << " " << n << " " << r << " " << N.rowdim() << " " << N.coldim() << std::endl;
	if (trial) report << " no problem in NSB Right" << std::endl;
#endif
} // right side 
	report << " done with right side" << std::endl;

/* left side stuff */
{
	// Solve consistent
	size_t b = (n/2 > 6) ? 6 : n/2;
	Matrix B(F,b,n), X(F, b, m), Y(F, b, n), Z(F, b, m);
	X.random();
	A.applyLeft(B, X); // B: B = XA
#if 1
	if (not F.isZero(dt)) { // check for zero apply
		if (MD.isZero(X)) { 
			X.write(report << "unlucky X" << std::endl, Tag::FileFormat::Plain) << std::endl;
			pass = false;
		}
		if (MD.isZero(B)) { 
			B.write(report << "bad apply, B" << std::endl, Tag::FileFormat::Plain) << std::endl;
			pass = false;
		}
	}
	A.solveLeft(Z, B); // Z: ZA = B = XA
	A.applyLeft(Y, Z); // Y: Y = ZA
	trial = MD.areEqual(Y, B);
	if (m == n and m == r) trial = trial and MD.areEqual(Z,X); // nonsing case
	pass = pass and trial;
	if (not trial) 
	{	report << "problem in solve Left" << std::endl;
	 	A.write(report << "FIBB A: " << std::endl) << std::endl;
	 	X.write(report << "X: random " << std::endl, Tag::FileFormat::Plain) << std::endl;
	 	B.write(report << "B: XA" << std::endl, Tag::FileFormat::Plain) << std::endl;
	 	Z.write(report << "Z: solve ZA = B " << std::endl, Tag::FileFormat::Plain) << std::endl;
	 	Y.write(report << "Y: ZA " << std::endl, Tag::FileFormat::Plain) << std::endl;
	}
	if (trial) report << " no problem in solve Left" << std::endl;
	// NSR
	A.nullspaceRandomLeft(X); // X:  XA = 0 
	A.applyLeft(Y, X); // Y: Y = XA
	trial = MD.isZero(Y);
	pass = pass and trial;
	if (not trial) report << "problem in NSR Left" << std::endl;
	if (trial) report << " no problem in NSR Left" << std::endl;
	// NSB
	Matrix N(F);
	A.nullspaceBasisLeft(N);
	trial = (N.rowdim() == m-r and N.coldim() == m); // proper num of rows.
	// todo: this rank should work on empty matrix, but doesn't
	if (m != r)
	{	trial = (trial and m-r == MD.rank(N)); // indep cols.
		Matrix Nim(F, N.rowdim(), N.coldim()); 
		A.applyLeft(Nim, N);
		trial = trial and MD.isZero(Nim); // in nullsp.
		if (n < 8) Nim.write(report<< "Nim " << std::endl, Tag::FileFormat::Plain);
	}
	if (n < 8) N.write(report << "N " << std::endl, Tag::FileFormat::Plain) << std::endl;
	pass = pass and trial;
	if (not trial) report << "problem in NSB Left" << std::endl;
	if (trial) report << " no problem in NSB Left" << std::endl;
#endif
} // left side
	report << " done with left side" << std::endl;

	commentator().stop (MSG_STATUS (pass));
	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;
	//srand(time(NULL));

	static size_t n = 10;
	static int32_t q = 101;
	//static uint32_t q = 2147483647U;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	commentator().start("FIBB test suite", "fibb");
	ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	typedef Givaro::Modular<double> Field; 
	Field F (q);
	//MatrixDomain<Field> MD(F);

	Diagonal<Field> D1(F, n); // random nonsing
	Diagonal<Field> D2(F, n); 
	for (size_t i = 0; i < n/2; ++i) D2.setEntry(i,i,F.zero); // rank half
	Diagonal<Field> D3(F, n); 
	for (size_t i = 0; i < n; ++i) D3.setEntry(i,i,F.zero); // zero matrix

	pass = pass and testFibb(D1, "D1 nonsing", F.zero, n); // nonsing
	pass = pass and testFibb(D2, "D2 sing", F.zero, n-n/2); // sing
	pass = pass and testFibb(D3, "D3 zero", F.zero, 0); // zero

	FIBBProduct<Field> Pr1(D1, D1); // nonsing product
	FIBBProduct<Field> Pr2(D1, D2); // sing product
	FIBBProduct<Field> Pr3(D3, D1); // zero product

	pass = pass and testFibb(Pr1, "Pr1 nonsing product", F.zero, n); // nonsing product
	pass = pass and testFibb(Pr2, "Pr2 sing product", F.zero, n-n/2); // sing product
	pass = pass and testFibb(Pr3, "Pr3 zero product", F.zero, 0); // zero product

	Permutation<Field> P1(F, n, n); // ident
	Permutation<Field> P2(F, n, n); P2.random(); 

	pass = pass and testFibb(P1, "P1 ident", F.one, n); // ident
	pass = pass and testFibb(P2, "P2 random perm", F.zero, n); // random perm

	FIBBProduct<Field> Pr4(P1, D1, P2); // nonsing product
	FIBBProduct<Field> Pr5(P1, D2, P2); // sing product
	FIBBProduct<Field> Pr6(P1, D3, P2); // zero product 

	pass = pass and testFibb(Pr4, "Pr4 nonsing product", F.zero, n); // nonsing product
	pass = pass and testFibb(Pr5, "Pr5 sing product", F.zero, n-n/2); // sing product
	pass = pass and testFibb(Pr6, "Pr6 zero product", F.zero, 0); // zero product

	report << "Done with diag and permutation products" << std::endl;
#if 1
	typedef DenseMatrix<Field> Matrix;
	typedef TriangularBlasMatrix<Matrix> TriangularMatrix;
	Matrix M(F, n, n); 
	//M.random();
	for (size_t i = 0; i < n; ++i)
	for (size_t j = 0; j < n; ++j)
	if ( i == 0 or i == j) M.setEntry(i, j, F.one);
	M.write(report << "base matrix " << std::endl) << std::endl;
	BlasMatrixDomain<Field> BMD(F);
	/*
	size_t r; 
	Matrix SM(M);
	r = BMD.rank(SM);
	report << "underlying rank " << r << std::endl;
	*/
	TriangularMatrix U(M, Tag::Shape::Upper, Tag::Diag::NonUnit); 
	TriangularMatrix L(M, Tag::Shape::Lower, Tag::Diag::Unit); 
	report << "Upper " << (int)Tag::Shape::Upper << ", Lower " << (int)Tag::Shape::Lower << std::endl;
	report << "Unit " << (int)Tag::Diag::Unit << ", NonUnit " << (int)Tag::Diag::NonUnit << std::endl;
	TriangularFIBB<Field> UU(U);
	TriangularFIBB<Field> LL(L);
	FIBBProduct<Field> Pr7(LL, UU); // nonsing product
	FIBBProduct<Field> Pr8(P1, LL, UU, P2); // nonsing product
	FIBBProduct<Field> Pr9(LL, P1, UU, P2); // nonsing product

	report << "Triangular, " << (int)U.getDiag() << " " << (int)U.getUpLo() << std::endl;
	pass = pass and testFibb(UU, "UU nonU", F.one, n); // Upper NonUnit
	report << "Triangular, " << (int)L.getDiag() << " " << (int)L.getUpLo() << std::endl;
	pass = pass and testFibb(LL, "LL unit", F.one, n); // Lower Unit
	report << "LU pattern nonsing" << std::endl;
	pass = pass and testFibb(Pr7, "Pr7 LU", F.one, n); // LU pattern nonsing
	report << "PLUQ pattern nonsing" << std::endl;
	pass = pass and testFibb(Pr8, "Pr8 PLUQ", F.one, n); // PLUQ pattern nonsing
	report << "LQUP pattern nonsing" << std::endl;
	pass = pass and testFibb(Pr9, "Pr9 LQUP", F.one, n); // LQUP pattern nonsing
#endif
	report << "Done with triangular and permutation products" << std::endl;
	//FactorizedMatrix<Field> F(...);

	commentator().stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
