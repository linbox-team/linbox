
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
#include <iostream>

#include "linbox-config.h"
#include "linbox/util/commentator.h"
#include "givaro/modular-double.h"
#include "test-blackbox.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/blackbox/fibb.h"
#include "linbox/blackbox/diagonal.h"
//#include "linbox/matrix/permutation-matrix.h"
#include "linbox/blackbox/permutation.h"
//#include "linbox/blackbox/triangular.h"


using namespace LinBox;

template<class Field>
bool testFibb(FIBB<Field>& A) 
{
/*
Verify solve with apply (first provide consistent RHS).
check NSR in NS (no check of randomness) and NSB 
*/
	typedef typename Field::Element Element;
	commentator().start("testFibb", "fibb");
	ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	bool pass = true, trial;
	typedef BlasSubmatrix<BlasMatrix<Field> > Matrix;
	const Field& F(A.field());
	BlasMatrixDomain<Field> MD(F);
	size_t m = A.rowdim();
	report << "rowdim " << m << std::endl;
	size_t n = A.coldim();
	report << "coldim " << n << std::endl;
	size_t r = A.rank(r);
	report << "rank " << r << std::endl;
	Element d; F.init(d, 1);
	F.write(report << "det ", d) << std::endl;
	if (m == n)
	{	typename Field::Element d; F.init(d);
		if (r < m and not F.isZero(A.det(d))) 
		{	report << "problem in  det" << std::endl;
			pass = false;
		}
	}
	if (pass) report << "dims, rank, and det run fine (but values were not checked)" << std::endl;

	/* right side stuff */
	{
	// Solve consistent  
	size_t b = (n/2 > 5) ? 5 : n/2;
	BlasMatrix<Field> Bn(F,m,b), Xn(F, n, b), Yn(F, m, b), Zn(F, n, b);
	Matrix B(Bn), X(Xn), Y(Yn), Z(Zn);
	X.random();
	A.applyRight(B, X); // B: B = AX
	A.solveRight(Z, B);  // Z: AZ = B = AX
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
	BlasMatrix<Field> N(F);
	A.nullspaceBasisRight(N);
	trial = N.rowdim() == n and N.coldim() == n-r; // proper num of cols.
	if (n != r)
	{	trial = trial and n-r == MD.rank(N); // indep cols.
		BlasMatrix<Field> Nim_base(F, N.rowdim(), N.coldim()); 
		Matrix Ns(N), Nim(Nim_base);
		A.applyRight(Nim, Ns);
		trial = trial and MD.isZero(Nim); // in nullsp.
	}
	pass = pass and trial;
	if (not trial) report << "problem in NSB Right " << m << " " << n << " " << r << " " << N.rowdim() << " " << N.coldim() << std::endl;
	if (trial) report << " no problem in NSB Right" << std::endl;
	}
	report << " done with right side" << std::endl;

	/* left side stuff */
	{
	// Solve consistent
	size_t b = (n/2 > 6) ? 6 : n/2;
	BlasMatrix<Field> Bn(F,b,n), Xn(F, b, m), Yn(F, b, n), Zn(F, b, m);
	Matrix B(Bn), X(Xn), Y(Yn), Z(Zn);
	X.random();
	A.applyLeft(B, X); // B: B = XA
	A.solveLeft(Z, B);  // Z: ZA = B = XA
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
	BlasMatrix<Field> N(F);
	A.nullspaceBasisLeft(N);
	trial = (N.rowdim() == m-r and N.coldim() == m); // proper num of rows.
	// todo: this rank should work on empty matrix, but doesn't
	if (m != r)
	{	trial = (trial and m-r == MD.rank(N)); // indep cols.
		BlasMatrix<Field> Nim_base(F, N.rowdim(), N.coldim()); 
		Matrix Ns(N), Nim(Nim_base);
		A.applyLeft(Nim, Ns);
		trial = trial and MD.isZero(Nim); // in nullsp.
	}
	pass = pass and trial;
	if (not trial) report << "problem in NSB Left" << std::endl;
	if (trial) report << " no problem in NSB Left" << std::endl;
	}

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

	Diagonal<Field> D1(F, n, true); // nonsing
	Diagonal<Field> D2(F, n, true); 
	for (size_t i = 0; i < n/2; ++i) D2.setEntry(i,i,F.zero); // rank half
	Diagonal<Field> D3(F, n, true); 
	for (size_t i = 0; i < n; ++i) D3.setEntry(i,i,F.zero); // zero matrix

	pass = pass and testFibb(D1); // nonsing
	pass = pass and testFibb(D2); // sing
	pass = pass and testFibb(D3); // zero

	FIBBProduct<Field> Pr1(D1, D1); // nonsing product
	FIBBProduct<Field> Pr2(D1, D2); // sing product
	FIBBProduct<Field> Pr3(D3, D1); // zero product

	pass = pass and testFibb(Pr1); // nonsing product
	pass = pass and testFibb(Pr2); // sing product
	pass = pass and testFibb(Pr3); // zero product

#if 0
	Permutation<Field> P1(F, n, n); // ident
	Permutation<Field> P2(F, n, n); P2.random(); 

	pass = pass and testFibb(P1); // ident
	pass = pass and testFibb(P2); // random perm

	FIBBProduct<Field> Pr4(P1, D1, P2); // nonsing product
	FIBBProduct<Field> Pr5(P1, D2, P2); // sing product
	FIBBProduct<Field> Pr6(P1, D3, P2); // zero product 

	pass = pass and testFibb(Pr4); // nonsing product
	pass = pass and testFibb(Pr5); // sing product
	pass = pass and testFibb(Pr6); // zero product

#endif
	report << "Done with diag products" << std::endl;
#if 0
	std::cout << "Triangular" << std::endl;
	BlasMatrix<Field> M(F, 3, 3); 
	//M.random();
	for (size_t i = 0; i < 3; ++i)
	for (size_t j = 0; j < 3; ++j)
	if ( i == 0 or i == j) M.setEntry(i, j, F.one);
	M.write(report << "base matrix " << std::endl) << std::endl;
	BlasMatrixDomain<Field> BMD(F);
	size_t r; 
	BlasSubmatrix<BlasMatrix<Field> > SM(M);
	r = BMD.rank(SM);
	report << "underlying rank " << r << std::endl;
	TriangularBlasMatrix<Field> U(M, Tag::Shape::Upper, Tag::Diag::NonUnit); 
	TriangularBlasMatrix<Field> L(M, Tag::Shape::Lower, Tag::Diag::Unit); 
	report << "Upper " << (int)Tag::Shape::Upper << ", Lower " << (int)Tag::Shape::Lower << std::endl;
	report << "Unit " << (int)Tag::Diag::Unit << ", NonUnit " << (int)Tag::Diag::NonUnit << std::endl;
	Triangular<Field> UU(U);
	Triangular<Field> LL(L);
	FIBBProduct<Field> Pr7(LL, UU); // nonsing product
	FIBBProduct<Field> Pr8(P1, LL, UU, P2); // nonsing product
	FIBBProduct<Field> Pr9(LL, P1, UU, P2); // nonsing product

	report << "Triangular, " << (int)U.getDiag() << " " << (int)U.getUpLo() << std::endl;
	pass = pass and testFibb(UU); // Upper NonUnit
	report << "Triangular, " << (int)L.getDiag() << " " << (int)L.getUpLo() << std::endl;
	pass = pass and testFibb(LL); // Lower Unit
	report << "LU pattern nonsing" << std::endl;
	pass = pass and testFibb(Pr7); // LU pattern nonsing
	report << "PLUQ pattern nonsing" << std::endl;
	pass = pass and testFibb(Pr8); // PLUQ pattern nonsing
	report << "LQUP pattern nonsing" << std::endl;
	pass = pass and testFibb(Pr9); // LQUP pattern nonsing
#endif
	//FactorizedMatrix<Field> F(...);

	commentator().stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
