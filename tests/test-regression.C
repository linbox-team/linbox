/* tests/test-regression.C
 * Copyright (C) LinBox
 * Written by Clement Pernet
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

/*! @file  tests/test-regression.C
 * @ingroup tests
 * @brief tests former bugs to check that no regression made them show up again.
 */
#include "linbox/linbox-config.h"
#include "givaro/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/ring/polynomial-ring.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/solutions/solve.h"
#include "linbox/solutions/charpoly.h"

using namespace LinBox;
typedef Givaro::ZRing<Givaro::Integer> ZRingInts;

bool writing=false;

#include "linbox/algorithms/smith-form-sparseelim-poweroftwo.h"
#include "test-smith-form.h"

bool testSolveSparse(){

    typedef Givaro::Modular<unsigned int> Field;
    typedef std::vector <std::pair <size_t, Field::Element> > SparseSeqVectorGFp;
    typedef SparseMatrix<Field, VectorTraits<SparseSeqVectorGFp>::SparseFormat> Matrix;
    Field F(127);
    Matrix A(F, 3, 3);
    DenseVector<Field> x(F, 3);
    DenseVector<Field> b(F, 3);
    A.setEntry(0,0, 1);
    A.setEntry(0,1, 2);
    A.setEntry(0,2, 3);
    A.setEntry(1,0, 126);
    A.setEntry(1,1, 2);
    A.setEntry(1,2, 5);
    A.setEntry(2,0, 2);
    A.setEntry(2,1, 3);
    A.setEntry(2,2, 1);
    for (int i=0;i<3;i++)
        F.assign(b[i],i+1);

    Matrix B(A);

    solve(x,B,b,  Method::DenseElimination());

    if (!F.areEqual (x[0],73)) return false;
    if (!F.areEqual (x[1],76)) return false;
    if (!F.areEqual (x[2],10)) return false;

	if (writing) std::cout << "TSSF: PASSED" << std::endl;

    return true;
}

#if __LINBOX_HAVE_SAGE
// CP: Include the .C file instead of linking to the linboxsage lib:
// avoid requiring to make install before make check.
#include "interfaces/sage/linbox-sage.C"
bool testSolveSparseSage(){
    size_t p = 127;
    c_vector_modint_linbox * A = new c_vector_modint_linbox[3];
    c_vector_modint_linbox  b;
    for (int i=0;i<3;++i){
        A[i].entries = new int[3];
        A[i].positions = new size_t[3];
        A[i].num_nonzero=3;
        for (int j=0;j<3;j++)
            A[i].positions[j]=j;
        A[i].p = p;
    }
    b.entries = new int[3];
    b.positions = new size_t[3];
    for (int j=0;j<3;j++)
        b.positions[j]=j;
    b.num_nonzero=3;
    b.p = p;
    for (int i=0;i<3;i++)
        b.entries[i]=i+1;

    std::vector<unsigned int> x;

    A[0].entries[0]=1;
    A[0].entries[1]=2;
    A[0].entries[2]=3;
    A[1].entries[0]=126;
    A[1].entries[1]=2;
    A[1].entries[2]=5;
    A[2].entries[0]=2;
    A[2].entries[1]=3;
    A[2].entries[2]=1;

    x = linbox_modn_sparse_matrix_solve(p,3,3,A,&b,1);

    if (x[0] != 73) return false;
    if (x[1] != 76) return false;
    if (x[2] != 10) return false;

	if (writing) std::cout << "TSSS: PASSED" << std::endl;

    return true;
}
#else
bool testSolveSparseSage(){return true;}
#endif

/**
 * Testing regresssion for issue 56 https://github.com/linbox-team/linbox/issues/56
 * reported by Vincent Delecroix
 */
template <class SolveMethod>
bool testFlatDixonSolver (const SolveMethod& m){
        // creating LinBox matrices and vectors
    Givaro::ZRing<Integer> ZZ;
    typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;

    SparseMatrix<Givaro::ZRing<Integer> > M (ZZ,1,2);
    Givaro::ZRing<Integer>::Element D; // denominator of the solution
    DenseVector A(ZZ, M.coldim()); // numerator of the solution
    DenseVector B(ZZ, M.rowdim()); // Right handside of the system

    M.setEntry(0,0,1);
    M.setEntry(0,1,1);
    ZZ.assign(B[0],1);

        // solving via Dixon lifting
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.one) || !ZZ.areEqual(A[1],ZZ.zero) || !ZZ.areEqual(D,ZZ.one)) {
        if (writing) std::cerr<<"**** ERROR **** Fail solving a flat system over QQ via Dixon Lifting"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TFDS: PASSED" << std::endl;


    return true;
}

template <class SolveMethod>
bool testFlatDixonSolver2 (const SolveMethod& m){
        // creating LinBox matrices and vectors
    Givaro::ZRing<Integer> ZZ;
    typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;

    SparseMatrix<Givaro::ZRing<Integer> > M (ZZ,1,2);
    Givaro::ZRing<Integer>::Element D; // denominator of the solution
    DenseVector A(ZZ, M.coldim()); // numerator of the solution
    DenseVector B(ZZ, M.rowdim()); // Right handside of the system

    M.setEntry(0,0,1);
    M.setEntry(0,1,0);
    ZZ.assign(B[0],1);

        // solving via Dixon Lifting
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.one) || !ZZ.areEqual(D,ZZ.one)) {
        if (writing) std::cerr<<"A = "<<A<<" D = "<<D<<std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail solving a flat system over QQ via Dixon Lifting"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TFD2: PASSED" << std::endl;


    return true;
}

template <class SolveMethod>
bool testTallDixonSolver (const SolveMethod& m){
        // creating LinBox matrices and vectors
    Givaro::ZRing<Integer> ZZ;
    typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;

    SparseMatrix<Givaro::ZRing<Integer> > M (ZZ,2,1);
    Givaro::ZRing<Integer>::Element D; // denominator of the solution
    DenseVector A(ZZ, M.coldim()); // numerator of the solution
    DenseVector B(ZZ, M.rowdim()); // Right handside of the system

    M.setEntry(0,0,1);
    M.setEntry(1,0,1);
    ZZ.assign(B[0],1);
    ZZ.assign(B[1],1);

        // solving via Dixon Lifting
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.one) || !ZZ.areEqual(D,ZZ.one)) {
        if (writing) std::cerr<<"**** ERROR **** Fail solving a tall system over QQ via Dixon Lifting"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TTDS: PASSED" << std::endl;


    return true;

}

template <class SolveMethod>
bool testSingularDixonSolver (const SolveMethod& m){
        // creating LinBox matrices and vectors
    Givaro::ZRing<Integer> ZZ;
    typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;

    SparseMatrix<Givaro::ZRing<Integer> > M (ZZ,2,2);
    Givaro::ZRing<Integer>::Element D; // denominator of the solution
    DenseVector A(ZZ, M.coldim()); // numerator of the solution
    DenseVector B(ZZ, M.rowdim()); // Right handside of the system

    M.setEntry(0,0,-2);
    ZZ.assign(B[0],-4);

        // solving via Sparse Elimination
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],Integer(2)) || !ZZ.areEqual(D,ZZ.one)) {
        if (writing) std::cerr<<"A = "<<A<<" D = "<<D<<std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail solving a singular system over QQ with a SparseMatrix"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TSDS: PASSED" << std::endl;

    return true;
}

template <class SolveMethod>
bool testZeroDixonSolver (const SolveMethod& m){
        // creating LinBox matrices and vectors
    Givaro::ZRing<Integer> ZZ;
    typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;

    SparseMatrix<Givaro::ZRing<Integer> > M (ZZ,1,1);
    Givaro::ZRing<Integer>::Element D; // denominator of the solution
    DenseVector A(ZZ, M.coldim()); // numerator of the solution
    DenseVector B(ZZ, M.rowdim()); // Right handside of the system

    M.setEntry(0,0,0);
    ZZ.assign(B[0],0);

        // solving via Dixon Lifting
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.zero) || !ZZ.areEqual(D,ZZ.one)) {
        if (writing) std::cerr<<"A = "<<A<<" D = "<<D<<std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail solving a zero over QQ via Dixon Lifting"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TZDS: PASSED" << std::endl;

    return true;
}



// test for bug #107 from Zhu
bool testDixonSolverWithMPrhs (){

    // creating LinBox matrices and vectors
    Givaro::ZRing<Integer> ZZ;
    typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;

    DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,1,1);
    Givaro::ZRing<Integer>::Element D; // denominator of the solution
    DenseVector X(ZZ, A.coldim()); // numerator of the solution
    DenseVector B(ZZ, A.rowdim()); // Right handside of the system

    A.setEntry(0,0,Integer("12345678901234567890"));
    ZZ.assign(B[0],Integer("12345678901234567890"));

    // solving via Dixon lifting
    solve (X, D, A, B, Method::DenseElimination());

    if (!ZZ.areEqual(X[0],ZZ.one) ||  !ZZ.areEqual(D,ZZ.one)) {
        if (writing) std::cerr<<"**** ERROR **** Fail solving a system over QQ with a DenseMatrix and a MP rhs"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TDSM: PASSED" << std::endl;


    return true;
}

bool testSparseRationalSolver() {
    typedef Givaro::QField<Givaro::Rational> Rats;
    Rats QQ;
    typedef DenseVector<Rats> RVector;
    SparseMatrix<Rats> A (QQ,1,3);
    RVector X(QQ, A.coldim()),B(QQ, A.rowdim()),L(QQ, A.rowdim());
    A.setEntry(0,1,1);
    A.setEntry(0,2,2);
    QQ.assign(B[0],1);

        // Directly solve of Q
    solve(X,A,B,Method::SparseElimination());

    MatrixDomain<Rats> MD(QQ);
    VectorDomain<Rats> VD(QQ);

    MD.vectorMul(L, A, X);

    if (! VD.areEqual(L, B)) {
        if (writing) A.write(std::cerr << "A:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
        if (writing) std::cerr<<"X:= "<< X << ';' << std::endl;
        if (writing) std::cerr<<"B:= "<< B << ';' << std::endl;
        if (writing) std::cerr<<"L:= "<< L << ';' << std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail solving a sparse system over QQ"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TSRS: PASSED" << std::endl;


    return true;
}

template<class SolveMethod, typename Matrix_t=SparseMatrix<ZRingInts>>
bool testDixonRectangularSolver() {
    SolveMethod m;

    ZRingInts ZZ;
    typedef DenseVector<ZRingInts> RVector;
    Matrix_t A (ZZ,1,3);
    RVector X(ZZ, A.coldim()),B(ZZ, A.rowdim()),L(ZZ, A.rowdim());
    ZRingInts::Element d;

    A.setEntry(0,1,1);
    A.setEntry(0,2,2);
    ZZ.assign(B[0],1);

        // Dixon Lifting
    solve(X,d,A,B,m);

    bool pass=true;

    if (ZZ.isZero(d))
        pass = false;
    else{
        MatrixDomain<ZRingInts> MD(ZZ);
        VectorDomain<ZRingInts> VD(ZZ);

        MD.vectorMul(L, A, X);
        VD.mulin(B, d);

        if (! VD.areEqual(L, B)) {
            pass = false;
        }
    }

    if (! pass) {
        if (writing) A.write(std::cerr << "A:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
        if (writing) std::cerr<<"X:= "<< X << ';' << std::endl;
        if (writing) std::cerr<<"d:= "<< d << ';' << std::endl;
        if (writing) std::cerr<<"d * B:= "<< B << ';' << std::endl;
        if (writing) std::cerr<<"L:= "<< L << ';' << std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail solving a sparse system over ZZ" <<std::endl;
        return false;
    } else
        if (writing) std::cout << "TDRS: PASSED" << std::endl;
    return true;
}

bool testSparseDiagDet(uint64_t n){
    ZRingInts ZZ;
    SparseMatrix<ZRingInts> A(ZZ, n, n);
    Integer x(2);
    for(size_t i=0; i<n; ++i)
        A.setEntry(i,i,x);
    Integer d, r(Integer(1)<<n);

    if (! ZZ.areEqual(r, det(d,A)) ) {
        if (writing) A.write(std::cerr << "A:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
        if (writing) std::cerr<<"r:= "<< r << ';' << std::endl;
        if (writing) std::cerr<<"d:= "<< d << ';' << std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail determinant 2 times Identity matrix" <<std::endl;
        return false;
    } else
        if (writing) std::cout << "TDD2: PASSED" << std::endl;

    return true;
}

bool testSparse1x1Det(uint64_t v){
    ZRingInts ZZ;
    SparseMatrix<ZRingInts> A(ZZ, 1, 1);
    Integer x(v);
    A.setEntry(0,0,x);
    Integer d;

    if (! ZZ.areEqual(x, det(d,A)) ) {
        if (writing) A.write(std::cerr << "A:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
        if (writing) std::cerr<<"x:= "<< x << ';' << std::endl;
        if (writing) std::cerr<<"d:= "<< d << ';' << std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail determinant 1x1 matrix" <<std::endl;
        return false;
    } else
        if (writing) std::cout << "TSD1: PASSED" << std::endl;

    return true;
}

bool testZeroDimensionalCharPoly(){
    Givaro::ZRing<Integer> ZZ;
    DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,0,0);
    DensePolynomial<Givaro::ZRing<Integer> > P (ZZ);
    PolynomialRing<Givaro::ZRing<Integer> > PR (ZZ);
    charpoly(P,A);
    return PR.isOne(P);
}

bool testZeroDimensionalMinPoly(){
    Givaro::ZRing<Integer> ZZ;
    DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,0,0);
    DensePolynomial<Givaro::ZRing<Integer> > P (ZZ);
    PolynomialRing<Givaro::ZRing<Integer> > PR (ZZ);
    minpoly(P,A);
    return PR.isOne(P);
}

bool testBigScalarCharPoly(){
    Givaro::ZRing<Integer> ZZ;
    DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,1,1);
    Integer x;
    ZZ.init(x);
    ZZ.assign(x,ZZ.one);
    x <<= 2000;
    A.setEntry(0,0,x);
    DensePolynomial<Givaro::ZRing<Integer> > P (ZZ);
    DensePolynomial<Givaro::ZRing<Integer> > Q (ZZ,2);
    PolynomialRing<Givaro::ZRing<Integer> > PR (ZZ);
    charpoly(P,A);
    ZZ.assign(Q[1],ZZ.one);
    ZZ.neg(Q[0],x);
    return PR.areEqual(P,Q);
}


template<typename Matrix_t=SparseMatrix<ZRingInts>, class SolveMethod>
bool testInconsistent (const SolveMethod& m){
        // creating LinBox matrices and vectors
    Givaro::ZRing<Integer> ZZ;
    typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;

    Matrix_t M (ZZ,1,1);
    Givaro::ZRing<Integer>::Element D; // denominator of the solution
    DenseVector A(ZZ, M.coldim()); // numerator of the solution
    DenseVector B(ZZ, M.rowdim()); // Right handside of the system

    M.setEntry(0,0,0);
    ZZ.assign(B[0],ZZ.one);

        // solving via Dixon Lifting
    try{
        solve (A, D, M, B, m);
    }
    catch (LinBox::LinboxMathInconsistentSystem& e){return true;}

    return false;
}


bool testLocalSmith(){
    typedef Givaro::ZRing<uint64_t> Ring; // characteristic 2
    typedef LinBox::SparseMatrix<Ring,
        LinBox::SparseMatrixFormat::SparseSeq > SparseMat;

    SmithList<Ring> local;
    Ring R;
    SparseMat A(R,2,3);
    A.setEntry(0,0, 2);
    A.setEntry(0,2, 1);
    A.setEntry(1,0, 2);

    LinBox::PowerGaussDomainPowerOfTwo< uint64_t > PGD;
    LinBox::GF2 F2;
    Permutation<GF2> Q(F2,A.coldim());

    PGD(local, A, Q, 5, PRESERVE_UPPER_MATRIX|PRIVILEGIATE_NO_COLUMN_PIVOTING);

    if (writing) {
        std::clog << "Local Smith: {";
        for(auto const& it:local) std::clog << '{' << it.first << ',' << it.second << '}';
        std::clog << '}';
    }

// > ([1,1] [1,2] )
// > [[1, 2, 0 ], [0, 1, 0 ]]
// > [[0,0,1], [1,0,0], [0,1,0]]

        // Smith form
    SmithList<Ring> correctSL{ {1U,1U},{2U,1U} };
    bool success( checkSNFExample(correctSL, local, R) );

    if (writing)
        A.write(std::clog << ", A:=", LinBox::Tag::FileFormat::Maple) << ';';

        // Upper triangular
    success &=
        (A[0].size() == 2) &&
        (A[0][0].first == 0U) &&
        (A[0][0].second == 1U) &&
        (A[0][1].first == 1U) &&
        (A[0][1].second == 2U) &&
        (A[1].size() == 1) &&
        (A[1][0].first == 1U) &&
        (A[1][0].second == 1U);

    if (writing)
        Q.write(std::clog << ", Q:=", LinBox::Tag::FileFormat::Maple) << ';';

        // Permutation
    success &=
        (Q[0] == 2) &&
        (Q[1] == 0) &&
        (Q[2] == 1);

    if (!success) {
        if (writing) std::clog<<"**** ERROR **** Fail Local Smith" <<std::endl;
        return false;
    } else
        if (writing) std::cout << "TSLS: PASSED" << std::endl;
    return success;
}

bool testDixonSmallFat() {
    bool success;
    using Ring = Givaro::ZRing<Integer>;
    using Matrix = DenseMatrix<Ring>;
    using Vector = DenseVector<Ring>;

    Ring ZZ;

    Matrix A(ZZ, 1, 3);
    Vector x(ZZ, A.coldim());
    Vector b(ZZ, A.rowdim());
    Vector r(ZZ, A.rowdim());
    Ring::Element d;

    A.setEntry(0, 1, 1);
    A.setEntry(0, 2, 2);
    ZZ.assign(b[0], 1);

    // Calling Dixon
    solve(x, d, A, b, Method::Dixon(SingularSolutionType::Deterministic) );

    A.apply(r,x);
    ZZ.divin(r[0],d);
    return success = ZZ.areEqual(r[0], b[0]);
}

template<typename CMethod=Method::Blackbox>
bool testZeroMatrixCharPoly(const CMethod& cMeth=CMethod()) {
    bool success;
	using Ring = Givaro::Modular<double>;
	using Matrix = SparseMatrix<Ring>;
    Ring R(3);

    Matrix A(R, 1, 1);
    A.setEntry(0, 0, R.zero);

    PolynomialRing<Ring>::Element c_A, Ex;

    charpoly(c_A, A, RingCategories::ModularTag(), cMeth );

    PolynomialRing<Ring> PZ(R,'X'); PZ.assign(Ex, Givaro::Degree(1), R.one);

    success = PZ.areEqual(c_A, Ex);

    if (!success) {
        if (writing) {
            std::clog<<"**** ERROR **** Fail ZMCP " <<std::endl;

            PZ.write(std::clog << "Ex: ", Ex) << std::endl;
            PZ.write(std::clog << "cA: ", c_A) << std::endl;
        }
        return false;
    } else
        if (writing) std::cout << "ZMCP" << cMeth.name() << cMeth.shapeFlags.flags << " : PASSED" << std::endl;

    return success;
}

bool testFourFourMatrix() {
    bool success;
    using Ring = Givaro::ZRing<Integer>;
    using Matrix = DenseMatrix<Ring>;
    Ring ZZ;

    Matrix A(ZZ, 4,4);
    for(size_t i=0; i<4; ++i) for(size_t j=0; j<4; ++j)
        A.setEntry(i,j, static_cast<uint64_t>(4*i+j+1));

    PolynomialRing<Ring>::Element c_A, Res;

    charpoly(c_A, A);

    PolynomialRing<Ring> PZ(ZZ,'X');
    PZ.assign(Res, Givaro::Degree(4), ZZ.one);
    Res[2] = -80;
    Res[3] = -34;

    success = PZ.areEqual(c_A, Res);

    if (!success) {
        if (writing) {
            std::clog<<"**** ERROR **** Fail tFFM " <<std::endl;

            PZ.write(std::clog << "Ex: ", Res) << std::endl;
            PZ.write(std::clog << "cA: ", c_A) << std::endl;
        }

        return false;
    } else
        if (writing) std::cout << "tFFM: PASSED" << std::endl;

    return success;
}

int main (int argc, char **argv)
{
    bool pass = true;

	// text is written to clog/cerr/cout iff a command line argument is present.
	if (argc > 1) writing = true;

    if (writing) {
        commentator().setReportStream(std::cout);
    }

    pass &= testSolveSparse  ();
    pass &= testSolveSparseSage ();
    pass &= testFlatDixonSolver (Method::SparseElimination());
    pass &= testFlatDixonSolver2 (Method::SparseElimination());
    pass &= testTallDixonSolver (Method::SparseElimination());
    pass &= testFlatDixonSolver (Method::DenseElimination());
    pass &= testFlatDixonSolver2 (Method::DenseElimination());
    pass &= testTallDixonSolver (Method::DenseElimination());
    pass &= testFlatDixonSolver (Method::Wiedemann());
    pass &= testFlatDixonSolver2 (Method::Wiedemann());
    pass &= testTallDixonSolver (Method::Wiedemann());
    pass &= testSingularDixonSolver (Method::SparseElimination());
    pass &= testZeroDixonSolver (Method::SparseElimination());
    pass &= testSingularDixonSolver (Method::DenseElimination());
    pass &= testZeroDixonSolver (Method::DenseElimination());
    pass &= testDixonSolverWithMPrhs ();
    pass &= testSparseRationalSolver ();
    pass &= testDixonRectangularSolver<Method::DenseElimination> ();
    pass &= testDixonRectangularSolver<Method::SparseElimination> ();
    pass &= testDixonRectangularSolver<Method::Wiedemann> ();
    pass &= testDixonRectangularSolver<Method::DenseElimination, DenseMatrix<ZRingInts>> ();
    pass &= testDixonRectangularSolver<Method::SparseElimination, DenseMatrix<ZRingInts>> ();
    pass &= testDixonRectangularSolver<Method::Wiedemann, DenseMatrix<ZRingInts>> ();
    pass &= testSparse1x1Det(1<<26);
    pass &= testSparseDiagDet(46);
    pass &= testZeroDimensionalCharPoly ();
    pass &= testZeroDimensionalMinPoly ();
    pass &= testZeroMatrixCharPoly<>();
    pass &= testZeroMatrixCharPoly(Method::Blackbox(Shape::Symmetric));
    pass &= testZeroMatrixCharPoly(Method::DenseElimination());
//     pass &= testZeroMatrixCharPoly(Method::SparseElimination());
    pass &= testBigScalarCharPoly ();
    pass &= testLocalSmith ();
    pass &= testInconsistent<DenseMatrix<ZRingInts>> (Method::DenseElimination());
    pass &= testDixonSmallFat();
    pass &= testFourFourMatrix();

        // Still failing: see https://github.com/linbox-team/linbox/issues/105
        //pass &= testInconsistent<> (Method::SparseElimination());
        //pass &= testInconsistent<> (Method::Wiedemann());

    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
