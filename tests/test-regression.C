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
#include "linbox/polynomial/dense-polynomial.h"
#include "linbox/ring/polynomial-ring.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/solutions/solve.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/algorithms/smith-form-sparseelim-poweroftwo.h"
using namespace LinBox;
typedef Givaro::ZRing<Givaro::Integer> ZRingInts;

bool writing=false;

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
    
    solve(x,B,b,  Method::BlasElimination());

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
bool testFlatDixonSolver (const Specifier& m){
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


bool testFlatDixonSolver2 (const Specifier& m){
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

bool testTallDixonSolver (const Specifier& m){
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

bool testSingularDixonSolver (const Specifier& m){
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
bool testZeroDixonSolver (const Specifier& m){
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
    solve (X, D, A, B, Method::BlasElimination());

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

template<typename Matrix_t=SparseMatrix<ZRingInts>>
bool testDixonRectangularSolver(const Specifier& m) {
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

/*
template<typename Matrix_t=SparseMatrix<ZRingInts>>
bool testInconsistent (const Specifier& m){
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
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.zero) || !ZZ.areEqual(D,ZZ.one)) {
        if (writing) std::cerr<<"A = "<<A<<" D = "<<D<<std::endl;
        if (writing) std::cerr<<"**** ERROR **** Fail solving an inconsistent system  over QQ via Dixon Lifting"<<std::endl;
        return false;
    } else
        if (writing) std::cout << "TICS: PASSED" << std::endl;

    return true;
}
*/


bool testLocalSmith(){
    typedef Givaro::ZRing<int64_t> Ring;
    typedef std::vector<std::pair<size_t,uint64_t> > Smith_t;
    typedef LinBox::SparseMatrix<Ring,
        LinBox::SparseMatrixFormat::SparseSeq > SparseMat;

    Smith_t local;
    Ring R;
    SparseMat A(R,2,3);
    A.setEntry(0,0, 2);
    A.setEntry(0,2, 1);
    A.setEntry(1,0, 2);

    LinBox::PowerGaussDomainPowerOfTwo< uint64_t > PGD;
    LinBox::GF2 F2;
    Permutation<GF2> Q(F2,A.coldim());

    PGD(local, A, Q, 5, PRESERVE_UPPER_MATRIX|PRIVILEGIATE_NO_COLUMN_PIVOTING);

    if (writing) std::cerr << "Local Smith: {";
    for(auto const& it:local) if (writing) std::cerr << it.first << ':' << it.second << ' ';
    if (writing) std::cerr << '}';

// > ([1,1] [1,2] )
// > [[1, 2, 0 ], [0, 1, 0 ]]
// > [[0,0,1], [1,0,0], [0,1,0]]

        // Smith form
    bool success =
        (local.size() == 2U) &&
        (local[0].first == 1U) &&
        (local[0].second == 1U) &&
        (local[1].first == 1U) &&
        (local[1].second == 2U) ;

    if (writing) A.write(std::cerr << ", A:=", LinBox::Tag::FileFormat::Maple) << ';';

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

    if (writing) Q.write(std::cerr << ", Q:=", LinBox::Tag::FileFormat::Maple) << ';';

        // Permutation
    success &=
        (Q[0] == 2) &&
        (Q[1] == 0) &&
        (Q[2] == 1);

    if (!success) {
        if (writing) std::cerr<<"**** ERROR **** Fail Local Smith" <<std::endl;
        return false;
    } else
        if (writing) std::cout << "TSLS: PASSED" << std::endl;
     return success;
}



int main (int argc, char **argv)
{
    bool pass = true;

	// text is written to cerr/cout iff a command line argument is present.
	if (argc > 1) writing = true;

    pass &= testSolveSparse  ();
    pass &= testSolveSparseSage ();
    pass &= testFlatDixonSolver (Method::SparseElimination());
    pass &= testFlatDixonSolver2 (Method::SparseElimination());
    pass &= testTallDixonSolver (Method::SparseElimination());
    pass &= testFlatDixonSolver (Method::BlasElimination());
    pass &= testFlatDixonSolver2 (Method::BlasElimination());
    pass &= testTallDixonSolver (Method::BlasElimination());
    pass &= testFlatDixonSolver (Method::Wiedemann());
    pass &= testFlatDixonSolver2 (Method::Wiedemann());
    pass &= testTallDixonSolver (Method::Wiedemann());
    pass &= testSingularDixonSolver (Method::SparseElimination());
    pass &= testZeroDixonSolver (Method::SparseElimination());
    pass &= testSingularDixonSolver (Method::BlasElimination());
    pass &= testZeroDixonSolver (Method::BlasElimination());
    pass &= testDixonSolverWithMPrhs ();
    pass &= testSparseRationalSolver ();
    pass &= testDixonRectangularSolver<> (Method::BlasElimination());
    pass &= testDixonRectangularSolver<> (Method::SparseElimination());
    pass &= testDixonRectangularSolver<> (Method::Wiedemann());
    pass &= testDixonRectangularSolver<DenseMatrix<ZRingInts>> (Method::BlasElimination());
    pass &= testDixonRectangularSolver<DenseMatrix<ZRingInts>> (Method::SparseElimination());
    pass &= testDixonRectangularSolver<DenseMatrix<ZRingInts>> (Method::Wiedemann());
    pass &= testSparse1x1Det(1<<26);
    pass &= testSparseDiagDet(46);
    pass &= testZeroDimensionalCharPoly ();
    pass &= testZeroDimensionalMinPoly ();
    pass &= testBigScalarCharPoly ();
    pass &= testLocalSmith ();
    /*
    pass &= testInconsistent<> (Method::BlasElimination());
    pass &= testInconsistent<> (Method::SparseElimination());
    pass &= testInconsistent<> (Method::Wiedemann());
    */

    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
