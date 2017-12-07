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
using namespace LinBox;

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

        // solving via Sparse Elimination
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.one) || !ZZ.areEqual(A[1],ZZ.zero) || !ZZ.areEqual(D,ZZ.one)) {
        std::cerr<<"Fail solving a flat system over QQ with a SparseMatrix"<<std::endl;
        return false;
    }

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

        // solving via Sparse Elimination
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.one) || !ZZ.areEqual(D,ZZ.one)) {
        std::cerr<<"A = "<<A<<" D = "<<D<<std::endl;
        std::cerr<<"Fail solving a flat system over QQ with a SparseMatrix"<<std::endl;
        return false;
    }

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

        // solving via Sparse Elimination
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.one) || !ZZ.areEqual(D,ZZ.one)) {
        std::cerr<<"Fail solving a tall system over QQ with a SparseMatrix"<<std::endl;
        return false;
    }

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
        std::cerr<<"A = "<<A<<" D = "<<D<<std::endl;
        std::cerr<<"Fail solving a singular system over QQ with a SparseMatrix"<<std::endl;
        return false;
    }
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

        // solving via Sparse Elimination
    solve (A, D, M, B, m);

    if (!ZZ.areEqual(A[0],ZZ.zero) || !ZZ.areEqual(D,ZZ.one)) {
        std::cerr<<"A = "<<A<<" D = "<<D<<std::endl;
        std::cerr<<"Fail solving a zero over QQ with a SparseMatrix"<<std::endl;
        return false;
    }
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

int main (int argc, char **argv)
{
    bool pass = true;

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
    pass &= testZeroDimensionalCharPoly ();
    pass &= testZeroDimensionalMinPoly ();
    pass &= testBigScalarCharPoly ();

    return pass ? 0 : -1;
}
