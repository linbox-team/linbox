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
#include "linbox-config.h"
#include "givaro/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/vector/blas-vector.h"
#include "solutions/solve.h"
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

#if LINBOX_HAVE_SAGE
#include "linbox/linbox-sage.h"
struct c_vector_modint_linbox {
	// copy of the declaration in vector_modn_sparse.pxi
	int *entries;
	int p;
	size_t *positions;
	size_t degree;
	size_t num_nonzero;
};
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

int main (int argc, char **argv)
{
    bool pass = true;

    if (!testSolveSparse  ()) pass = false;
    if (!testSolveSparseSage  ()) pass = false;

    return pass ? 0 : -1;
}
