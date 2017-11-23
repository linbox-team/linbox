/* linbox-sage.C
 * Copyright (C) 2007 Martin Albrecht
 *               2008 Clement Pernet
 *
 * Written by Martin Albrecht
 *            Clement Pernet <clement.pernet@gmail.com>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#include <iostream>

#include <cstdlib>
#include <vector>
#include <list>
#include "linbox/linbox-config.h"
#include "linbox/util/commentator.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/vector/sparse.h"

#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/solve.h"
#include "linbox/solutions/methods.h"
#include "linbox/integer.h"

#include <gmp.h>

#include "linbox-sage.h"

using namespace LinBox;


/*************************************************************************
  sparse modulo Z/nZ
 *************************************************************************/

struct c_vector_modint_linbox {
	// copy of the declaration in vector_modn_sparse.pxi
	int *entries;
	int p;
	size_t *positions;
	size_t degree;
	size_t num_nonzero;
};

typedef unsigned int mod_int;
typedef Givaro::Modular<unsigned int> GFp;
typedef GFp::Element  GFpElement;
typedef std::vector <std::pair <size_t, GFpElement> > SparseSeqVectorGFp;
typedef SparseMatrix<GFp, VectorTraits<SparseSeqVectorGFp>::SparseFormat> SparseMatrixGFp;

SparseMatrixGFp* linbox_new_modn_sparse_matrix(mod_int modulus, size_t numrows, size_t numcols, void *rows)
{
	GFp *F=new GFp(modulus);
	SparseMatrixGFp* M=new SparseMatrixGFp(*F, numrows, numcols);

	struct c_vector_modint_linbox *A = static_cast<struct c_vector_modint_linbox *>(rows);

	for(size_t i = 0; i < numrows; ++i) {
		for(size_t j = 0; j < A[i].num_nonzero; ++j) {
                        M->setEntry(i, A[i].positions[j], A[i].entries[j]);
		}
	}
	return M;
}

void linbox_delete_modn_sparse_matrix (SparseMatrixGFp * A) {
        const GFp * F = &A->field();
        delete A;
        delete F;
}

static std::vector<GFpElement> linbox_new_modn_sparse_vector(mod_int modulus, size_t len, void *_vec)
{
	std::vector<GFpElement> A(len);

	if (_vec==NULL) {
		return A;
	}

	struct c_vector_modint_linbox *vec = static_cast<struct c_vector_modint_linbox*>(_vec);
	for(size_t i = 0; i < vec->num_nonzero; ++i) {
		A[vec->positions[i]] = vec->entries[i];
	}
	return A;
}

unsigned long linbox_modn_sparse_matrix_rank(mod_int modulus,
					     size_t numrows, size_t numcols,
					     void *rows, int gauss)
{
	unsigned long M_rank;
	GFpElement M_det;

	SparseMatrixGFp *M = linbox_new_modn_sparse_matrix(modulus, numrows, numcols, rows) ;
	const GFp &F = M->field();
	GaussDomain<GFp> dom(F);

	if(!gauss) {
		dom.InPlaceLinearPivoting(M_rank, M_det, *M, numrows, numcols);
	}
	else {
		dom.NoReordering(M_rank, M_det, *M, numrows, numcols);
	}

	//*pivots = (int*)calloc(sizeof(int), dom.pivots.size());

	//   int j=0;
	//   for(std::vector<int>::const_iterator i= dom.pivots.begin(); i!= dom.pivots.end(); ++i, ++j){
	//     (*pivots)[j] = *i;
	//   }
        linbox_delete_modn_sparse_matrix(M);
	return M_rank;
}

std::vector<mod_int> linbox_modn_sparse_matrix_solve(mod_int p, size_t numrows, size_t numcols,
						     void *_a, void *b, int method)
{
	// solve ax = b, for x, a matrix, b vector, x vector

	SparseMatrixGFp *A =linbox_new_modn_sparse_matrix(p, numrows, numcols, _a);

	const GFp & F = A->field();

        DenseVector<GFp> X(F, numrows);
        DenseVector<GFp> B(F, linbox_new_modn_sparse_vector(p, numcols, b) );

	switch(method) {
	case 1:
		solve(X, *A, B, Method::BlasElimination());
		break;

	case 2:
		solve(X, *A, B, Method::Blackbox());
		break;

	case 3:
		solve(X, *A, B, Method::Wiedemann());
		break;

	default:
		solve(X, *A, B);
	}

        linbox_delete_modn_sparse_matrix(A);

	return X.refRep();
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

