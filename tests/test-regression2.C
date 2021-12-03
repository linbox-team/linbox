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
#include "linbox/blackbox/zo-gf2.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/solutions/solve/solve-wiedemann.h"

bool writing=false;


using namespace LinBox;

template<class Ring, class Matrix>
bool testWiedemannSingular() {
// A = [
// 1, 0, 1, 0;
// 0, 1, 1, 1;
// 0, 0, 0, 0;
// 1, 1, 0, 1;
// 1, 0, 1, 0
// ]
// b = [1; 1; 0; 0; 1]
// r = [0; 0; 1; 0]
    Ring gf2(2);

    Matrix A(gf2, 5,4);
    A.setEntry(0,0,gf2.one);A.setEntry(0,0,gf2.one);A.setEntry(0,2,gf2.one);
    A.setEntry(1,1,gf2.one);A.setEntry(1,2,gf2.one);A.setEntry(1,3,gf2.one);
    A.setEntry(3,0,gf2.one);A.setEntry(3,1,gf2.one);A.setEntry(3,3,gf2.one);
    A.setEntry(4,0,gf2.one);A.setEntry(4,2,gf2.one);

    DenseVector<Ring> X(gf2,4), B(gf2, 5);
    gf2.assign(B[0],gf2.one);
    gf2.assign(B[1],gf2.one);
    gf2.assign(B[4],gf2.one);

	solve (X, A, B, RingCategories::ModularTag(), Method::Wiedemann());

	DenseVector<Ring> r(gf2, A.rowdim());
    A.apply(r,X);
	VectorDomain<Ring> VD(gf2);
	if (VD.areEqual (r,B))
        return true;
    else
        return false;
}



int main (int argc, char **argv)
{
    bool pass = true;

	// text is written to clog/cerr/cout iff a command line argument is present.
	if (argc > 1) writing = true;

    if (writing) {
        commentator().setReportStream(std::cout);
    }

    using GMd = Givaro::Modular<double>;
    using MMatrix = SparseMatrix<GMd>;


    pass &= testWiedemannSingular<GMd,MMatrix>();

        // Still failing: see https://github.com/linbox-team/linbox/issues/233
    // using GFt = GF2;
    // using ZMatrix=ZeroOne<GFt>;
    // pass &= testWiedemannSingular<GFt,ZMatrix>();

    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
