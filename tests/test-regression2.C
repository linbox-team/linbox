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
#include "linbox/solutions/solve/solve-dixon.h"
#include "linbox/solutions/solve/solve-wiedemann.h"

bool writing=false;


using namespace LinBox;

template<class Ring, class Matrix>
bool testWiedemannSingular() {
    std::clog << "Test Wiedemann Singular ... ";
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
	if (VD.areEqual (r,B)) {
        std::clog << "PASSED.\n";
        return true;
    } else {
        std::clog << "ERROR.\n";
        return false;
    }
}

bool testDixonDetOne(size_t count) {
    std::clog << "Test Dixon Det 1 ... ";
    using Ring = Givaro::ZRing<Integer>;
    using IArray = std::vector<Integer>;
    using IVector = DenseVector<Ring>;
    using IMatrix = DenseMatrix<Ring>;

    Ring Z;

    IMatrix A(Z,2,2, IArray{26,5,5,1}.begin() );
    A.write(std::clog, Tag::FileFormat::Maple) << std::endl;

    bool pass(true);
    IVector x(Z,2), v(Z,2); Integer d;

    for(size_t i=0; i<count; ++i) {
        IVector b(Z, Ring::RandIter(Z), 2);
//         b.write(std::clog, Tag::FileFormat::Maple) << std::endl;

        solve(x,d,A,b,RingCategories::IntegerTag(), Method::Dixon());

//         x.write(std::clog, Tag::FileFormat::Maple) << std::endl;
        A.apply(v,x);

        pass &= (d == Z.one) && VectorDomain<Ring>(Z).areEqual(v,b) ;
    }

    IVector b(Z, IArray{-959558580,-451007454}.begin(), 2);
    b.write(std::clog, Tag::FileFormat::Maple) << std::endl;

    solve(x,d,A,b,RingCategories::IntegerTag(), Method::Dixon());

    x.write(std::clog, Tag::FileFormat::Maple) << std::endl;
    A.apply(v,x);

    pass &= (d == Z.one) && VectorDomain<Ring>(Z).areEqual(v,b);

	if (pass) {
        std::clog << "PASSED.\n";
    } else {
        std::clog << "ERROR.\n";
    }
    return pass;
}



int main (int argc, char **argv)
{
    bool pass(true);

	// text is written to clog/cerr/cout iff a command line argument is present.
	if (argc > 1) writing = true;

    if (writing) {
        commentator().setReportStream(std::cout);
        commentator().setMaxDepth(-1);
        commentator().setMaxDetailLevel(-1);
    }

    using GMd = Givaro::Modular<double>;
    using MMatrix = SparseMatrix<GMd>;


    pass &= testWiedemannSingular<GMd,MMatrix>();

        // Still failing: see https://github.com/linbox-team/linbox/issues/233
    // using GFt = GF2;
    // using ZMatrix=ZeroOne<GFt>;
    // pass &= testWiedemannSingular<GFt,ZMatrix>();

    pass &= testDixonDetOne(10);

    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
