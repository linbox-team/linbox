/* Copyright (C) LinBox
 *
 * Author: Zhendong Wan
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file  tests/test-rational-solver.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */

#include "linbox/linbox-config.h"
#include "linbox/ring/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/vector/stream.h"
#include "linbox/util/commentator.h"

#include "test-common.h"

#include <iostream>

using namespace LinBox;

/// Testing Nonsingular Random Diagonal solve.
template <class Ring, class Field, class Vector>
bool testRandomSolve (const Ring& R,
              const Field& f,
              LinBox::VectorStream<Vector>& stream1,
              LinBox::VectorStream<Vector>& stream2)
{
    commentator().start("Testing Nonsingular Random Diagonal solve ",
                        "testNonsingularRandomDiagonalSolve",
                        stream1.size());

    bool ret = true;

    VectorDomain<Ring> VD(R);
    Vector d(R), b(R), x(R), y(R);

    VectorWrapper::ensureDim (d, stream1.n ());
    VectorWrapper::ensureDim (b, stream1.n ());
    VectorWrapper::ensureDim (x, stream1.n ());
    VectorWrapper::ensureDim (y, stream1.n ());

    int n = (int)d.size();

    while (stream1 && stream2) {
        commentator().startIteration ((unsigned)stream1.j ());
        bool iter_passed = true;

        bool zeroEntry;
        do {
            stream1.next (d);
            zeroEntry = false;
            for (size_t i=0; i<stream1.n(); i++)
            zeroEntry |= R.isZero(d[(size_t)i]);
        } while (zeroEntry);

        stream2.next (b);

        std::ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        report << "Diagonal entries: ";
        VD.write (report, d);
        report << endl;

        report << "Right-hand side:  ";
        VD.write (report, b);
        report << endl;

        //Diagonal<Ring> D(R, d);
        BlasMatrix<Ring> D(R, n, n);

        for(int i = 0; i < n; ++i) R.init (D[(size_t)i][(size_t)i],  d[(size_t)i]);

        typedef DixonSolver<Ring, Field, PrimeIterator<IteratorCategories::HeuristicTag> > RSolver;
        RSolver rsolver;

        BlasVector<Ring> num(R,(size_t)n);
        typename Ring::Element den;

        auto solveResult = rsolver.solve(num, den, D, b, 30); // often 5 primes are not enough

        if (solveResult == SS_OK) {
          D. apply (y, num);
          VD. mulin(b, den);

          if (!VD.areEqual (y, b)) {
            ret = iter_passed = false;
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
              << "ERROR: Computed solution is incorrect" << endl;
          }
        }
        else {
            ret = iter_passed = false;
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
              << "ERROR: Did not return OK solving status" << endl;
        }

        commentator().stop ("done");
                commentator().progress ();
    }

    stream1.reset ();
    stream2.reset ();
    commentator().stop (MSG_STATUS (ret), (const char *) 0, "testNonsingularRandomDiagonalSolve");

    return ret;
}

int main(int argc, char** argv)
{
    bool pass = true;
    static size_t n = 10;
    static int iterations = 1;
    static Argument args[] = {
        { 'n', "-n N", "Set order of test matrices to N.", TYPE_INT, &n},
        { 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
        END_OF_ARGUMENTS };

    parseArguments (argc, argv, args);

    using Field = Givaro::Modular<int32_t>;
    using Ring  = Givaro::IntegerDom;

    Ring R; Ring::RandIter gen(R);
    Field F(101);

    RandomDenseStream<Ring> s1 (R, gen, n, (unsigned int)iterations), s2 (R, gen, n, (unsigned int)iterations);
    if (!testRandomSolve(R, F, s1, s2)) pass = false;

    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
