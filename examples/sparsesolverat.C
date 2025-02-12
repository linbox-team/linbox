/*
 * examples/sparsesolverat.C
 *
 * Copyright (C) The LinBox Group
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
 */

/**\file examples/sparsesolverat.C
 @example examples/sparsesolverat.C
 @author Jean-Guillaume.Dumas@univ-grenoble-alpes.fr
 * \brief Direct sparse solver over the rationals 
 * \ingroup examples
 */
#include <iostream>

#include "givaro/modular.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/methods.h"

using namespace LinBox;


template<typename DVector, typename EDom>
int rhs(DVector& B, const EDom& DD, bool createB, std::ifstream& invect) {
    if (createB) {
        std::cerr << "Creating a random {-1,1} vector " << std::endl;
        for(auto it=B.begin(); it != B.end(); ++it)
            if (drand48() <0.5)
                *it = -1;
            else
                *it = 1;
    } else {
        for(auto&& it:B) invect >> it;
        invect.close();
    }

    std::clog << "B is [";
    for(auto it:B) DD.write(std::clog, it) << ' ';
    std::clog << ']' << std::endl;

    return 0;
}

int main (int argc, char **argv)
{

        // set 2 to see Q L U P factorization; 
        // see fill-in with 2 and __LINBOX_ALL__ or __LINBOX_FILLIN__ defined
	commentator().setMaxDetailLevel (1); 
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::clog);
//     commentator().setDefaultReportFile("/dev/stdout"); // to see activities

	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: solve <matrix-file-in-supported-format> [<dense-vector-file>] [0/1 <integer solve>]" << std::endl;
		return 0;
	}
	std::ifstream input (argv[1]);
	if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }
	std::ifstream invect;

	bool createB = false;
	if (argc == 2) {
		createB = true;
	}

    bool integralsolve=false;

	if (argc >= 3) {
		invect.open (argv[2], std::ifstream::in);
		if (!invect) {
			createB = true;
            integralsolve = atoi(argv[2]);
		}
		else {
			createB = false;
		}
	}

	if (argc >= 4)
        integralsolve = atoi(argv[3]);

    typedef Givaro::QField<Givaro::Rational> Rats;
    Rats QQ;
    typedef DenseVector<Rats> RVector;

    if (integralsolve) {
        std::clog << "Integral solving" << std::endl;

        typedef Givaro::ZRing<Givaro::Integer> Ints;
        Ints ZZ;
        typedef DenseVector<Ints> ZVector;

        MatrixStream<Rats> ms( QQ, input );
        size_t nrow, ncol; ms.getDimensions(nrow,ncol);

        SparseMatrix<Ints> A(ZZ,nrow,ncol);
        ZVector X(ZZ, A.coldim()),B(ZZ, A.rowdim());

            // Read rational matrix and rhs, then compute denominator LCM
        {
            SparseMatrix<Rats> RA ( ms );
            std::clog << "A is " << RA.rowdim() << " by " << RA.coldim() << std::endl;
            Givaro::Integer ABlcm(1);

            for(auto iterow = RA.rowBegin() ; iterow != RA.rowEnd(); ++iterow) {
                for(auto iter = iterow->begin(); iter != iterow->end(); ++iter) {
                    lcm(ABlcm, ABlcm, iter->second.deno());
                }
            }

            std::clog << "A denominator lcm: " << ABlcm << std::endl;

            RVector RB(QQ, RA.rowdim());
            rhs(RB, QQ, createB, invect);
            for(auto iter = RB.begin(); iter != RB.end(); ++iter) {
                lcm(ABlcm, ABlcm, iter->deno());
            }

            std::clog << "A & B denominator lcm: " << ABlcm << std::endl;


                // A x = b is equivalent to (l.A) x = (l.b)
            auto iterow = RA.rowBegin();
            auto iterit = A.rowBegin();

            for(; iterow != RA.rowEnd(); ++iterow, ++iterit) {
                for(auto iter = iterow->begin(); iter != iterow->end(); ++iter) {
                    iterit->emplace_back( iter->first, (iter->second.nume() * ABlcm) / iter->second.deno() );
                }
            }
            auto iter = RB.begin();
            auto itez = B.begin();
            for( ; iter != RB.end(); ++iter, ++itez) {
                *itez = (iter->nume() * ABlcm) / iter->deno();
            }

        }

        Givaro::ZRing<Integer>::Element d;

        Timer chrono;

        std::clog << "Integral Sparse Elimination" << std::endl;
        chrono.start();
        solveInPlace (X, d, A, B, Method::SparseElimination());
        chrono.stop();

        std::cout << "(SparseElimination) Solution is [";
        for(auto it:X) ZZ.write(std::cout, it) << ' ';
        std::cout << "] / ";
        ZZ.write(std::cout, d)<< std::endl;

        std::clog << "CPU time (seconds): " << chrono.usertime() << std::endl;

    } else {

        MatrixStream<Rats> ms( QQ, input );

        SparseMatrix<Rats> A ( ms );
		std::clog << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

		RVector X(QQ, A.coldim()),B(QQ, A.rowdim());

		// Sparse Elimination
        rhs(B, QQ, createB, invect);

        Timer chrono;

        std::clog << "Direct Rational Sparse Elimination" << std::endl;
        chrono.start();
            // @fixme Can't pass a Randiter anymore, we need an API through Method to set seed or such
            // typename Rats::RandIter generator(QQ,0,BaseTimer::seed() );
            // solveInPlace (X, A, B, Method::SparseElimination(), generator);
        solveInPlace (X, A, B, Method::SparseElimination());
        chrono.stop();

        std::clog << "(SparseElimination) Solution is [";
        for(auto it:X) QQ.write(std::cout, it) << ' ';
        std::clog << ']' << std::endl;

        std::clog << "CPU time (seconds): " << chrono.usertime() << std::endl;
	}

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
