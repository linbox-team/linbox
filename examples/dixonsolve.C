/* Copyright (C) The LinBox group
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

/**\file examples/dixonsolve.C
 @example examples/dixonsolve.C
 @author Jean-Guillaume.Dumas@univ-grenoble-alpes.fr
 @author Romain.Lebreton@umontpellier.fr
 * \brief Dixon System Solving via Lifting using dense LU or sparse LU
 * \ingroup examples
 */
#include <iostream>

#include "linbox/algorithms/rational-solver.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/args-parser.h"
#include "linbox/util/error.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/randiter/random-prime.h"
#include <givaro/givrandom.h>

// DenseLU
#include "linbox/matrix/dense-matrix.h"

// SparseElim
#include "linbox/matrix/sparse-matrix.h"

using namespace LinBox;
typedef Givaro::ZRing<Givaro::Integer> Ints;
typedef DenseVector<Ints> ZVector;


template<typename _Matrix, typename _EliminationMethod>
int test(_Matrix A, std::string vector_file, bool inv, bool pp, bool sparse_elim) {
    // TODO : git rm dixondenseLU dixonsparseelim

    std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

    if (pp)
    {
            // Print Matrix

            // Matrix Market
            // std::cout << "A is " << A << std::endl;

            // Maple
        A.write(std::cout << "A:=", Tag::FileFormat::Maple) << ';' << std::endl;
    }

        // Vector File
    Ints ZZ;
    std::ifstream invect;

    ZVector B(ZZ, A.rowdim());

    bool createB = vector_file.empty();
    if (!createB) {
        invect.open (vector_file, std::ifstream::in);
        if (!invect) {
            createB = true;
        } else {
            for(ZVector::iterator it=B.begin(); it != B.end(); ++it)
                invect >> *it;
        }
    }

        // Vectors
    ZVector X(ZZ, A.coldim());

    if (createB) {
        ZVector U(ZZ, A.coldim());
        Givaro::GivRandom bgen( BaseTimer::seed() );
        if (inv) {
            std::cerr << "Creating a random {-1,1} vector " << std::endl;
            for(auto& it:B) it = (bgen.brand()?1:-1);
        } else {
            std::cerr << "Creating a random consistant {-1,1} vector " << std::endl;
            for(FFPACK::rns_double::integer& it:U) it = (bgen.brand()?1:-1);

            // B = A U
            A.apply(B,U);
        }
    }

    if(pp)
    {
            // Print RHS
        B.write(std::cout << "B:=", Tag::FileFormat::Maple) << ';' << std::endl;
    }

    std::cout << "B is " << B.size() << "x1" << std::endl;

    Timer chrono;

        // BlasElimination
    // TODO : à vérifier si cela marche avec Method::DenseElimination
    _EliminationMethod M;
    if (inv){
        M.singularity = Singularity::NonSingular;
    }


        //====================================================
        // BEGIN Replacement solve with fixed prime

    Ints::Element d;
    Method::Dixon m(M);
    typedef Givaro::Modular<double> Field;

    const size_t bitsize((size_t) FieldTraits<Field>::bestBitSize(A.rowdim()));
    Givaro::Integer randomPrime( *(PrimeIterator<>(bitsize)) );

    FixedPrimeIterator fixedprime( randomPrime );
    DixonSolver<Ints, Field, FixedPrimeIterator, _EliminationMethod> rsolve(A.field(), fixedprime);
    std::cout << "Using: " << *fixedprime << " as the fixed p-adic." << std::endl;

    chrono.start();
    if (!sparse_elim){
            // Dense Elimination
        if (inv)
        {
            std::cout << "Solving using Dense Elimination for non singular system" << std::endl;
            SolverReturnStatus ss = rsolve.solveNonsingular(X, d, A, B);
            if (ss != SS_OK) {
                std::cerr << "Error during solveNonsingular (possibly singular matrix or p-adic precision too small)" << std::endl;
                exit(-1);
            }
        }
        else
        {
            std::cout << "Solving using Dense Elimination for any system" << std::endl;
            SolverReturnStatus ss = rsolve.solve(X, d, A, B);
            if (ss == SS_FAILED){
                std::cerr << "Error during solve (all primes used were bad)" << std::endl;
                exit(-1);
            }
            if (ss == SS_INCONSISTENT){
                std::cerr << "Error: system appeared inconsistent" << std::endl;
                exit(-1);
            }
        }
    } else {
            // Sparse Elimination
        try
        {
            std::cout << "Solving using Sparse Elimination for any system" << std::endl;
            rsolve.solve(X, d, A, B);
        }
        catch(LinboxError& e)
        {
            std::cerr << e << '\n';
            exit(-1);
        }
    }
    chrono.stop();

    std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;

    {
            // Solution size

        std::cout<<"Reduced solution: \n";
        size_t maxbits=0;
        for (size_t i=0;i<A.coldim();++i){
            maxbits=(maxbits > X[i].bitsize() ? maxbits: X[i].bitsize());
        }
        std::cout<<" numerators of size   "<<maxbits<<" bits" << std::endl
                 <<" denominators hold over "<<d.bitsize()<<" bits\n";
    }


    {
			// Check Solution

        VectorDomain<Ints> VD(ZZ);
        MatrixDomain<Ints> MD(ZZ);
        ZVector LHS(ZZ, A.rowdim()), RHS(ZZ, B);
            // check that Ax = d.b
        MD.vectorMul(LHS, A, X);
        VD.mulin(RHS, d);
        if (VD.areEqual(LHS, RHS))
            std::cout << "Ax=d.b : Yes" << std::endl;
        else
            std::cout << "Ax=d.b : No" << std::endl;
    }

    {
            // Print Solution

        std::cout << "Solution is [";
        for(auto it:X) ZZ.write(std::cout, it) << " ";
        std::cout << "] / ";
        ZZ.write(std::cout, d)<< std::endl;
    }
    return 0;
}

int main (int argc, char **argv) {

    // TODO : seed ?
    std::string matrix_file = "";
    std::string vector_file = "";
    bool inv = false;
    bool pp = false;
    bool sparse_elim = false;


    Argument as[] = {
        { 'm', "-m FILE", "Set the input file for the matrix.",  TYPE_STR , &matrix_file },
        { 'v', "-v FILE", "Set the input file for the vector.",  TYPE_STR , &vector_file },
        { 'i', "-i"     , "whether the matrix is known to be invertible.",  TYPE_BOOL , &inv },
        { 'p', "-p"     , "whether you want to pretty print the matrix.",  TYPE_BOOL , &pp },
        { 's', "-s"     , "whether to use sparse elimination.",  TYPE_BOOL , &sparse_elim },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    // Matrix File
    if (matrix_file.empty()) {
        std::cerr << "You must specify an input file for the matrix with -m" << std::endl;
        exit(-1);
    }
    std::ifstream input (matrix_file);
    if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; exit(-1); }


        // Read Integral matrix from File
    Ints ZZ;
    MatrixStream< Ints > ms( ZZ, input );

    if (sparse_elim){
        SparseMatrix<Ints> A(ms);
        return test<SparseMatrix<Ints>,Method::SparseElimination>
            (A, vector_file, inv , pp, sparse_elim);
    } else {
        DenseMatrix<Ints> A(ms);
        return test<DenseMatrix<Ints>,Method::DenseElimination>
            (A, vector_file, inv , pp, sparse_elim);
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
