/*
 * Copyright (C) LinBox
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
 */

#include <linbox/matrix/random-matrix.h>
#include <linbox/solutions/solve.h>

using namespace LinBox;

namespace {
    template <class Matrix>
    struct TypeToString {
        static constexpr const char* value = "Blackbox";
    };

    template <class Field>
    struct TypeToString<DenseMatrix<Field>> {
        static constexpr const char* value = "DenseMatrix";
    };

    template <class... MatrixArgs>
    struct TypeToString<SparseMatrix<MatrixArgs...>> {
        static constexpr const char* value = "SparseMatrix";
    };

    template <class Matrix>
    inline const char* type_to_string(const Matrix& A)
    {
        return TypeToString<Matrix>::value;
    }

    template <class SolveMethod, class ResultVector, class Matrix, class Vector>
    void print_error(const ResultVector& x, const Matrix& A, const Vector& b, std::string reason)
    {
        std::cerr << "/!\\ " << SolveMethod::name() << " on " << type_to_string(A) << " over ";
        A.field().write(std::cerr);
        std::cerr << " of size " << A.rowdim() << "x" << A.coldim() << " FAILS (" << reason << ")" << std::endl;
	A.write(std::cerr)<<std::endl;
	x.write(std::cerr)<<std::endl;
	b.write(std::cerr)<<std::endl;
    }

    template <class Matrix, class Domain>
    void generateMatrix(Domain& D, Matrix& A, int bitSize, int seed)
    {
        // @note RandomDenseMatrix can also fill a SparseMatrix, it has just a sparsity of 0.
        Givaro::Integer samplesize(1); samplesize <<= bitSize;
        typename Domain::RandIter randIter(D, seed, samplesize);
        LinBox::RandomDenseMatrix<typename Domain::RandIter, Domain> RDM(D, randIter);

        // @note Within our test m == n implies invertible,
        // but for the rectangular case, we pass whatever.
        if (A.rowdim() == A.coldim()) {
            RDM.randomFullRank(A);
        }
        else {
            RDM.random(A);
        }
    }

    template <class Matrix, class Vector, class Domain>
    void generateVector(Domain& D, const Matrix& A, Vector& b, int vectorBitSize, int seed)
    {
        Givaro::Integer samplesize(1); samplesize <<= vectorBitSize;
        typename Domain::RandIter vectorRandIter(D, seed, samplesize);

        // @note For the rectangular cases,
        // we set b so that the system is not inconsistent
        if (A.rowdim() == A.coldim()) {
            b.random(vectorRandIter);
        }
        else {
            Vector x(D, A.coldim());
            x.random(vectorRandIter);
            A.apply(b, x); // b <- Ax
        }
    }
}

template <class SolveMethod, class Matrix, class Vector, class ResultMatrix, class ResultVector>
bool check_result(ResultVector& x, Matrix& A, Vector& b, ResultMatrix& RA, ResultVector& Rb)
{
    ResultVector RAx(RA.field(), Rb.size());
    RA.apply(RAx, x);

    VectorDomain<typename ResultMatrix::Field> VD(RA.field());
    if (!VD.areEqual(RAx, Rb)) {
        print_error<SolveMethod>(x, A, b, "Ax != b");
        return false;
    }

    return true;
}

template <class SolveMethod, class Matrix, class Vector, class ResultDomain>
bool test_solve(const SolveMethod& method, Matrix& A, Vector& b, ResultDomain& RD, bool verbose)
{
    using ResultVector = DenseVector<ResultDomain>;

    if (verbose) {
        std::cout << "--- Testing " << SolveMethod::name() << " on " << type_to_string(A) << " over ";
        A.field().write(std::cout) << " of size " << A.rowdim() << "x" << A.coldim() << std::endl;
    }

    //
    // Copying data for final check
    //

    using ResultMatrix = typename Matrix::template rebind<ResultDomain>::other;

    // @note The following const_cast prevents this ambiguity warning:
    // warning: ISO C++ says that these are ambiguous
    // - linbox/matrix/sparsematrix/sparse-sequence-vector.h:660:3: note: candidate 1
    //      SparseMatrix (const SparseMatrix<_Tp1, _Rw1> &Mat, const Field& F)
    // - linbox/matrix/sparsematrix/sparse-sequence-vector.h:643:3: note: candidate 2
    //      SparseMatrix (const Field &F, VectStream &stream)
    ResultMatrix RA(A, const_cast<const ResultDomain&>(RD));
    ResultVector Rb(RD, b);

    //
    // Solving
    //

    ResultVector x(RD, A.coldim());

    bool ok = true;
    try {
        solve(x, A, b, method);
        ok = ok && check_result<SolveMethod>(x, A, b, RA, Rb);

        solveInPlace(x, A, b, method);
        ok = ok && check_result<SolveMethod>(x, A, b, RA, Rb);
    } catch (...) {
        print_error<SolveMethod>(x, A, b, "throws error");
        return false;
    }

    return ok;
}

template <class SolveMethod, class Domain, class ResultDomain>
bool test_dense_solve(const SolveMethod& method, Domain& D, ResultDomain& RD, int m, int n, int bitSize, int vectorBitSize, int seed,
                      bool verbose)
{
    using Vector = DenseVector<Domain>;

    DenseMatrix<Domain> A(D, m, n);
    Vector b(D, A.rowdim());

    generateMatrix(D, A, bitSize, seed);
    generateVector(D, A, b, vectorBitSize, seed + 1);

    return test_solve(method, A, b, RD, verbose);
}

template <class SolveMethod, class Domain, class ResultDomain>
bool test_sparse_solve(const SolveMethod& method, Domain& D, ResultDomain& RD, int m, int n, int bitSize, int vectorBitSize, int seed,
                       bool verbose)
{
    using Vector = DenseVector<Domain>;

    SparseMatrix<Domain> A(D, m, n);
    Vector b(D, A.rowdim());

    generateMatrix(D, A, bitSize, seed);
    generateVector(D, A, b, vectorBitSize, seed + 1);

    return test_solve(method, A, b, RD, verbose);
}

template <class SolveMethod, class Domain, class ResultDomain>
bool test_blackbox_solve(const SolveMethod& method, Domain& D, ResultDomain& RD, int m, int n, int bitSize, int vectorBitSize, int seed,
                         bool verbose)
{
    using Vector = DenseVector<Domain>;

    DenseMatrix<Domain> AData(D, n, m);
    generateMatrix(D, AData, bitSize, seed);

    // We're using a Transpose to use our blackbox interface
    Transpose<DenseMatrix<Domain>> A(AData);
    Vector b(D, A.rowdim());
    generateVector(D, A, b, vectorBitSize, seed + 1);

    return test_solve(method, A, b, RD, verbose);
}

int main(int argc, char** argv)
{
    Integer q = 131071;
    bool verbose = false;
    bool loop = false;
    int seed = -1;
    int bitSize = 10;
    int vectorBitSize = -1;
    int m = 32;
    int n = 24;
    std::string dispatchString = "Auto";

    static Argument args[] = {
        {'q', "-q", "Field characteristic.", TYPE_INTEGER, &q},
        {'v', "-v", "Enable verbose mode.", TYPE_BOOL, &verbose},
        {'l', "-l", "Infinite loop of tests.", TYPE_BOOL, &loop},
        {'s', "-s", "Seed for randomness.", TYPE_INT, &seed},
        {'b', "-b", "Bit size for rational solve tests.", TYPE_INT, &bitSize},
        {'B', "-B", "Vector bit size for rational solve tests (defaults to -b if not specified).", TYPE_INT, &vectorBitSize},
        {'m', "-m", "Row dimension of matrices.", TYPE_INT, &m},
        {'n', "-n", "Column dimension of matrices.", TYPE_INT, &n},
        {'d', "-d", "Dispatch mode (either Auto, Sequential, SMP or Distributed).", TYPE_STR, &dispatchString},
        END_OF_ARGUMENTS};

    parseArguments(argc, argv, args);

    // Setting up context

    Communicator communicator(0, nullptr);

    MethodBase method;
    method.pCommunicator = &communicator;
    method.dispatch = Dispatch::Auto;
    if (dispatchString == "Distributed")
        method.dispatch = Dispatch::Distributed;
    else if (dispatchString == "Sequential")
        method.dispatch = Dispatch::Sequential;
    else if (dispatchString == "SMP")
        method.dispatch = Dispatch::SMP;
    else if (dispatchString != "Auto") {
        std::cerr << "-d Dispatch mode should be either Auto, Sequential, SMP or Distributed" << std::endl;
        return EXIT_FAILURE;
    }

    if (vectorBitSize < 0) {
        vectorBitSize = bitSize;
    }

    if (seed < 0) {
        seed = time(nullptr);
    }

    if (verbose) {
        commentator().setReportStream(std::cout);
    }

    Givaro::ZRing<Givaro::Integer> ZZ;
    Givaro::QField<Givaro::Rational> QQ;
    Givaro::Modular<double> F(q);

    bool ok = true;
    do {
        // ----- Rational Auto
        ok = ok && test_dense_solve(Method::Auto(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
#if 0
        ok = ok && test_sparse_solve(Method::Auto(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        // @fixme Dixon<Wiedemann> does not compile
        // ok = ok && test_blackbox_solve(Method::Auto(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);

        ok = ok && test_dense_solve(Method::Auto(method), QQ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        ok = ok && test_sparse_solve(Method::Auto(method), QQ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        // ok = ok && test_blackbox_solve(Method::Auto(method), QQ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);

        // ----- Rational CRA
        // @fixme @bug When bitSize = 5 and vectorBitSize = 50, CRA fails
        ok = ok && test_dense_solve(Method::CRAAuto(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        ok = ok && test_sparse_solve(Method::CRAAuto(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        // ok = ok && test_blackbox_solve(Method::CRAAuto(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);

        ok = ok && test_dense_solve(Method::CRAAuto(method), QQ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        ok = ok && test_sparse_solve(Method::CRAAuto(method), QQ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        // ok = ok && test_blackbox_solve(Method::CRAAuto(method), QQ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);

        // ----- Rational Dixon
        ok = ok && test_dense_solve(Method::Dixon(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        ok = ok && test_sparse_solve(Method::Dixon(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);
        // @fixme Dixon<Wiedemann> does not compile
        // ok = ok && test_blackbox_solve(Method::Dixon(method), ZZ, QQ, m, n, bitSize, vectorBitSize, seed, verbose);

        // ----- Rational SymbolicNumeric
        // @note SymbolicNumeric methods are only implemented on DenseMatrix
        // @fixme Singular case fails
        // ok = ok && test_dense_solve(Method::SymbolicNumericOverlap(method), ZZ, QQ, m, n, bitSize, vectorBitSize,
        // seed, verbose);
        // @fixme Fails
        // ok = ok && test_sparse_solve(Method::SymbolicNumericNorm(method), ZZ, QQ, m, n, bitSize, vectorBitSize,
        // seed, verbose);

        // ----- Modular Auto
        ok = ok && test_dense_solve(Method::Auto(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_sparse_solve(Method::Auto(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_blackbox_solve(Method::Auto(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular Blackbox
        ok = ok && test_dense_solve(Method::Blackbox(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_sparse_solve(Method::Blackbox(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_blackbox_solve(Method::Blackbox(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular DenseElimination
        ok = ok && test_dense_solve(Method::DenseElimination(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_sparse_solve(Method::DenseElimination(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_blackbox_solve(Method::DenseElimination(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular SparseElimination
        ok = ok && test_dense_solve(Method::SparseElimination(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_sparse_solve(Method::SparseElimination(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_blackbox_solve(Method::SparseElimination(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular Wiedemann
        ok = ok && test_dense_solve(Method::Wiedemann(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_sparse_solve(Method::Wiedemann(method), F, F, m, n, 0, 0, seed, verbose);
        ok = ok && test_blackbox_solve(Method::Wiedemann(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular Lanczos
        // @fixme Dense is segfaulting
        // ok = ok && test_dense_solve(Method::Lanczos(method), F, F, m, n, 0, 0, seed, verbose);
        // @fixme Singular case fails
        // ok = ok && test_sparse_solve(Method::Lanczos(method), F, F, m, n, 0, 0, seed, verbose);
        // ok = ok && test_sparse_blackbox_solve(Method::Lanczos(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular BlockLanczos
        // @fixme Dense does not compile
        // ok = ok && test_dense_solve(Method::BlockLanczos(method), F, F, m, n, 0, 0, seed, verbose);
        // @fixme Sparse is segfaulting
        // ok = ok && test_sparse_solve(Method::BlockLanczos(method), F, F, m, n, 0, 0, seed, verbose);
        // ok = ok && test_blackbox_solve(Method::BlockLanczos(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular BlockWiedemann
        // @deprecated These do not compile anymore
        // ok = ok && test_dense_solve(Method::BlockWiedemann(method), F, F, m, n, 0, 0, seed, verbose);
        // ok = ok && test_sparse_solve(Method::BlockWiedemann(method), F, F, m, n, 0, 0, seed, verbose);
        // ok = ok && test_blackbox_solve(Method::BlockWiedemann(method), F, F, m, n, 0, 0, seed, verbose);

        // ----- Modular Coppersmith
        // @deprecated These do not compile anymore
        // ok = ok && test_dense_solve(Method::Coppersmith(method), F, F, m, n, 0, 0, seed, verbose);
        // ok = ok && test_sparse_solve(Method::Coppersmith(method), F, F, m, n, 0, 0, seed, verbose);
        // ok = ok && test_blackbox_solve(Method::Coppersmith(method), F, F, m, n, 0, 0, seed, verbose);
#endif
        if (!ok) {
            std::cerr << "Failed with seed: " << seed << std::endl;
        }

        seed += 1;
    } while (ok && loop);

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
