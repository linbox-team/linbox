// @fixme LICENCE

#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/solve-wip.h>

using namespace LinBox;

template <class Field, class Method>
void run_integer(Communicator& communicator, size_t dimension) {
    Field F;

    BlasVector<Field> b(F, dimension);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    DenseMatrix<Field> A(F, dimension, dimension);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    Method method;
    method.pCommunicator = &communicator;
    // @fixme Dixon fails with dimension = 3
    // when solution type is kept default (that is to say Determinist)
    method.solutionType = SolutionType::Diophantine;

    // (num, den) interface
    {
        std::cout << "--------------- " << Method::name() << " (" << dimension << ") [(num, den) API]" << std::endl;

        BlasVector<Field> xNum(F, dimension);
        typename Field::Element xDen(1);
        solve(xNum, xDen, A, b, method);

        if (xDen != 2 || xNum[0] != 2 || xNum[1] != 3) {
            A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
            std::cout << "b: " << b << std::endl;
            std::cout << "x: " << xNum << "/" << xDen << std::endl;
            std::cerr << "===> OUCH" << std::endl;
        }
    }

    // Rational vector interface
    {
        Givaro::ZRing<Givaro::Rational> F;
        BlasVector<Givaro::ZRing<Givaro::Rational>> x(F, dimension);

        std::cout << "--------------- " << Method::name() << " (" << dimension << ") [Rational vector API]" << std::endl;

        solve(x, A, b, method);

        if (x[0].nume() != 1 || x[0].deno() != 1 || x[1].nume() != 3 || x[1].deno() != 2) {
            A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
            std::cout << "b: " << b << std::endl;
            std::cout << "x: " << x << std::endl;
            std::cerr << "===> OUCH" << std::endl;
        }
    }
}

template <class Matrix, class Method>
void run_modular() {
    using Field = typename Matrix::Field;

    Field F(101);

    BlasVector<Field> b(F, 2);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    Matrix A(F, 2, 2);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    Method method;
    method.blockingFactor = 1;

    std::cout << "--------------- " << Method::name() << std::endl;

    BlasVector<Field> x(F, 2);
    solve(x, A, b, method);

    if (x[0] != 1 || x[1] != 52) {
        A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "x: " << x << std::endl;
        std::cerr << "===> OUCH" << std::endl;
    }
}

int main(void)
{
    // @fixme -v for verbose mode (commentator on)

    Communicator communicator(0, nullptr);

    commentator().setReportStream(std::cout);

    // @fixme Uncomment these working tests (currently commented to speed up compilation times)
    run_integer<Givaro::ZRing<Integer>, MethodWIP::Auto>(communicator, 2);
    run_integer<Givaro::ZRing<Integer>, MethodWIP::CraAuto>(communicator, 2);
    run_integer<Givaro::ZRing<Integer>, MethodWIP::CraAuto>(communicator, 3);
    run_integer<Givaro::ZRing<Integer>, MethodWIP::DixonAuto>(communicator, 2);
    run_integer<Givaro::ZRing<Integer>, MethodWIP::DixonAuto>(communicator, 3);

    run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::Auto>();
    run_modular<SparseMatrix<Givaro::Modular<double>>, MethodWIP::Auto>();
    run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::DenseElimination>();
    run_modular<SparseMatrix<Givaro::Modular<double>>, MethodWIP::DenseElimination>();
    run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::SparseElimination>();
    run_modular<SparseMatrix<Givaro::Modular<double>>, MethodWIP::SparseElimination>();
    // run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::Wiedemann>(); // @fixme Can't compile
    run_modular<SparseMatrix<Givaro::Modular<double>>, MethodWIP::Wiedemann>();
    // run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::Lanczos>(); // @fixme Segmentation fault
    // run_modular<SparseMatrix<Givaro::Modular<double>>, MethodWIP::Lanczos>(); // @fixme Segmentation fault
    // run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::BlockLanczos>(); // @fixme Can't compile
    // run_modular<SparseMatrix<Givaro::Modular<double>>, MethodWIP::BlockLanczos>(); // @fixme Segmentation fault

    // @deprecated These do not compile anymore
    // run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::BlockWiedemann>();
    // run_modular<DenseMatrix<Givaro::Modular<double>>, MethodWIP::Coppersmith>();

    return 0;
}
