// @fixme LICENCE

#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/solutions/solve-wip.h> // @fixme Just testing it compiles for now

using namespace LinBox;

template <class Field, class SolveMethod>
void run_rational_2x2() {
    Field F;

    BlasVector<Field> b(F, 2);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    BlasMatrix<Field> A(F, 2, 2);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    BlasVector<Field> xNum(F, 2);
    typename Field::Element xDen(1);

    solve(xNum, xDen, A, b, SolveMethod());

    std::cout << "---------------" << std::endl;
    A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "x: " << xNum << "/" << xDen << std::endl;

    if (xDen != 2 || xNum[0] != 2 || xNum[1] != 3) {
        std::cerr << "===> OUCH" << std::endl;
    }
}

template <class Field, class SolveMethod>
void run_2x2() {
    Field F(101);

    BlasVector<Field> b(F, 2);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    BlasMatrix<Field> A(F, 2, 2);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    BlasVector<Field> x(F, 2);

    solve(x, A, b, SolveMethod());

    std::cout << "---------------" << std::endl;
    A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "x: " << x << std::endl;

    if (x[0] != 1 || x[1] != 52) {
        std::cerr << "===> OUCH" << std::endl;
    }
}

int main(void)
{
    // @fixme Test high-level? Auto, Elimination, Blackbox

    // @fixme Test Cra with different underlying Method
    // @note Alias to Method::CRAWIP<Method::Auto, Dispatch::Auto> m;
    run_rational_2x2<Givaro::ZRing<Integer>, Method::Cra>();
    run_rational_2x2<Givaro::ZRing<Integer>, Method::Dixon>();

    run_2x2<Givaro::Modular<float>, Method::DenseElimination>();
    run_2x2<Givaro::Modular<float>, Method::SparseElimination>();
    run_2x2<Givaro::Modular<float>, Method::Wiedemann>();
    run_2x2<Givaro::Modular<float>, Method::BlockWiedemann>();
    run_2x2<Givaro::Modular<float>, Method::Coppersmith>();

    return 0;
}
