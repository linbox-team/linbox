// @fixme LICENCE

#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/solve-wip.h> // @fixme Just testing it compiles for now

using namespace LinBox;

template <class Field, class Method>
void run_rational_2x2() {
    Field F;

    BlasVector<Field> b(F, 3);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    BlasMatrix<Field> A(F, 3, 3);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    BlasVector<Field> xNum(F, 3);
    typename Field::Element xDen(1);

    solve(xNum, xDen, A, b, Method());

    std::cout << "--------------- " << Method::name() << std::endl;
    A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "x: " << xNum << "/" << xDen << std::endl;

    if (xDen != 2 || xNum[0] != 2 || xNum[1] != 3) {
        std::cerr << "===> OUCH" << std::endl;
    }
}

template <class Field, class Method>
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

    solve(x, A, b, Method());

    std::cout << "--------------- " << Method::name() << std::endl;
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
    run_rational_2x2<Givaro::ZRing<Integer>, MethodWIP::Cra>();
    run_rational_2x2<Givaro::ZRing<Integer>, MethodWIP::Dixon>();

    run_2x2<Givaro::Modular<float>, MethodWIP::DenseElimination>();
    run_2x2<Givaro::Modular<float>, MethodWIP::SparseElimination>();
    run_2x2<Givaro::Modular<float>, MethodWIP::Wiedemann>();
    run_2x2<Givaro::Modular<float>, MethodWIP::BlockWiedemann>();
    run_2x2<Givaro::Modular<float>, MethodWIP::Coppersmith>();

    return 0;
}
