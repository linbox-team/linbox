// @fixme LICENCE

#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/solve-wip.h> // @fixme Just testing it compiles for now

using namespace LinBox;

template <class Field, class Method>
void run_rational(Communicator& communicator, size_t dimension) {
    Field F;

    BlasVector<Field> b(F, dimension);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    BlasMatrix<Field> A(F, dimension, dimension);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    BlasVector<Field> xNum(F, dimension);
    typename Field::Element xDen(1);

    Method method;
    method.pCommunicator = &communicator;
    // @fixme Dixon fails with dimension = 3
    // solution type is kept default (that is to say Determinist)
    method.solutionType = SolutionType::Diophantine;
    solve(xNum, xDen, A, b, method);

    std::cout << "--------------- " << Method::name() << " (" << dimension << ")" << std::endl;

    if (xDen != 2 || xNum[0] != 2 || xNum[1] != 3) {
        A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "x: " << xNum << "/" << xDen << std::endl;
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

    if (x[0] != 1 || x[1] != 52) {
        A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "x: " << x << std::endl;
        std::cerr << "===> OUCH" << std::endl;
    }
}

int main(void)
{
    Communicator communicator(0, nullptr);

    // @fixme Test high-level? Auto, Elimination, Blackbox

    // @fixme Test Cra with different underlying Method
    run_rational<Givaro::ZRing<Integer>, MethodWIP::Cra>(communicator, 2);
    run_rational<Givaro::ZRing<Integer>, MethodWIP::Cra>(communicator, 3);
    run_rational<Givaro::ZRing<Integer>, MethodWIP::Dixon>(communicator, 2);
    run_rational<Givaro::ZRing<Integer>, MethodWIP::Dixon>(communicator, 3);

    run_2x2<Givaro::Modular<float>, MethodWIP::DenseElimination>();
    run_2x2<Givaro::Modular<float>, MethodWIP::SparseElimination>();
    // run_2x2<Givaro::Modular<float>, MethodWIP::Wiedemann>();
    // run_2x2<Givaro::Modular<float>, MethodWIP::BlockWiedemann>();
    // run_2x2<Givaro::Modular<float>, MethodWIP::Coppersmith>();

    return 0;
}
