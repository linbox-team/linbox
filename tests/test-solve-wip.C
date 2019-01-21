// @fixme LICENCE

#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/solutions/solve-wip.h> // @fixme Just testing it compiles for now

using namespace LinBox;

using Field = Givaro::ZRing<Integer>;

template <class SolveMethod>
void run_2x2() {
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
    Field::Element xDen(1);

    // @note Alias to Method::CRAWIP<Method::Auto, Dispatch::Auto> m;
    solve(xNum, xDen, A, b, SolveMethod());

    std::cout << "---------------" << std::endl;
    A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
    b.write(std::cout << "b: ") << std::endl;
    std::cout << xNum << "/" << xDen << std::endl;

    if (xDen != 2 || xNum[0] != 2 || xNum[1] != 3) {
        std::cerr << "===> OUCH" << std::endl;
    }
}

int main(void)
{
    run_2x2<Method::Cra>();
    run_2x2<Method::Dixon>();

    return 0;
}
