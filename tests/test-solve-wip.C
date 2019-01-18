// @fixme LICENCE

#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/solutions/solve-wip.h> // @fixme Just testing it compiles for now

using namespace LinBox;

int main(void)
{
    using Field = Givaro::ZRing<Integer>;
    Field F;

    BlasVector<Field> b(F, 2);
    BlasMatrix<Field> A(F, 2, 2);

    BlasVector<Field> xNum(F, 2);
    Field::Element xDen(1);

    // @note Alias to Method::CRAWIP<Method::Auto, Dispatch::Auto> m;
    solve(xNum, xDen, A, b, Method::Cra());

    return 0;
}
