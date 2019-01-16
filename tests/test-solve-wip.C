// @fixme LICENCE

#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/solutions/solve-wip.h> // @fixme Just testing it compiles for now

using namespace LinBox;

int main(void)
{
    using Field = Givaro::ZRing<Integer>;
    Field F;

    MdrVector<BlasVector<Field>> x(F, 2);
    BlasVector<Field> b(F, 2);
    BlasMatrix<Field> A(F, 2, 2);

    Method::CRAWIP<Method::Hybrid, Dispatch::Auto> m;

    solve(x, A, b, m);

    return 0;
}
