#include <iostream>
#include <fstream>
#include <linbox/config-blas.h>
#include <linbox/linbox-config.h>
#include <linbox/field/modular-double.h>
#include <linbox/fflas/fflas.h>
#include <linbox/util/timer.h>

#include "test-common.h"

using namespace LinBox;
int main (int argc, char ** argv) {
    size_t n=300, nmax=1000, prec=256;
   
    static Argument args[] = {
        { 'n', "-n n", "Operate over the \"field\" GF(Q) [1] for integer modulus.", TYPE_INT, &n },
        { 'm', "-m m", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INT, &nmax },
        { 'p', "-p p", "Operate over the \"field\" GF(Q) [1] for uint16 modulus.", TYPE_INT, &prec },
        { '\0' }
    };

    parseArguments (argc, argv, args);

    commentator.start("Optimization suite", "Optim");
    std::ostream& report = commentator.report();

    Modular<double> F(17);
    LinBox::Timer chrono;

    double *A, *C;
    A = new double[nmax*nmax];
    C = new double[nmax*nmax];
    for (size_t i=0; i<nmax*nmax;++i){
        A[i]=random() % 17;
    }

    do {
        chrono.start();
        FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                     n, n, n, 1., A, n, A, n, 0., C, n);
        chrono.stop();
        report << std::endl
                  << "fgemm " << FFLAS::WinoSteps(n) << "Wino: " << n << "x" << n << ": "
                  << chrono.usertime() << " s, "
                  << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
                  << std::endl;

            n+=prec;
    } while ((n < nmax));

    delete[] A;
    delete[] C;

    return 0;
}
