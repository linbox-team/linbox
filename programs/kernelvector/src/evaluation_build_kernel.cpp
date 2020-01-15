/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include <fflas-ffpack/utils/args-parser.h>

#include <givaro/modular-general.h>
#include <linbox/blackbox/diagonal.h>
#include <linbox/blackbox/permutation.h>

#include "common_define.hpp"
#include "common_read.hpp"
#include "common_spmm.hpp"
#include "common_dv.hpp"

#include "common_kernel.hpp"

#include <recint/rmint.h>

#ifdef __USE_128bits
using FieldElement_t = RecInt::rmint<7u, RecInt::MG_INACTIVE>;
#else
using FieldElement_t = RecInt::rmint<6u, RecInt::MG_INACTIVE>;
#endif


//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::cerr << "command-line:";
    for(int i = 0; i < argc; ++i) 
       std::cerr << ' ' << argv[i];
    std::cerr << std::endl;

    std::string folder = "output/";
    std::string matrixFile ="A.sms";
    Givaro::Integer ip(101);

    static Argument as[] =
    {
        { 'd', "-d FOLDER",  "The directory where results (kernel.dv, *.mat?.sms) are stored.", TYPE_STR, &folder },
        { 'f', "-f FILE",  "The matrix file name.", TYPE_STR, &matrixFile },
        { 'p', "-p MOD", "The (Integer) modulus.", TYPE_INTEGER,   &ip },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);
    Modulus_t p((Modulus_t)ip);
    //std::string filename=removeExtension(matrixFile);
    std::string filename=folder+std::string("/")+removeExtension(basename(matrixFile));
    
    std::cerr << "[FCHK] Checking that the computed vector for B is correct..." << std::endl;

    // Get values
    // Init matrix B from stream
    SMatrix_t M;
    getMatrix(p, filename+".matB.sms", M);

    std::cerr << "[FCHK] Modulo " << p << ", B Matrix read." << std::endl;
    //----- Final checks -----//
    
    // Check that A * X == 0 and X != 0
    DVector_t X;
    std::ifstream kernelStream(folder + "/kernel.dv");
    readDV(kernelStream, X);
    kernelStream.close();
    
    DVector_t Y(X.size());
    
    if (isNullVector(p,X)) {
        std::cerr << "/!\\ The computed vector is the null one... Sorry. Try again..." << std::endl;
        return EXIT_FAILURE;
    }
    
    PAR_BLOCK { spmv(p, Y, M, X); }

    
    if (! isNullVector(p,Y)) {
        writeDV(std::cout,Y);
        std::cerr << "/!\\ Nah, the computed vector is not from kernel of B... Sorry. Game over..." << std::endl;
        return EXIT_FAILURE;
    }
    
    std::cerr << "[FCHK] ALL IS OK! It's a non-null kernel vector for B!" << std::endl;

    
    // Get values
    // Init matrix C from stream
    getMatrix(p, filename+".matC.sms", M, false);
    std::cerr << "[REKE] C Matrix read." << std::endl;
    Y.resize(M.nRows);

    PAR_BLOCK { spmv(p, Y, M, X); }
    
    std::cerr << "[REKE] Cv computed." << std::endl;
    typedef  Givaro::ZRing<Modulus_t> ZDom;
    ZDom Zi;

    LinBox::Diagonal< ZDom > D(Zi);
    std::ifstream matDiag(filename+".matD.sms");
    try { D.read(matDiag); } catch(const LinBox::MatrixStreamError& e) { 
        if ((D.rowdim()==0) || (D.coldim()==0))
            std::cerr << " But it is fine, diagonal matrix is the zero matrix, continuing ... " << std::endl;
        else
            throw e;
    }
    std::cerr << "[REKE] D matrix read." << std::endl;
    
    LinBox::Diagonal< ZDom >::Vector_t& vD = D.getData();
    
    FieldElement_t::init_module(p);
    for(size_t ii=0; ii<D.rowdim();++ii) {
        FieldElement_t Dii(vD[ii]); 
// Zi.write(std::cerr << "1/", vD[ii]) << " mod " << p << '=';
        inv(Dii);
// Zi.write(std::cerr, Dii) << " mod " << p  << ';' << std::endl;
        neg(Dii);
        FieldElement_t Yii(Y[ii]);
// Zi.write(std::cerr, Y[ii]) << "*(";
// Zi.write(std::cerr, Dii) << ") mod " << p << '=';
	Yii *= Dii;
        Zi.convert(Y[ii], Yii);
// Zi.write(std::cerr, Y[ii]) << " mod " << p  << ';' << std::endl;
    }    
    
//         F.invin(vD[ii]);
//         F.negin(vD[ii]);
//         F.mulin(Y[ii],vD[ii]);
    
    std::cerr << "[REKE] -D^{-1}Cv computed." << std::endl;

    // Get values
    // Init matrix A from stream
    getMatrix(p, matrixFile, M, false);

    std::cerr << "[REKE] Initial Matrix read." << std::endl;


    Y.insert(Y.end(),X.begin(),X.end());
    Y.insert(Y.end(),M.nCols-X.size()-D.rowdim(),Zi.zero);

 //    writeDV(std::cout << "Y:" ,Y);

    std::ifstream matPerm(filename+".matP.sms");
    LinBox::Permutation<ZDom> P(Zi, M.nCols);
    P.read(matPerm);
    std::cerr << "[REKE] Permutation read: " << P.rowdim() << 'x' << P.coldim() <<std::endl;
    

    DVector_t Z(Y.size());
    P.applyTranspose(Z,Y);

    std::ofstream output(filename+".kernel.dv");
    writeDV(output,Z);
    std::cerr << "[REKE] Kernel vector produced." << std::endl;

    PAR_BLOCK { spmv(p, Y, M, Z); }

     if (! isNullVector(p,Y)) {
         writeDV(std::cout << "ERROR:\n",Y);
        std::cerr << "/!\\ Nah, the computed vector is not from kernel of B... Sorry. Game over..." << std::endl;
        return EXIT_FAILURE;
    }
    
     //std::cerr << "[FINAL]/!\\ ALL IS OK! It's a non-null kernel vector !!!" << std::endl;

   

    return EXIT_SUCCESS;
}

