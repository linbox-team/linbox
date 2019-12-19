/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include "common_define.hpp"
#include "common_read.hpp"
#include "common_spmm.hpp"
#include "common_dv.hpp"

#include <fflas-ffpack/utils/args-parser.h>
#include <fflas-ffpack/paladin/parallel.h>

#include <cstdlib>

//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::string folder = "output/";

    static Argument as[] =
    {
        { 'd', "-d FOLDER",  "The directory where kernel.dv is stored.", TYPE_STR, &folder },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);
    
    std::cerr << "[FCHK] Checking that the computed vector is correct..." << std::endl;

    // Get values
    Modulus_t p;
    Size_t m, n, s1, s2, d, K;
    std::string matrixFile;
    loadInfo(folder + "/info.txt", matrixFile, p, m, n, s1, s2, d, K);

    // Init matrix from stream
    SMatrix_t A;
    getMatrix(p, matrixFile, A);

    std::cerr << "[FCHK] Matrix read." << std::endl;
    //----- Final checks -----//
    
    // Check that A * X == 0 and X != 0
    DVector_t X;
    std::ifstream kernelStream(folder + "/kernel.dv");
    readDV(kernelStream, X);
    
    DVector_t Y(X.size());

    if (isNullVector(p,X)) {
        std::cerr << "[ERROR] The computed vector is the null one. Sorry. Try again." << std::endl;
        return EXIT_FAILURE;
    }
    
    PAR_BLOCK { spmv(p, Y, A, X); }

    if (!isNullVector(p,Y)) {

        std::cerr << "[FCHK] First attempt failed : the computed vector is not from the kernel..." << std::endl;
	std::cerr<<"[FCHK] Second chance: trying to apply A on the vector"<<std::endl;
	DVector_t Z(X.size());
	PAR_BLOCK { spmv(p, Z, A, Y); }
	if (!isNullVector(p,Z)) {
	  std::ofstream verifStream(folder + "/outvect.dv");
	  writeDV(verifStream,Y);
	  //        writeDV(std::cout,Y);
	  std::cerr<<"[ERROR] still not a vector from the Kernel. Sorry. Try again."<<std::endl;
	  return EXIT_FAILURE;
	} else
	  {
	    std::string cmd= "cp "; cmd+= folder + "/kernel.dv "; cmd+= folder + "/kernel1.dv "; 
	    system(cmd.c_str());
	    std::ofstream kernelvecStream(folder + "/kernel.dv");
	    writeDV(kernelvecStream,Y);
	    std::cerr<<"[SUCCESS] A * kernelvec is a vector from the kernel of A"<<std::endl;
	    std::cerr<<"[INFO] Saving the result in kernel.dv"<<std::endl;
	  }
    }

    std::cerr << "/!\\ ALL IS OK! a non-null kernel vector has been found !" << std::endl;

    return EXIT_SUCCESS;
}

