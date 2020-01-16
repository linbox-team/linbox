/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */


#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/ffpack/ffpack.h>
#include <fflas-ffpack/utils/args-parser.h>

#include "common_sdm.hpp"
#include "verif_define.hpp"

#include <givaro/montgomery-ruint.h>

inline bool checkG0KernelVector(const Field_t& field, std::istream& polyStream, std::istream& g0Stream, size_t s2)
{
    std::cerr << "[G0CK] Verifying 1x" << s2 << " g0 kernel... " << std::endl;

    TTimer timer; timer.clear(); timer.start();
    DMatrix_t e(1 * s2);
    SDM_READ_VECTOR(g0Stream, e, 1u);
    

    DMatrix_t G0(s2 * s2);
    DMatrix_t res(1 * s2);
    DMatrix_t nullVector(1 * s2, field.zero);
    DMatrix_t nullMatrix(s2 * s2, field.zero);

    // No nullspace vector
    if (e == nullVector) {
        std::cerr << "[G0CK] /!\\ ERROR kernel vector of G0 is null." << std::endl;
        exit(EXIT_FAILURE);
    }

    SDM_READ_MATRIX(field, polyStream, G0);

    // If G0 is full of zeros, fixed e to e1
    if (G0 == nullMatrix)
    {
        timer.stop();
        std::cerr << "[G0CK] OK, G0 is zero, in " << timer << std::endl;
        return EXIT_SUCCESS;
    }

    FFLAS::fgemm(field, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s2, 1, s2, 
                 field.one, res.data(), 1, G0.data(), s2, 
                 field.zero, e.data(), 1);

    // No nullspace vector
    if (res != nullVector) {
        std::cerr << "[G0CK] /!\\ ERROR vector is not in the kernel of G0." << std::endl;
        exit(EXIT_FAILURE);
    }

    timer.stop();
    std::cerr << "[G0CK] G0 kernel vector OK in " << timer << std::endl;
    
    return EXIT_SUCCESS;
}

//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::string folder = "output/";

    static Argument as[] =
    {
        { 'd', "-d FOLDER", "The directory where poly / g0kernel." SDM_EXTENSION " are stored.", TYPE_STR, &folder },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);

    TTimer seqTimer;

    std::cerr << "[G0CK] Importing initial values... ";
    seqTimer.clear(); seqTimer.start();

    // Get values
    Modulus_t signedp;
    Size_t m, n, s1, s2, gs2, d, K;
    std::string matrixFile;
    loadInfo(folder + "/info.txt", matrixFile, signedp, m, n, s1, s2, d, K);
    const Field_t::Residu_t p(signedp);

    // Construct field
    std::ios_base::sync_with_stdio(false);
    ::srand(time(nullptr));
    RecInt::srand(time(nullptr));
    Field_t field(p);

    //----- Real computing -----//
    
    std::ifstream polyStream(folder + "/poly." SDM_EXTENSION, std::ios::in | std::ios::binary);
    SDM_READ_HEADER(polyStream, d, s2, s2);

    if (!polyStream.is_open()) {
        std::cerr << std::endl << "[G0CK] /!\\ Polynomial file not found or readable." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::ifstream g0Stream(folder + "/g0kernel.dv", std::ios::in | std::ios::binary);
    SDM_READ_HEADER(g0Stream, gs2);

    if (!g0Stream.is_open()) {
        std::cerr << std::endl<< "[G0CK] /!\\ Kernel file not found or readable." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (s2 != gs2) {
        std::cerr << std::endl<< "[G0CK] /!\\ Polynomial file and g0kernel diemnsions incompatible." << std::endl;
        exit(EXIT_FAILURE);
    }
    seqTimer.stop();

    std::cerr << "done in " << seqTimer << std::endl;

    return checkG0KernelVector(field, polyStream, g0Stream, gs2);
}

