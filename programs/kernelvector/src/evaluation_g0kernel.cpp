/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */


#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/ffpack/ffpack.h>
#include <fflas-ffpack/utils/args-parser.h>

#include "common_define.hpp"
#include "common_sdm.hpp"
#include "common_read.hpp"
#include "common_spmm.hpp"
#include "common_dv.hpp"
using Index_t = uint64_t;       // Sufficiently large to store an array index.

#include <givaro/montgomery-ruint.h>

#ifdef __USE_128bits
using FieldElement_t = RecInt::ruint128;
#else
using FieldElement_t = RecInt::ruint64;
#endif
using Field_t = Givaro::Montgomery<FieldElement_t>;

//------------------//
//----- Kernel -----//

// Should be called directly after header is read from stream
// e should be full of field.zero
inline void getG0KernelVector(const Field_t& field, std::istream& polyStream, std::vector<FieldElement_t>& e)
{
    const Size_t s2 = e.size();
    
    std::vector<FieldElement_t> G0(s2 * s2);
    std::vector<FieldElement_t> nullMatrix(s2 * s2, field.zero);
    std::vector<FieldElement_t> nullVector(s2, field.zero);

    SDM_READ_MATRIX(field, polyStream, G0);

    // If G0 is full of zeros, fixed e to e1
    if (G0 == nullMatrix)
    {
        e[0] = field.one;
        return;
    }

    // Find a nullspace vector
    int loop = 0;
    while (e == nullVector && ++loop < 5000)
    {
        FFPACK::RandomNullSpaceVector(field, FFLAS::FFLAS_SIDE::FflasRight, s2, s2,
                                      G0.data(), s2, e.data(), 1);
    }

    // No nullspace vector
    if (e == nullVector) {
        std::cerr << "[G0KL] /!\\ Cannot find not null kernel vector of G0." << std::endl;
        exit(EXIT_FAILURE);
    }
}

//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::string folder = "output/";

    static Argument as[] =
    {
        { 'd', "-d FOLDER", "The directory where poly." SDM_EXTENSION " is stored.", TYPE_STR, &folder },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);

    TTimer seqTimer;

    std::cerr << "[G0KL] Importing initial values... ";
    seqTimer.clear(); seqTimer.start();

    // Get values
    Modulus_t signedp;
    Size_t m, n, s1, s2, d, K;
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
        std::cerr << "[G0KL] /!\\ Polynomial file not found or readable." << std::endl;
        exit(EXIT_FAILURE);
    }
    seqTimer.stop();

    std::cerr << "done in " << seqTimer << std::endl;

    std::cerr << "[G0KL] Computing g0 kernel... ";
    seqTimer.clear(); seqTimer.start();

    // Compute kernel
    std::vector<FieldElement_t> e(s2, field.zero);
    getG0KernelVector(field, polyStream, e);

    seqTimer.stop();
    std::cerr << "done in " << seqTimer << std::endl;

    // Export kernel
    std::ofstream kernelStream(folder + "/g0kernel.dv");
    writeDV(field, kernelStream, e);
        
    std::cerr << "[G0KL] Exported G0 kernel vector to " << folder << "/g0kernel.dv" << std::endl;

    return EXIT_SUCCESS;
}

