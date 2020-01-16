/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include <givaro/givrandom.h>
#include <fflas-ffpack/utils/args-parser.h>

#include "verif_read.hpp"
#include "verif_define.hpp"
#include "common_sdm.hpp"

//---------------------------//
//----- Testing minpoly -----//

// \sum_i UT * W_i * G_i * V == 0
TTimer testMinpoly(const std::string& seqFile, const std::string& polyFile,
                   const Field_t& field, const DMatrix_t& U, const DMatrix_t& V, int seed)
{
    // Load streams
    std::ifstream seqStream(seqFile);
    std::ifstream polyStream(polyFile);
    if (!seqStream.is_open()) throw std::logic_error("Sequence file not found.");
    if (!polyStream.is_open()) throw std::logic_error("Minpoly file not found.");

    // Import headers
    Index_t d, s1, s2, ng;
    SDM_READ_HEADER(seqStream, d, s1, s2);
    SDM_READ_HEADER(polyStream, ng, s2, s2);

    // Matrices
    FieldElement_t r;
    DMatrix_t Wi(s1 * s2);
    DMatrix_t UWi(1u * s2);
    DMatrix_t Gi(s2 * s2);
    DMatrix_t GiV(s2 * 1u);
    
    field.assign(r, field.zero);

    std::cerr << "[TPLY] " << s1 << 'x' << s2 << " sequence of length " << d <<  " ... "  << std::endl;

    Givaro::GivRandom gen(seed);
    const size_t shift=gen() % (d-ng);

    // Real computation
    TTimer timer;
    timer.clear();
    timer.start();
    for(size_t j=0;j<shift;++j) 
        SDM_READ_MATRIX(seqStream, Wi);
    timer.stop();

    std::cerr << "[TPLY] Sequence shifted by " << shift << " in " << timer << std::endl;

    std::cerr << "[TPLY] Starting " << s2 << 'x' << s2 << " minpoly verif of degree " << ng << " ..." << std::endl;

    timer.clear();
    timer.start();
    
    for (Index_t i = 0u; i < ng; ++i) {
        // U * W_i
        // SDM_READ_MATRIX(field, seqStream, Wi);
        SDM_READ_MATRIX(seqStream, Wi); // Note: field is Montgomery, just a constant multiplication missing
        FFLAS::fgemm(field, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 1u, s2, s1,
                     field.one, U.data(), s1, Wi.data(), s2, field.zero, UWi.data(), s2);
    
        // G_i * V
        // SDM_READ_MATRIX(field, polyStream, Gi);
        SDM_READ_MATRIX(field, polyStream, Gi); // Note: field is Montgomery, just a constant multiplication missing
        FFLAS::fgemm(field, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s2, 1u, s2,
                     field.one, Gi.data(), s2, V.data(), 1u, field.zero, GiV.data(), 1u);
                     
        // r += U * W_i * G_i * V
        for (Index_t j = 0u; j < s2; ++j)
            field.axpyin(r, UWi.at(j), GiV.at(j));
    }

    timer.stop();

    if (field.isZero(r)) std::cerr << "[TPLY] Verif finished OK in " << timer << "." << std::endl;
    else std::cerr << "[TPLY] /!\\ Verif finished WITH ERRORS in " << timer << "." << std::endl;

    return timer;
}

//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::string folder = "output/";
    int seed = -1;

    static Argument as[] =
    {
        { 'd', "-d FOLDER", "The directory where ." SDM_EXTENSION " are stored.", TYPE_STR, &folder },
        { 'g', "-g SEED", "Seed for randomness.", TYPE_INT, &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);
    if (seed < 0) seed = time(nullptr);
    FFLAS::writeCommandString(std::cerr << "Command line: ", as, argv[0]) << std::endl;

    // Get values
    Modulus_t signedq;
    Index_t m, n, s1, s2, d, K;
    std::string matrixFile;
    loadInfo(folder + "/info.txt", matrixFile, signedq, m, n, s1, s2, d, K);
//     std::cerr << "matrixFile: " << matrixFile << std::endl;
//     std::cerr << "modulus: " << signedq << std::endl;
//     std::cerr << "nRows: " << m << std::endl;
//     std::cerr << "nCols: " << n << std::endl;
//     std::cerr << "s1: " << s1 << std::endl;
//     std::cerr << "s2: " << s2 << std::endl;
//     std::cerr << "d: " << d << std::endl;
//     std::cerr << "K: " << K << std::endl;

    const Field_t::Residu_t q(signedq);

    // Construct field
    std::ios_base::sync_with_stdio(false);
    Field_t field(q);
    Field_t::RandIter g(field, seed);

    // File name bases
    std::string seqFileBase = folder + "/seq";
    std::string chkFileBase = folder + "/chk";

    //----- Minpoly test -----//

    // Generate U
    DMatrix_t U(1u * s1);
    for (auto& u : U)
        g.random(u);

    // Generate V
    DMatrix_t V(s2 * 1u);
    for (auto& v : V)
        g.random(v);

    // Check that \sum_i UT * W_i * G_i * V == 0
    testMinpoly(folder + "/seq." SDM_EXTENSION, folder + "/poly." SDM_EXTENSION,
                field, U, V, seed);

    return EXIT_SUCCESS;
}

