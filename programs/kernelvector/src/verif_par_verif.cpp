/*  Copyright (c) 2015 HPAC
 *  Written by  Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dums@imag.fr>
 */

#include <fflas-ffpack/utils/args-parser.h>
#include <fflas-ffpack/paladin/parallel.h>

#include "verif_read.hpp"
#include "verif_define.hpp"
#include "common_sdm.hpp"

//-----------------//
//----- Utils -----//

void readMAT(const Field_t& field, const std::string& matrixFile, DMatrix_t& matrix)
{
    std::ifstream matrixStream(matrixFile);
    Size_t m = 0u, n = 0u;
    std::string tmp;

    // Header
    matrixStream >> m >> n >> tmp ;
    matrix.resize(m * n);

    // Content
    for (auto& x : matrix)
        field.read(matrixStream, x);
}

void readList(const Field_t& field, const std::string& filename, std::vector<FieldElement_t>& r)
{
    std::ifstream ifs(filename);
    std::string tmp;
    Size_t rsize;

    // Header
    ifs >> tmp >> rsize;
    r.resize(rsize);

    // Content
    for (auto& ri : r)
        field.read(ifs, ri);
}

//----------------------------//
//----- Testing sequence -----//

//! for all j ; \sum_k r_K * W_{j*K + k}  == T * Y_j
TTimer testSequence(const std::string& seqFileBase, const std::string& chkFileBase,
                    const DMatrix_t& r, const DMatrix_t& T,
                    const Index_t K, const int64_t iteration,
                    const Field_t& field)
{
    TTimer timer;
    timer.clear();
    timer.start();

    bool failed = false;
    std::ifstream seqStream, chkStream;

    // Import headers
    Index_t d = 0u, dK = 0u, m = 0u, s1 = 0u, s2 = 0u;
    chkStream.open(chkFileBase + "_header." SDM_EXTENSION, std::ios::in | std::ios::binary);
    seqStream.open(seqFileBase + "_header." SDM_EXTENSION, std::ios::in | std::ios::binary);
    if (!chkStream.is_open()) throw std::logic_error("Checkpoint header file not found.");
    if (!seqStream.is_open()) throw std::logic_error("Sequence header file not found.");
    SDM_READ_HEADER(chkStream, dK, m, s2);
    SDM_READ_HEADER(seqStream, d, s1, s2);
    chkStream.close();
    seqStream.close();

    // Checking chk_i and seq_i
    if (iteration != -1u) {
        std::cerr << "[TSEQ] Testing pair " << iteration << "..." << std::endl;
        auto sNumber = formatNumber(iteration);
        chkStream.open(chkFileBase + sNumber + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        seqStream.open(seqFileBase + sNumber + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        if (!chkStream.is_open()) throw std::logic_error("Checkpoint file not found.");
        if (!seqStream.is_open()) throw std::logic_error("Sequence file not found.");
        dK = 1u;
    }
    // Choose the files containing all
    else {
        std::cerr << "[TSEQ] Testing " << dK << " elements in sequence..." << std::endl;
        chkStream.open(chkFileBase + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        seqStream.open(seqFileBase + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        
        // Discard headers
        Size_t d, m, n;
        SDM_READ_HEADER(chkStream, d, m, n);
        SDM_READ_HEADER(seqStream, d, m, n);
    }

    // Vectors
    DMatrix_t Yj(m * s2);
    DMatrix_t TYj(1u * s2);

    DMatrix_t WjK_k(s1 * s2);
    DMatrix_t SWr(1u * s2);

    timer.stop();
    
    std::cerr << "[TCHK] Sequence elements read OK in " << timer << "." << std::endl;

    timer.clear();
    timer.start();
// FFLAS::MMHelper<Field_t,FFLAS::MMHelperAlgo::Winograd,
//     typename FFLAS::ModeTraits<Field_t>::value,
//     FFLAS::ParSeqHelper::Parallel> 
//     WH (field, 0, SPLITTER(MAX_THREADS, FFLAS::RECURSIVE, FFLAS::TWO_D_ADAPT));

    // Real computation
    for (Index_t j = 0u; j < dK; ++j)
    {
        // Right part of equality
        // SDM_READ_MATRIX(field, chkStream, Yj);
        SDM_READ_MATRIX(chkStream, Yj); // Note: field is Montgomery, just a constant multiplication missing

// PAR_BLOCK{
//         FFLAS::fgemm(field, FFLAS::FflasTrans, FFLAS::FflasNoTrans, s1, s2, m,
//                      field.one, T.data(), s1, Yj.data(), s2, field.zero, TYj.data(), s2, WH);
// } 
        FFLAS::fgemm(field, FFLAS::FflasTrans, FFLAS::FflasNoTrans, 1, s2, m,
                     field.one, T.data(), 1, Yj.data(), s2, field.zero, TYj.data(), s2);

        // Left part of equality
        FFLAS::fzero(field, 1, s2, SWr.data(), s2);

        for (Index_t k = 0u; k < K; ++k)
        {
            // SDM_READ_MATRIX(field, seqStream, WjK_k);
            SDM_READ_MATRIX(seqStream, WjK_k); // Note: field is Montgomery, just a constant multiplication missing
            
            FFLAS::fgemm(field, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 1, s2, s1,
                         field.one, r.data()+s1*k, 1, WjK_k.data(), s2, 
                         field.one, SWr.data(), s2);
       }

        // Check
        if (SWr != TYj)
        {
            std::cerr << "[TSEQ] /!\\ Check " << j << " FAILED." << std::endl;
            failed = true;
        }
        else {
            std::cerr << "[TSEQ] Check " << j << " is OK." << std::endl;
        }
    }

    timer.stop();

    if (!failed) std::cerr << "[TSEQ] Verif finished OK in " << timer << "." << std::endl;
    else         std::cerr << "[TSEQ] /!\\ Verif finished WITH ERRORS in " << timer << "." << std::endl;


    return timer;
}

//-------------------------------//
//----- Testing checkpoints -----//

// for all j ; P * Y_j == Z * Y_{j-1}
TTimer testCheckpoints(const std::string& chkFileBase, const Field_t& field, const DMatrix_t& P, const DMatrix_t& Z,
                       const int64_t iteration)
{
    TTimer timer;
    timer.clear();
    timer.start();

    bool failed = false;
    std::ifstream chkStream;

    // Import header
    Index_t dK = 0u, m = 0u, s2 = 0u;
    chkStream.open(chkFileBase + "_header." SDM_EXTENSION, std::ios::in | std::ios::binary);
    if (!chkStream.is_open()) throw std::logic_error("Checkpoint header file not found.");
    SDM_READ_HEADER(chkStream, dK, m, s2);
    chkStream.close();

    // Matrices
    DMatrix_t Yi(m * s2);
    DMatrix_t PYi(1u * s2);
    DMatrix_t ZYim1(1u * s2);

    // Checking chk_i vs chk_{i+1}
    if (iteration != -1u) {
        std::cerr << "[TCHK] Testing checkpoint " << iteration + 1u << "..." << std::endl;
        chkStream.open(chkFileBase + formatNumber(iteration) + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        if (!chkStream.is_open()) throw std::logic_error("Checkpoint file not found.");
        
        // Read initial
        // SDM_READ_MATRIX(field, chkStream, Yi);
        SDM_READ_MATRIX(chkStream, Yi); // Note: field is Montgomery, just a constant multiplication missing
        chkStream.close();
        
        // Open next checkpoint file
        chkStream.open(chkFileBase + formatNumber(iteration + 1u) + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        if (!chkStream.is_open()) throw std::logic_error("Checkpoint file not found.");
        dK = 2u;
    }
    // Choose the file containing all
    else {
        std::cerr << "[TCHK] Testing " << dK << " checkpoints..." << std::endl;
        chkStream.open(chkFileBase + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        if (!chkStream.is_open()) throw std::logic_error("Checkpoint file not found.");
        
        // Discard header
        Size_t d, m, n;
        SDM_READ_HEADER(chkStream, d, m, n);
        
        // Read initial
        // SDM_READ_MATRIX(field, chkStream, Yi);
        SDM_READ_MATRIX(chkStream, Yi); // Note: field is Montgomery, just a constant multiplication missing
        // Read Y_j and Y_{j-1} will have the same constant factor, so egality check is correct.
        // with SDM_READ_MATRIX(field, chkStream, Yi);, we read Y_j and Y_{j-1} correctly
        // with SDM_READ_MATRIX(chkStream, Yi);       , we read Y_j.B and Y_{j-1}.B, it is faster, but it does not matter ..!
    }
    timer.stop();
    
    std::cerr << "[TCHK] Checkpoint read OK in " << timer << "." << std::endl;

    timer.clear();
    timer.start();

    // Real computation

    for (Index_t i = 0u; i < dK - 1u; ++i) {
        // Z * Y_{i-1}
        FFLAS::fgemm(field, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 1u, s2, m,
                     field.one, Z.data(), m, Yi.data(), s2, field.zero, ZYim1.data(), s2);

        // SDM_READ_MATRIX(field, chkStream, Yi);
        SDM_READ_MATRIX(chkStream, Yi); // Note: field is Montgomery, just a constant multiplication missing

        // P * Yi
        FFLAS::fgemm(field, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 1u, s2, m,
                     field.one, P.data(), m, Yi.data(), s2, field.zero, PYi.data(), s2);

        if (PYi != ZYim1) {
            std::cerr << "[TCHK] /!\\ Check " << i << " FAILED." << std::endl;
            failed = true;
        }
        else {
            std::cerr << "[TCHK] Check " << i << " is OK." << std::endl;
        }
    }

    timer.stop();

    if (!failed) std::cerr << "[TCHK] Verif finished OK in " << timer << "." << std::endl;
    else std::cerr << "[TCHK] /!\\ Verif finished WITH ERRORS in " << timer << "." << std::endl;

    return timer;
}

//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::string folder = "output/";
    int64_t nb_checkpoints = -1u;
    size_t checks = 3;
    int seed = -1;

    static Argument as[] =
    {
        { 'd', "-d FOLDER", "The directory where ." SDM_EXTENSION " are stored.", TYPE_STR, &folder },
        { 'i', "-i NBCHECKPOINTS", "Test for nbcheckpoints iterations.", TYPE_INT, &nb_checkpoints },
        { 'c', "-c CHECKS", "[0,1,2,3]: check checkpoints 0/1 and/or sequence 00/10.", TYPE_INT, &checks },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);
    if (seed < 0) seed = time(nullptr);
    FFLAS::writeCommandString(std::cerr << "Command line: ", as, argv[0]) << std::endl;
    PAR_BLOCK { std::cerr << "[TCHK] NUM_THREADS=" << NUM_THREADS << std::endl; }

    // Get values
    Modulus_t signedq;
    Index_t m, n, s1, s2, d, K;
    std::string matrixFile;
    loadInfo(folder + "/info.txt", matrixFile, signedq, m, n, s1, s2, d, K);
    const Field_t::Residu_t q(signedq);

    // Construct field
    std::ios_base::sync_with_stdio(false);
    Field_t field(q);

    // File name bases
    std::string seqFileBase = folder + "/seq";
    std::string chkFileBase = folder + "/chk";

    TTimer chrono; 

    if (checks & 0x1) {
    //----- Checkpoints test -----//
        chrono.clear(); chrono.start();

            // Get P
        DMatrix_t P;
        readMAT(field, folder + "/verif_p.mat", P);
        if (P.size() != m) throw std::logic_error("P size is not what expected.");

            // Get Z
        DMatrix_t Z;
        readMAT(field, folder + "/verif_z.mat", Z);
        if (Z.size() != m) throw std::logic_error("Z size is not what expected.");
   
        chrono.stop();
        std::cerr << "[TCHK] P and Z read OK in " << chrono << "." << std::endl;

        chrono.clear(); chrono.start();
            // Check that for all j ; P * Y_j == Z * Y_{j-1}
        PARFOR1D(iteration,nb_checkpoints,FFLAS::ParSeqHelper::Parallel<>(),
                 testCheckpoints(chkFileBase, field, P, Z, iteration);
                 );
        chrono.stop();
        std::cerr << "[TCHK] done in " << chrono << "." << std::endl;
    }
    
    if (checks & 0x2) {
            //----- Sequence test -----//
        chrono.clear(); chrono.start();

            // Get r
        DMatrix_t r;
        readMAT(field, folder + "/verif_r.mat", r);
        if (r.size() != s1 * K) throw std::logic_error("r size is not what expected.");

            // Get T
        DMatrix_t T;
        readMAT(field, folder + "/verif_t.mat", T);
        if (T.size() != s1 * m) throw std::logic_error("T size is not what expected.");

        chrono.stop();
        std::cerr << "[TSEQ] R and T read OK in " << chrono << "." << std::endl;

        chrono.clear(); chrono.start();
            // Check that for all j ; \sum_k W_{j*K + k} * r_k == T * Y_j
        PARFOR1D(iterseq,nb_checkpoints,FFLAS::ParSeqHelper::Parallel<>(),
                 testSequence(seqFileBase, chkFileBase, r, T, K, iterseq, field);
                 );
        chrono.stop();
        std::cerr << "[TSEQ] done in " << chrono << "." << std::endl;
    }
    
    return EXIT_SUCCESS;
}

