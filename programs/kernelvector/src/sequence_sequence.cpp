/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include "common_define.hpp"
#include "common_sdm.hpp"
#include "common_read.hpp"
#include "common_spmm.hpp"

#include <fflas-ffpack/utils/args-parser.h>
#include <fflas-ffpack/paladin/parallel.h>

//------------------//
//----- Savers -----//

void saveLeftProj(const std::string& filename, const LeftProj_t& U)
{
    std::ofstream ofs(filename);

    ofs << U.size() << std::endl;
    for (const auto u : U)
        ofs << u << std::endl;
}

void readLeftProj(LeftProj_t& U, size_t b1, const std::string& filename)
{
    std::ifstream ifs(filename);
    getProj(ifs, U, b1);
}

//! Prints Y_i to a stream
// We export it row-major
inline
void printCheckpoint(std::ostream& YOut, const DMatrix_t& Y)
{
    SDM_WRITE_MATRIX(YOut, Y.data);
}

//! Prints U.Y_i to a stream
// We export it row-major
inline
void printProjection(std::ostream& WOut, const LeftProj_t& U, const DMatrix_t& Y)
{
    DVector_t W;

    for (const auto r : U)
        for (Index_t j = 0u; j < Y.nCols; ++j)
            W.emplace_back(Y.data.at(r * Y.nCols + j));

    SDM_WRITE_MATRIX(WOut, W);
}

//--------------------//
//----- Sequence -----//

//! Y = A * Y, using tY as temporary
inline
void nextElement(const Modulus_t p, const SMatrix_t& A, DMatrix_t& Y, DMatrix_t& tY)
{
    spmm(p, tY, A, Y);  // tY = A * Y
    std::swap(Y, tY);   // Just 3 pointers
}

/*! Name conventions:
 *  A   The matrix of size m x m
 *  U   The left projection of size s2 x m (in fact, just a vector of size s2 with row indexes)
 *  V   The right matrix of size m x s1
 *  W   The sequence    (U * A^i * V)
 *  Y   The checkpoints (A^k * Y)
 */
void generateSequence(const std::string& seqFileBase, const std::string& chkFileBase, Modulus_t p, const size_t restart, const Size_t d, const Size_t K,
                      const SMatrix_t& A, const LeftProj_t& U, const DMatrix_t& V,
                      const int64_t stopAfter, const bool breaking)
{
    std::cerr << "[SEQC] Starting sequence generation..." << std::endl;
    std::cerr << "[SEQC] SIMD optimization uses SIMD_STEP=" << __SIMD_STEP << std::endl;

    std::ofstream WOut, YOut;
    std::string sNumber = (breaking)? formatNumber(restart) : "";
    WOut.open(seqFileBase + sNumber + "."+ SDM_EXTENSION, std::ios::out | std::ios::binary);
    YOut.open(chkFileBase + sNumber + "."+ SDM_EXTENSION, std::ios::out | std::ios::binary);
    
    if (!breaking) {
        SDM_WRITE_HEADER(WOut, d, U.size(), V.nCols);
        SDM_WRITE_HEADER(YOut, d/K, V.nRows, V.nCols);
    }

    // Y_0 = V
    DMatrix_t Y(V);
    DMatrix_t tY(A.nRows, V.nCols);

        // Prevent Invalid access if A.nRows is not divisible by __SIMD_STEP
    tY.data.reserve( (A.nRows+(A.nRows%__SIMD_STEP))*V.nCols );

    TTimer seqTimer;

    std::cerr << "[SEQC] Exporting initial values... ";

    seqTimer.clear(); seqTimer.start();
    printCheckpoint(YOut, Y);
    seqTimer.stop();

    std::cerr << "done in " << seqTimer.usertime() << "s." << std::endl;

    double globalTime = 0.;
    Index_t j = restart;
    Index_t k = 0u;

    for (Index_t i = restart*K; i < d; ++i)
    {
        // Sequence
        std::cerr << "[SEQC] SPMM " << i << "... ";

        seqTimer.clear(); seqTimer.start();
        printProjection(WOut, U, Y);    // W_i = U * Y_i
        nextElement(p, A, Y, tY);       // Y_i = A * Y_(i-1) [= A^i * V]
        seqTimer.stop();

        std::cerr << "done in " << seqTimer.usertime() << "s." << std::endl;
        globalTime += seqTimer.usertime() / 60.;

        // Checkpoints
        if (++k == K)
        {
            k = 0u;
            ++j;

            // If breaking, create new files
            if (breaking) {
                sNumber = formatNumber(j);
                WOut.close();
                YOut.close();
                WOut.open(seqFileBase + sNumber + "." SDM_EXTENSION, std::ios::out | std::ios::binary);
                YOut.open(chkFileBase + sNumber + "." SDM_EXTENSION, std::ios::out | std::ios::binary);
            }

            std::cerr << "[SEQC] Exporting checkpoint " << j << "... ";

            seqTimer.clear(); seqTimer.start();
            printCheckpoint(YOut, Y);
            seqTimer.stop();

            std::cerr << "done in " << seqTimer.usertime() << "s." << std::endl;
            globalTime += seqTimer.usertime() / 60.;

        }

        if (stopAfter > 0 && globalTime >= stopAfter) {
            std::cerr << "/!\\ Stopping program prematurely due to -m option." << std::endl;
            break;
        }
    }

    std::cerr << "[SEQC] Generated full sequence in " << globalTime << " minutes total." << std::endl;
}

//-------------------//
//----- Getters -----//

void getRandomRightMatrix(Modulus_t p, Size_t nRows, Size_t nCols, DMatrix_t& V)
{
    std::cerr << "[GRRM] Generating random right matrix with s2=" << nCols << std::endl;

    if ((nCols % __SIMD_STEP) != 0u) {
        std::cerr << "/!\\ The SIMD step constant should divide s2." << std::endl;
        std::cerr << "    Please, use option -t to change s2, or change __USE_SIMD value." << std::endl;
        exit(EXIT_FAILURE);
    }

    V.nRows = nRows;
    V.nCols = nCols;
    V.data.resize(nRows * nCols);

    for (auto& dat : V.data) {
        dat = rand();
        dat %= p;
    }
    
}

void getRandomLeftProj(Size_t s1, Size_t m, LeftProj_t& U)
{
    // Left matrix creation
    if (s1 > m) {
        std::cerr << "/!\\ Left projection size s1 requires to be lower than matrix number of rows." << std::endl;
        std::cerr << "    Please, use option -s to change s1." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "[GRLP] Generating left projection with s1=" << s1 << std::endl;

    // Note: std::set will avoid duplication
    while (U.size() != s1)
        U.insert(rand() % m);
}

void getLengthInfo(Size_t m, Size_t n, Size_t s1, Size_t s2, Size_t& d, Size_t& K)
{
    if (d == 0u) {
        const Size_t Delta(10u);
    	d = std::ceil((double)m / s1) + std::ceil((double)n / s2) + Delta;
    	std::cerr << "[GLIN] Length of sequence is " << d << ". (Automatic, Delta=" << Delta << ')' << std::endl;
    } else {
    	std::cerr << "[GLIN] Length of sequence is " << d << ". (Forced)" << std::endl;
    }


    if (K == 0u) {
            //  K = sqrt(m);
      K = 4*sqrt(ceil((double)s1/s2)*d);
        std::cerr << "[GLIN] Generating a checkpoint each " << K << " elements. (Automatic)" << std::endl;
    } else {
        std::cerr << "[GLIN] Generating a checkpoint each " << K << " elements. (Forced)" << std::endl;
    }

}

//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::string modulus;
    std::string matrixFile, leftProjFile;
    std::string folder = "output/";
    Size_t s1 = 64u;
    Size_t s2 = 16u;
    Size_t K = 0u;
    Size_t d = 0u;
    size_t restart = 0u;
    int seed = -1;
    int64_t stopAfter = -1;
    bool breaking = false;

    static Argument as[] =
    {
        { 'q', "-q Q", "Set the field characteristic.",                     TYPE_STR,   &modulus },
        { 's', "-s S1", "Set the size of the left block.",                  TYPE_INT,   &s1 },
        { 't', "-t S2", "Set the size of the right block.",                 TYPE_INT,   &s2 },
        { 'd', "-d DIR", "The output directory.",                           TYPE_STR,   &folder },
        { 'f', "-f FILE", "Set input SRD matrix file.",                     TYPE_STR,   &matrixFile },
        { 'l', "-l FILE", "Set input left projection.",                     TYPE_STR,   &leftProjFile },
        { 'r', "-r restart", "Set the index of restart with a checkpoint.", TYPE_INT,   &restart },
        { 'g', "-g SEED", "Set the seed for randomness.",                   TYPE_INT,   &seed },
        { 'm', "-m MIN", "Will stop the program after MIN minutes.",        TYPE_INT,   &stopAfter },
        { 'K', "-K K", "Override the distance between two checkpoints.",    TYPE_INT,   &K },
        { 'S', "-S S", "Override the sequence length.",    TYPE_INT,   &d },
        { 'B', "-B", "Enables breaking checkpoints and sequence files.",    TYPE_BOOL,  &breaking },
        END_OF_ARGUMENTS
    };

    assert((s2 % (2u * __SIMD_STEP)) == 0u && "Needed as we want to exploit SIMD units.");

    FFLAS::parseArguments(argc, argv, as);
    if (seed < 0) seed = time(nullptr);
    FFLAS::writeCommandString(std::cerr << "Command line: ", as, argv[0]) << std::endl;

    // Initializing random
    std::ios_base::sync_with_stdio(false);
    ::srand(seed);

    // Timer
    if (stopAfter > 0) {
        std::cerr << "/!\\ Timer is set with -m option." << std::endl;
        std::cerr << "    The program will stop prematurely after " << stopAfter << " minutes." << std::endl;
    }

    // Init modulus
    Modulus_t p;
    getModulus(modulus, p);

    // Init matrix from stream
    SMatrix_t matrix;
    getMatrix(p, matrixFile, matrix);

    if (matrix.nRows % __SIMD_STEP) 
        std::cerr << "/!\\ Number of rows is not an exact multiple of the SIMD step." << std::endl;
    

    // Right matrix
    DMatrix_t V;
    if (restart == 0)
        getRandomRightMatrix(p, matrix.nCols, s2, V);
    else {
        std::string chkFileBase = folder + "/chk";
        std::ifstream chkStream;
        std::cerr << "[SEQ] Restarting at step " << restart << ", reading checkpoint..." << std::endl;
        chkStream.open(chkFileBase + formatNumber(restart) + "." SDM_EXTENSION, std::ios::in | std::ios::binary);
        if (!chkStream.is_open()) throw std::logic_error("Checkpoint file not found.");
        V.nRows = matrix.nCols;
        V.nCols = s2;
        V.data.resize(V.nRows * V.nCols);
        SDM_READ_MATRIX(chkStream, V.data);
        chkStream.close();
        std::cerr << "[SEQ] Checkpoint read." << std::endl;
    }
    
        

    // Left matrix
    LeftProj_t U;
    if (leftProjFile == "")
        getRandomLeftProj(s1, matrix.nRows, U);
    else
        readLeftProj(U,s1,leftProjFile);
    saveLeftProj(folder + "/leftproj.txt", U);

    // Parameters and infos
    getLengthInfo(matrix.nRows, matrix.nCols, s1, s2, d, K);
    saveInfo(folder + "/info.txt", matrixFile, p, matrix.nRows, matrix.nCols, s1, s2, d, K);

    // Sequence generation
    std::string seqFileBase = folder + "/seq";
    std::string chkFileBase = folder + "/chk";

    std::ofstream seqStream(seqFileBase + "_header." SDM_EXTENSION, std::ios::out | std::ios::binary);
    std::ofstream chkStream(chkFileBase + "_header." SDM_EXTENSION, std::ios::out | std::ios::binary);

    SDM_WRITE_HEADER(seqStream, d, s1, s2);
    SDM_WRITE_HEADER(chkStream, d/K, matrix.nRows, s2);

    PAR_BLOCK {
      std::clog<<"Generating sequence with "<<NUM_THREADS<<" threads"<<std::endl;
      generateSequence(seqFileBase, chkFileBase, p, restart, d, K, matrix, U, V, stopAfter, breaking);
    }

    return EXIT_SUCCESS;
}

