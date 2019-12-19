/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include <fflas-ffpack/utils/args-parser.h>

#include "verif_read.hpp"
#include "verif_define.hpp"
#include "common_sdm.hpp"

//-----------------//
//----- Utils -----//

void saveMAT(const Field_t& field, const std::string& matrixFile, const DMatrix_t& matrix,
             const Index_t m, const Index_t n)
{
    std::ofstream matrixStream(matrixFile);
    matrixStream << m << " " << n << " M" << std::endl;
    for (const auto& x : matrix)
        field.write(matrixStream, x) << std::endl;
}

void saveList(const Field_t& field, const std::string& filename, const std::vector<FieldElement_t>& r)
{
    std::ofstream ofs(filename);

    ofs << "Size: " << r.size() << std::endl;
    for (const auto& ri : r)
        field.write(ofs, ri) << std::endl;
}

//------------------------//
//----- Precomputing -----//

TTimer precomputeZ(DMatrix_t& Z, const DMatrix_t& P, const SMatrix_t& A,
                 Index_t K, const Field_t& field, const Index_t m)
{
    TTimer timer;
    timer.clear();
    timer.start();

    Z = P; // Copy
    DMatrix_t tZ(Z.size());

    PAR_BLOCK {
        std::cerr << "[PREZ] NUMTHREADS=" << NUM_THREADS << " Precomputing Z = P.A^K..." << std::endl;    
        
        
            // Computation
        for (Index_t i = 0u; i < K; ++i) {
//             pfspmm(field, A, 1u, Z.data(), 1u, field.zero, tZ.data(), 1u);
            pfspmv(field, A, Z.data(), field.zero, tZ.data());
            std::swap(Z, tZ); // Z <- tZ
        }
    }
    

    timer.stop();
    std::cerr << "[PREZ] /!\\ Precomputing Z done in " << timer << "." << std::endl;

    return timer;
}

DMatrix_t& applyProj(DMatrix_t& V, const Proj_t& U, const DMatrix_t& r,
                     const Index_t k, const Index_t s1, const Field_t& field) {
    const FieldElement_t * rk = r.data()+s1*k;
    size_t i=0; for(const auto& Uj : U) {
        field.assign(V[Uj], rk[i]); ++i;
    }
    return V;
}

DMatrix_t& resetProj(DMatrix_t& V, const Proj_t& U, const Field_t& field) {
    for(const auto& Uj : U) {
        field.assign(V[Uj], field.zero);
    }
    return V;
}

//! T = AT x UT is a 1 x m  row-major
TTimer precomputeT(DMatrix_t& T, const Proj_t& U, const SMatrix_t& A, const DMatrix_t& r,
                   const Index_t K, const Index_t m, const Index_t s1,
                   const Field_t& field)
{
    TTimer timer;
    timer.clear();
    timer.start();

    
    PAR_BLOCK {
        std::cerr << "[PRET] NUMTHREADS=" << NUM_THREADS << " Precomputing T = \\sum_k U * A^k * r_k..." << std::endl;
        DMatrix_t tk(T.size(),field.zero);
        DMatrix_t tU(T.size(),field.zero);
        
        applyProj(T, U, r, K-1, s1, field); // T = U (r[K-1])^T 
        
        for(size_t k=1 ; k<K; ++k) {
//             pfspmm(field, A, 1, T.data(), 1, field.zero, tU.data(), 1); // tU = A * T
            pfspmv(field, A, T.data(), field.zero, tU.data()); // tU = A * T
            
            applyProj(tk, U, r, K-1-k, s1, field);
            FFLAS::pfaddin(field, m, 1, tk.data(), 1, tU.data(), 1,NUM_THREADS); // tU += tk
            resetProj(tk, U, field);
            
            std::swap(tU,T);
        }
        
    }

    timer.stop();
    std::cerr << "[PRET] /!\\ Precomputing T done in " << timer << "." << std::endl;
    return timer;
}

//-------------------//
//----- Getters -----//

void getMatrix(const Field_t& field, const std::string& matrixFile, SMatrix_t& matrix,
               Index_t& rowDim, Index_t& colDim, Index_t& nnz, bool transposed)
{

    // Select stream
    std::ifstream matrixStream(matrixFile);
    if (!matrixStream.is_open()) {
        std::cerr << "/!\\ This matrix file is not valid: " << matrixFile << "." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "[GMAT] Reading matrix from file " << matrixFile << "..." << std::endl;

    TTimer timer;
    timer.clear(); timer.start();

    auto extension = matrixFile.substr(matrixFile.find_last_of(".") + 1);
    if (extension == "srd")
        readSRDMatrix(field, matrixStream, matrix, rowDim, colDim, nnz, transposed);
    else if (extension == "sms")
        readSMSMatrix(field, matrixStream, matrix, rowDim, colDim, nnz, transposed);
    else {
        std::cerr << "/!\\ Unknown extension: " << extension << std::endl;
        std::cerr << "    Please use lowercase and sms or srd matrix file format." << std::endl;
        exit(EXIT_FAILURE);
    }

    timer.stop();

    std::cerr << "[GMAT] Read matrix " << rowDim << "x" << colDim << " in ";
    std::cerr << timer.usertime() << "s." << std::endl;
    std::cerr << "[GMAT] Matrix has " << nnz << " nnz." << std::endl;
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
        { 'g', "-g SEED", "Set the seed for randomness.",                   TYPE_INT,   &seed },
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
    const Field_t::Residu_t q(signedq);

    // Construct field
    std::ios_base::sync_with_stdio(false);
    ::srand(seed);
    RecInt::srand(seed);
    Field_t field(q);
    Field_t::RandIter g(field);

    // Read the matrix
    SMatrix_t AT;
    Size_t nnz = 0u;
    Index_t rowDim, colDim;
    getMatrix(field, matrixFile, AT, rowDim, colDim, nnz, true);
    assert(rowDim == n && colDim == m);

    //----- Checkpoints test -----//

    // Generate P of size 1 * m fully random
    DMatrix_t P(1u * m);
    for (auto& p : P)
        g.random(p);

    saveMAT(field, folder + "/verif_p.mat", P, 1u, m);

    // Precompute Z = P * A^K
    DMatrix_t Z(P.size());
    precomputeZ(Z, P, AT, K, field, m); // ZT = AT^K.PT
    // This works with the transpose because P.nRows == 1u

    saveMAT(field, folder + "/verif_z.mat", Z, 1u, m);

    //----- Sequence test -----//

    // Generate (r) of fully random scalars
    DMatrix_t r(K*s1);
    for (auto& ri : r)
        g.random(ri);
    saveMAT(field, folder + "/verif_r.mat", r, K, s1);

    // Read U
    Proj_t U;
    std::ifstream projStream(folder + "/leftproj.txt");
    if (!projStream.is_open())
        throw std::runtime_error("Cannot open projection file.");
    getProj(projStream, U, s1);

    // Precompute T = \sum_k U * A^k * r_k
    DMatrix_t T(1u * m, 0u);

    precomputeT(T, U, AT, r, K, m, s1, field);

    saveMAT(field, folder + "/verif_t.mat", T, m, s1);

    return EXIT_SUCCESS;
}

