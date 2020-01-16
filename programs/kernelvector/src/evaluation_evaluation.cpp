/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include "evaluation_gemv.hpp"
//------------------//
//----- Kernel -----//

inline void getNextGeGeneric(std::istream& polyStream, DVector_t& Ge, const DVector_t& e)
{
    const Size_t s2 = e.size();
    DMatrix_t G(s2, s2);

    // Get G
    SDM_READ_MATRIX(polyStream, G.data);
    
    gemv(Ge, G, e);
}

//------------------//
//----- Kernel -----//

inline void getNextGeFirstCanonical(std::istream& polyStream, DVector_t& Ge)
{
    // Get G (G is square so stride is also its size)
    SDM_READ_VECTOR(polyStream, Ge, Ge.size());
}

template<bool isFirstCanonical>
inline void getNextGe(std::istream& polyStream, DVector_t& Ge, const DVector_t& e);

template<>
inline void getNextGe<true>(std::istream& polyStream, DVector_t& Ge, const DVector_t& e) { return getNextGeFirstCanonical(polyStream,Ge); }

template<>
inline void getNextGe<false>(std::istream& polyStream, DVector_t& Ge, const DVector_t& e) { return getNextGeGeneric(polyStream,Ge,e); }




// Requires X to be null matrix
template<bool g0isFirstCanonical>
void kernel(const Modulus_t p, DVector_t& X, const SMatrix_t& A,
            std::istream& chkStream, std::istream& polyStream,
            const Index_t K, const DVector_t& e, const bool distributed)
{
    std::cerr << "[KRNL] Computing kernel vector..." << std::endl;

    // Import headers
    Size_t ng, dK, m, s2;

    SDM_READ_HEADER(chkStream, dK, m, s2);
    SDM_READ_HEADER(polyStream, ng, s2, s2);
    
    std::cerr << "[KRNL] Polynomial length is " << ng << std::endl;

    // Storage
    DVector_t Ge(s2);   // Always G_{j*K + k} * e
    DMatrix_t AkYj(m, s2);

    // Temporaries
    DMatrix_t tAkYj(m, s2);
    DVector_t tX(X.size());

    // Real computation
    TTimer krnlTimer, globalTimer;
    globalTimer.clear(); globalTimer.start();

    // Skipping G0
    if (!distributed) {
        DMatrix_t G0(s2, s2);
        SDM_READ_MATRIX(polyStream, G0.data);
        --ng;
    }

    // Compute
    Size_t nK = static_cast<Size_t>(std::ceil(static_cast<double>(ng) / K));
    std::cerr << "[KRNL] There are " << nK << " blocks to compute." << std::endl;
    
    for (Index_t j = 0u; j < nK; ++j)
    {
        krnlTimer.clear(); krnlTimer.start();
    
        // AkYj = A^0 * Y_{j*K} = A^{j*K} * V
        SDM_READ_MATRIX(chkStream, AkYj.data);
        
        krnlTimer.stop();
        std::cerr << "[KRNL] Checkpoint " << j << " read in " << krnlTimer << "." << std::endl;
        krnlTimer.clear(); krnlTimer.start();

        // As K might not divide ng...
        const Index_t lK = (j == nK)? (ng - nK * K) : K;
        std::cerr << "[KRNL] There are " << lK << " iterations in this block." << std::endl;
        for (Index_t k = 0u; k < lK; ++k)
        {
            // Get next first column of G
            getNextGe<g0isFirstCanonical>(polyStream, Ge, e);
// writeDV(std::cerr << "VGe[" << k << "]: ", Ge);

            // tX = AkYj * G_{j*K + k+1} * e
            gemv(tX, AkYj, Ge);

            // X += tX
            for (Index_t i = 0u; i < m; ++i) {
                X[i] += tX[i];
                if (X[i] >= p)
                    X[i] -= p;
            }
            
            krnlTimer.stop();
            std::cerr << "[KRNL] GEMV: " << krnlTimer << ". ";
            krnlTimer.clear(); krnlTimer.start();

            // AkYj <- A^{k+1} * Y_{jK} = A * AkYj
            spmm(p, tAkYj, A, AkYj);  // tY = A * Y
            std::swap(AkYj, tAkYj);
            
            krnlTimer.stop();
            std::cerr << "SPMM: " << krnlTimer << "." << std::endl;
            krnlTimer.clear(); krnlTimer.start();
        }
    }

    globalTimer.stop();
    std::cerr << "[KRNL] Done in " << globalTimer << "." << std::endl;
}


#include <deque>

// Requires X to be null matrix
template<bool g0isFirstCanonical>
void horner(const Modulus_t p, DVector_t& X, const SMatrix_t& A,
            std::istream& chkStream, std::istream& polyStream,
            const Index_t K, const DVector_t& e, const bool distributed)
{
    std::cerr << "[KRNL] Computing kernel vector via Horner's scheme..." << std::endl;

    // Import headers
    Size_t ng, dK, m, s2;

    SDM_READ_HEADER(chkStream, dK, m, s2);
    SDM_READ_HEADER(polyStream, ng, s2, s2);
    
    std::cerr << "[KRNL] Polynomial length is " << ng << std::endl;

    // Storage
    std::deque<DVector_t> VGe(K);   // Always G_{j*K + k} * e

    DMatrix_t AkYj(m, s2);

    // Temporaries
    DMatrix_t tAkYj(m, s2);
    DVector_t tX(X.size());

    // Real computation
    TTimer krnlTimer, globalTimer;
    globalTimer.clear(); globalTimer.start();

    // Skipping G0
    if (!distributed) {
        DMatrix_t G0(s2, s2);
        SDM_READ_MATRIX(polyStream, G0.data);
        --ng;
    }

    // Compute
    Size_t nK = static_cast<Size_t>(std::ceil(static_cast<double>(ng) / K));
    std::cerr << "[KRNL] There are " << nK << " blocks to compute." << std::endl;
    
    for (Index_t j = 0u; j < nK; ++j)
    {
        krnlTimer.clear(); krnlTimer.start();
    
        // AkYj = A^0 * Y_{j*K} = A^{j*K} * V
        SDM_READ_MATRIX(chkStream, AkYj.data);
        
        krnlTimer.stop();
        std::cerr << "[KRNL] Checkpoint " << j << " read in " << krnlTimer << "." << std::endl;
        krnlTimer.clear(); krnlTimer.start();

        // As K might not divide ng...
        const Index_t lK = (j == nK)? (ng - nK * K) : K;
        std::cerr << "[KRNL] There are " << lK << " iterations in this block." << std::endl;

        VGe.resize(0);
        for (Index_t k = 0u; k < lK; ++k)
        {
            DVector_t Ge(s2);
                // Get next first column of G
            getNextGe<g0isFirstCanonical>(polyStream, Ge, e);
            VGe.push_front(Ge);
        }
        krnlTimer.stop();
        std::cerr << "[KRNL] first " << lK << " columns read in " << krnlTimer << "." << std::endl;
        krnlTimer.clear(); krnlTimer.start();

        tX=VGe.front(); // tX = G_{j*K + K-1} * e
// writeDV(std::cerr << "VGe[" << (lK-1) << "]: ", tX);
        DVector_t tU(X.size());
            // tU = AkYj tX 
        gemv(tU, AkYj, tX);

        krnlTimer.stop();
        std::cerr << "[KRNL] GEMV: " << krnlTimer << '.' << std::endl;
        krnlTimer.clear(); krnlTimer.start();

        for(Index_t k = 1u; k < lK; ++k) {
            
            spmv(p, tX, A, tU);  // tX = A * tU
            
            krnlTimer.stop();
            std::cerr << "[KRNL] SPMV: " << krnlTimer << ". ";
            krnlTimer.clear(); krnlTimer.start();

            const DVector_t & VGek = VGe[k];
// writeDV(std::cerr << "VGe[" << (lK-1-k) << "]: ", tX);
            // tU = AkYj VGek 
            gemv(tU, AkYj, VGek);

            krnlTimer.stop();
            std::cerr << "GEMV: " << krnlTimer << ". " << std::endl;
            krnlTimer.clear(); krnlTimer.start();

            // tU += tX
            for (Index_t i = 0u; i < m; ++i) {
                tU[i] += tX[i];
                if (tU[i] >= p)
                    tU[i] -= p;
            }
        }
        
            // X += tU
        for (Index_t i = 0u; i < m; ++i) {
            X[i] += tU[i];
            if (X[i] >= p)
                X[i] -= p;
        }
        
        
        krnlTimer.stop();
        std::cerr << "[KRNL] ADDV: " << krnlTimer << "." << std::endl;

    }

    globalTimer.stop();
    std::cerr << "[KRNL] Done in " << globalTimer << "." << std::endl;
}

//----------------//
//----- Main -----//

int main(int argc, char **argv)
{
    std::string folder = "output/";
    bool distributed = false;

    static Argument as[] =
    {
        { 'd', "-d FOLDER",  "The directory where .sdm are stored.", TYPE_STR, &folder },
        { '+', "-+",  "Enable distributed mode.", TYPE_BOOL, &distributed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);
    
    if (distributed) std::cerr << "[EVAL] Distributed mode ENABLED!" << std::endl;
    else             std::cerr << "[EVAL] Distributed mode DISABLED!" << std::endl;

    // Get values
    Modulus_t p;
    Size_t m, n, s1, s2, d, K;
    std::string matrixFile;

    std::cerr << "[EVAL] Reading info from '" << folder << "/info.txt'" << std::endl;
    loadInfo(folder + "/info.txt", matrixFile, p, m, n, s1, s2, d, K);

    // Construct field
    std::ios_base::sync_with_stdio(false);
    ::srand(time(nullptr));
    RecInt::srand(time(nullptr));
    FieldElement_t::init_module(p);
    
    // Get G0
    DVector_t e; // e is such that G_0 * e == 0 and e != 0
    std::ifstream g0kernelStream(folder + "/g0kernel.dv");
    if (!g0kernelStream.is_open()) throw std::runtime_error("Kernel for G0 file g0kernel.dv not found.");
    readDV(g0kernelStream, e);

    // Init matrix from stream
    SMatrix_t A;
    getMatrix(p, matrixFile, A);

    //----- Real computing -----//
    
    std::ifstream chkStream (folder + "/chk." SDM_EXTENSION, std::ios::in | std::ios::binary);
    std::ifstream polyStream(folder + "/poly." SDM_EXTENSION, std::ios::in | std::ios::binary);
    if (!chkStream.is_open() || !polyStream.is_open()) throw std::runtime_error("Polynomial file not found.");

    // Compute kernel
    DVector_t X(m, 0u);


    if (isFirstCanonical(e)) {      
        std::cerr << "[EVAL] g0 kernel is the first canonical vector." << std::endl;
    	PAR_BLOCK {
            horner<true>(p, X, A, chkStream, polyStream, K, e, distributed);
    	}
    } else {
        std::cerr << "[EVAL] g0 kernel generic." << std::endl;
    	PAR_BLOCK {
            horner<false>(p, X, A, chkStream, polyStream, K, e, distributed);
    	}
    }
    

    // Export kernel
    std::ofstream kernelStream(folder + "/kernel.dv");
    writeDV(kernelStream, X);
    
    std::cerr << "[EVAL] The computed vector has been saved to " << folder << "/kernel.dv" << std::endl;

    return EXIT_SUCCESS;
}

