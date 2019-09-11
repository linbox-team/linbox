#include "linbox/linbox-config.h"
#include "linbox/util/commentator.h"
#include "test-common.h"

#include <iostream>
#include <fstream>
#include <cstdio>

#include "linbox/matrix/matrix-domain.h"
#include "linbox/blackbox/toeplitz.h"

using namespace LinBox;

typedef NTL_zz_pX PolynomialRing;
typedef typename PolynomialRing::Element Polynomial;
typedef typename PolynomialRing::Coeff Coeff;

typedef typename PolynomialRing::CoeffField Field;
typedef typename Field::Element Element;
typedef typename Field::RandIter RandIter;

typedef BlasMatrixDomain<Field> MatrixDom;
typedef typename MatrixDom::OwnMatrix Matrix;

class Helper {
    Field _F;
    PolynomialRing _R;
    RandIter _RI;
    
public:
    Helper(const Field &F, const PolynomialRing &R) : _F(F), _R(R), _RI(F) {}
    
    void printMatrix(Matrix &M) const {
        for (size_t i = 0; i < M.rowdim(); i++) {
            for (size_t j = 0; j < M.coldim(); j++) {
                Element e;
                M.getEntry(e, i, j);
                _F.write(std::cout, e) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    void randomPolynomial(Polynomial &f, size_t degree) const {
        _R.assign(f, _R.zero);
		
		for (size_t i = 0; i <= degree; i++) {
			Coeff c;
			_RI.random(c);
			_R.setCoeff(f, i, c);
		}
    }
    
    void randomVector(BlasVector<Field> &V) const {
        for (size_t i = 0; i < V.size(); i++) {
            Element e;
            _RI.random(e);
            V.setEntry(i, e);
        }
    }
    
    void initToeplitzMatrix(Matrix &M, const Polynomial &f) const {        
        for (size_t i = 0; i <  M.rowdim(); i++) {
            for (size_t j = 0; j < M.coldim(); j++) {
                Coeff c;
                _R.getCoeff(c, f, M.coldim() + i - j - 1);
                M.setEntry(i, j, c);
            }
        }
    }
    
    bool testApply(size_t rowdim, size_t coldim) const {
        Polynomial f;
        randomPolynomial(f, rowdim + coldim - 2);
        
        Toeplitz<Field, PolynomialRing> T(_R, f, rowdim, coldim);
        Matrix M(_F, rowdim, coldim);
        initToeplitzMatrix(M, f);
        
        BlasVector<Field> vIn(_F, coldim);
        BlasVector<Field> vOut1(_F, rowdim), vOut2(_F, rowdim);
        randomVector(vIn);
        
        M.apply(vOut1, vIn);
        T.apply(vOut2, vIn);
        
        VectorDomain<Field> VD(_F);
        if (!VD.areEqual(vOut1, vOut2)) {
            return false;
        }
        
        return true;
    }
    
    bool testApplyTranspose(size_t rowdim, size_t coldim) const {
        Polynomial f;
        randomPolynomial(f, rowdim + coldim - 2);
        
        Toeplitz<Field, PolynomialRing> T(_R, f, rowdim, coldim);
        Matrix M(_F, rowdim, coldim);
        initToeplitzMatrix(M, f);
        
        BlasVector<Field> vIn(_F, rowdim);
        BlasVector<Field> vOut1(_F, coldim), vOut2(_F, coldim);
        randomVector(vIn);
        
        M.applyTranspose(vOut1, vIn);
        T.applyTranspose(vOut2, vIn);
        
        VectorDomain<Field> VD(_F);
        if (!VD.areEqual(vOut1, vOut2)) {
            return false;
        }
        
        return true;
    }
};

int main(int argc, char **argv) {
    
    size_t rowdim = 3, coldim = 4;
    uint64_t p = 10007;
    
    static Argument args[] = {
        { 'n', "-n N", "Set row dimension of test matrices.", TYPE_INT, &rowdim },
        { 'm', "-m M", "Set column dimension of test matrices.", TYPE_INT, &coldim },
        { 'p', "-p P", "Set field cardinality.", TYPE_INT, &p },
        END_OF_ARGUMENTS
    };
    parseArguments (argc, argv, args);
    
    commentator().start("Toeplitz test suite", "toeplitz");
    bool pass = true;
    
    Field F(p);
    PolynomialRing R(F);
    
    Helper H(F, R);
    
    pass = pass && H.testApply(rowdim, coldim);
    pass = pass && H.testApplyTranspose(rowdim, coldim);
    
	commentator().stop(MSG_STATUS (pass),"Toeplitz test suite");
    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
