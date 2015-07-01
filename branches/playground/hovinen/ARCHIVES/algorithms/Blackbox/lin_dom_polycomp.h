#ifndef __B_B_DOM_POLY_COMPOSE_H_
#define __B_B_DOM_POLY_COMPOSE_H_

// ======================================================
// LinBox Project 1999
//
// Time-stamp: <08 Sep 99 19:12:40 Jean-Guillaume.Dumas@imag.fr> 
//
// Black Blox Polynomial Composition
// An example of generic algorithm
// file : lin_dom_polycomp.h
// ======================================================

template <class BBA, class Inter, class Poly_t, class Diag_t>
class PolynomialCompositionDom {
public:
    typedef typename BBA::PreferredInMatrix_t                           PreferredInMatrix_t;
    typedef typename BBA::PreferredOutMatrix_t                          PreferredOutMatrix_t;
    typedef          PolynomialCompositionDom<BBA,Inter,Poly_t,Diag_t>  Self_t;
    typedef typename BBA::Domain_t                                      Domain_t;
    typedef typename BBA::Type_t                                        Type_t;
    typedef typename BBA::coefficientSpace                              coefficientSpace;
private:
    Domain_t    _domain;
        // PolynomialComposition is  _poly(BBA) where _poly[0] * Id = _diag
    BBA         * _amat;
    Inter       * _vect1;
    Inter       * _vect2;
    Poly_t      * _poly;
    Diag_t      * _diag;

public:
        // Transpose information
//     typedef BB_with_transpose BB_info;
   
        //-- Default cstors:
    PolynomialCompositionDom() {};

    PolynomialCompositionDom(BBA * a, Inter * v1, Inter * v2, Poly_t * pol, Diag_t * diag, const Domain_t& D) 
        : _amat(a),  _vect1(v1),  _vect2(v2), _poly(pol), _diag(diag), _domain(D) {
                // Check sizes

        long n = _amat->n_col();
        if ( (_amat->n_row() != n) || 
             (_diag->size()  != n) ||
             (_vect1->size() != n) ||
             (_vect2->size() != n) 
             )
            cerr << "Exception : matrices are not of same size\n"
                 << _vect2->size() 
                 << " times " 
                 << _amat->n_row() << "x" << _amat->n_col() 
                 << " times " 
                 << _vect1->size() 
                 << " plus " 
                 << _diag->size() 
                 << endl;
    }
    
        // Size
    long n_row() const { return _amat->n_row();};
    long n_col() const { return _amat->n_col();};
    
    Domain_t getdomain() { return _domain; }

    template<class OutMatrix, class InMatrix>
    inline OutMatrix& Apply(OutMatrix& outM,  const InMatrix& inM ) {
        *_vect1 = inM;
        for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
            _domain.mul(outM[ll],(*_diag)[ll], inM[ll] );
        Indice j = 0;
        for(++j;j<_poly->size();++j) {
            _amat->Apply(*_vect2, *_vect1);
            for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
                _domain.axpyin(outM[ll],(*_poly)[j],(*_vect2)[ll]);
            if (++j >= _poly->size() ) break;
            _amat->Apply(*_vect1, *_vect2);
            for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
                _domain.axpyin(outM[ll],(*_poly)[j],(*_vect1)[ll]);   
        }
        return outM;
    }
    
    template<class OutMatrix, class InMatrix>
    inline OutMatrix& ApplyTrans(OutMatrix& outM,  const InMatrix& inM ) {
        *_vect1 = inM;
        for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
            outM[ll] = (*_diag)[ll];
        Indice j = 0;
        for(++j;j<_poly->size();++j) {
            _amat->ApplyTrans(*_vect2, *_vect1);
            for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
                _domain.axpyin(outM[ll],(*_poly)[j],(*_vect2)[ll]);
            if (++j >= _poly->size() ) break;
            _amat->ApplyTrans(*_vect1, *_vect2);
            for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
                _domain.axpyin(outM[ll],(*_poly)[j],(*_vect1)[ll]);   
        }
        return outM;
    }
    
};



#endif __B_B_DOM_POLY_COMPOSE_H_
