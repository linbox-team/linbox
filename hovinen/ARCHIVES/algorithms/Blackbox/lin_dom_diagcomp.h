#ifndef __B_B_DOM_Diag_COMPOSE_H_
#define __B_B_DOM_Diag_COMPOSE_H_

// ======================================================
// LinBox Project 1999
//
// Black Blox Diagonal Composition
// An example of generic algorithm
// file : lin_dom_compminus.C
// ======================================================

template <class BBA, class BBB, class Inter, class Diag_t>
class DiagonalCompositionDom {
public:
    typedef typename BBA::PreferredInMatrix_t                           PreferredInMatrix_t;
    typedef typename BBA::PreferredOutMatrix_t                          PreferredOutMatrix_t;
    typedef          DiagonalCompositionDom<BBA,BBB,Inter,Diag_t>       Self_t;
    typedef typename BBA::Domain_t                                      Domain_t;
    typedef typename BBA::Type_t                                        Type_t;
    typedef typename BBA::coefficientSpace                              coefficientSpace;
private:
    Domain_t    _domain;
        // DiagonalComposition is BBA x BBB - diag
    BBA         * _amat;
    BBB         * _bmat;
    Inter       * _vect;
    Diag_t      * _diag;

public:
        // Transpose information
//     typedef BB_with_transpose BB_info;
   
        //-- Default cstors:
    DiagonalCompositionDom() {};

    DiagonalCompositionDom(BBA * a, BBB * b, Inter * v, Diag_t * diag, const Domain_t& D) 
        : _amat(a),_bmat(b),  _vect(v), _diag(diag), _domain(D) {
                // Check sizes
        if (_amat->n_col() != _bmat->n_row())
            cerr << "Exception : matrices must be compatible\n"
                 << _amat->n_row() << "x" << _bmat->n_col() 
                 << endl;
    }
    
        // Size
    long n_row() const { return _amat->n_row();};
    long n_col() const { return _amat->n_col();};
    
    Domain_t getdomain() { return _domain; }

    template<class OutMatrix, class InMatrix>
    inline OutMatrix& Apply(OutMatrix& outM,  const InMatrix& inM ) {
        _amat->Apply(outM, _bmat->Apply(*_vect, inM) );
        for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
            _domain.amxyin(outM[ll],(*_diag)[ll],inM[ll]);   
        return outM;
    }
    
    template<class OutMatrix, class InMatrix>
    inline OutMatrix& ApplyTrans(OutMatrix& outM,  const InMatrix& inM ) {
        _bmat->ApplyTrans(outM, _amat->ApplyTrans(*_vect, inM) );
        for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
            _domain.amxyin(outM[ll],(*_diag)[ll],inM[ll]);   
        return outM;
    }
    
};



#endif __B_B_DOM_Diag_COMPOSE_H_
