#ifndef __B_B_DOM_COMPOSE_AND_SUB_H_
#define __B_B_DOM_COMPOSE_AND_SUB_H_

// ======================================================
// LinBox Project 1999
//
// Black Blox Identity and Composition
// An example of generic algorithm
// file : lin_dom_compminus.C
// ======================================================

template <class BBA, class BBB, class Inter>
class IdentityCompositionDom {
public:
    typedef typename BBB::PreferredInMatrix_t           PreferredInMatrix_t;
    typedef typename BBA::PreferredOutMatrix_t          PreferredOutMatrix_t;
    typedef          IdentityCompositionDom<BBA,BBB,Inter>      Self_t;
        // Coefficient space, BBA's for instance
    typedef typename BBA::Domain_t                      Domain_t;
    typedef typename BBA::Type_t                        Type_t;
    typedef typename BBA::coefficientSpace              coefficientSpace;
private:
        // IdentityComposition is  BBA x BBB - lambda.I
    BBA         * _amat;
    BBB         * _bmat;
    Inter       * _vect;
    Domain_t    _domain;
    Type_t      _lambda;

public:
        // Transpose information
//     typedef BB_with_transpose BB_info;
   
        //-- Default cstors:
    IdentityCompositionDom() {};

    IdentityCompositionDom(BBA * a, BBB* b, Inter * v, const Domain_t& D, Type_t l) 
        : _amat(a), _bmat(b), _vect(v), _domain(D), _lambda(l) {
                // Check sizes
        if ( (_amat->n_col() != _bmat->n_row()) || 
             (_amat->n_col() != _vect->size()))
            cerr << "Exception : matrices are not of same size\n"
                 << _amat->n_row() << "x" << _amat->n_col() 
                 << " times " <<  _vect->size() 
                 << " times " << _bmat->n_row() << "x" << _bmat->n_col() << endl;
            // Check coefficient spaces
        typename BBA::coefficientSpace xxx( typename BBB::coefficientSpace() );
    }
    
        // Size
    long n_row() const { return _amat->n_row();};
    long n_col() const { return _bmat->n_col();};
    
    Domain_t getdomain() { return _bmat->getdomain(); }
    Type_t getlambda() { return _lambda; }

    template<class OutMatrix, class InMatrix>
    inline OutMatrix& Apply(OutMatrix& outM,  const InMatrix& inM ) {
            // Need the preferred intermediate matrix type 
            // allocate intermediate with the good size
//         Domain_t _domain = _bmat->getdomain();
        _amat->Apply (outM, _bmat->Apply ( *_vect, inM) );
        for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
            _domain.amxyin( outM[ll], _lambda, inM[ll]);
// //         { 
// //             cerr << "ll:" << ll << endl;
// //             cerr << "o:" << outM[ll] << "  l:" << _lambda << "  i:" << inM[ll] << endl;
// //             cerr << "fo:" << _domain.access( outM[ll] ) << "  fl:" << _domain.access(_lambda) << "  fi:" << _domain.access(inM[ll]) << endl;

//             Type_t ktmp, utmp;
//             _domain.mul(utmp, _lambda, inM[ll]);
//             _domain.sub(ktmp, outM[ll], utmp);
//             outM[ll] = ktmp;
// //             cerr << "r:" << outM[ll] << endl;
// //             cerr << "fr:" << _domain.access(outM[ll]) << endl;
//         }
        return outM;
    }
    
    template<class OutMatrix, class InMatrix>
    inline OutMatrix& ApplyTrans(OutMatrix& outM,  const InMatrix& inM ) {
            // Compiles if BBA and BBB have ApplyTrans, or if no call to ApplyTrans is made
//         Domain_t _domain = _amat->getdomain();
        _bmat->ApplyTrans (outM, _amat->ApplyTrans ( *_vect, inM) ) ;
        for(Indice ll=GIVMIN(outM.size(), inM.size()); ll--;)
            _domain.amxyin( outM[ll], _lambda, inM[ll]);
//         { 
// //             cerr << "ll:" << ll << endl;
// //             cerr << "o:" << outM[ll] << "  l:" << _lambda << "  i:" << inM[ll] << endl;
// //             cerr << "fo:" << _domain.access( outM[ll] ) << "  fl:" << _domain.access(_lambda) << "  fi:" << _domain.access(inM[ll]) << endl;
//             Type_t ktmp, utmp;
//             _domain.mul(utmp, _lambda, inM[ll]);
//             _domain.sub(ktmp, outM[ll], utmp);
//             outM[ll] = ktmp;
// //             cerr << "r:" << outM[ll] << endl;
// //             cerr << "fr:" << _domain.access(outM[ll]) << endl;
//         }
        
        return outM;
    }
};



#endif __B_B_DOM_COMPOSE_AND_SUB_H_
