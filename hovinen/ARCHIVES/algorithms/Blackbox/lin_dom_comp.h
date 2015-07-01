#ifndef __B_B_DOM_COMPOSE_H_
#define __B_B_DOM_COMPOSE_H_

// ======================================================
// LinBox Project 1999
//
// Black Blox Transposition
// An example of generic algorithm
// file : lin_dom_trans.h
// ======================================================

template <class BBA, class BBB, class Inter>
class CompositionDom {
public:
    typedef typename BBB::PreferredInMatrix_t           PreferredInMatrix_t;
    typedef typename BBA::PreferredOutMatrix_t          PreferredOutMatrix_t;
    typedef          CompositionDom<BBA,BBB,Inter>      Self_t;
        // Coefficient space, BBA's for instance
    typedef typename BBA::Domain_t                      Domain_t;
    typedef typename BBA::Type_t                        Type_t;
    typedef typename BBA::coefficientSpace              coefficientSpace;
private:
        // Composition is  BBA x BBB
    BBA   * _amat;
    BBB   * _bmat;
    Inter * _vect;

public:
        // Transpose information
//     typedef BB_with_transpose BB_info;
   
        //-- Default cstors:
    CompositionDom() {};

    CompositionDom(BBA * a, BBB* b, Inter * v) 
        : _amat(a), _bmat(b), _vect(v) {
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
    
    Domain_t getdomain() { return _amat->getdomain(); }


    template<class OutMatrix, class InMatrix>
    inline OutMatrix& Apply(OutMatrix& outM,  const InMatrix& inM ) {
            // Need the preferred intermediate matrix type 
            // allocate intermediate with the good size
        return _amat->Apply (outM, _bmat->Apply ( *_vect, inM) );
    }
    
    template<class OutMatrix, class InMatrix>
    inline OutMatrix& ApplyTrans(OutMatrix& outM,  const InMatrix& inM ) {
            // Compiles if BBA and BBB have ApplyTrans, or if no call to ApplyTrans is made
        return _bmat->ApplyTrans (outM, _amat->ApplyTrans ( *_vect, inM) ) ;
    }
};



#endif __B_B_DOM_COMPOSE_H_
