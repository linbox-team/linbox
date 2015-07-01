#ifndef __B_B_DOM_TRANSPOSE_H_
#define __B_B_DOM_TRANSPOSE_H_

// ======================================================
// LinBox Project 1999
//
// Black Blox Transposition
// An example of generic algorithm
// file : lin_dom_trans.h
// ======================================================

template <class BlackBoxDomain>
class TransposeDom {
public:
    typedef          BlackBoxDomain                             BlackBoxDomain_t;
    typedef typename BlackBoxDomain::PreferredInMatrix_t        PreferredInMatrix_t;
    typedef typename BlackBoxDomain::PreferredOutMatrix_t       PreferredOutMatrix_t;
    typedef typename BlackBoxDomain::Domain_t                   Domain_t;
    typedef typename BlackBoxDomain::Type_t                     Type_t;
    typedef typename BlackBoxDomain::coefficientSpace           coefficientSpace;
    typedef          TransposeDom< BlackBoxDomain >             Self_t;
private:
    BlackBoxDomain_t * _BB_domain;

public:
        // Transpose information
//     typedef BB_with_transpose BB_info;
   
        //-- Default cstors:
    TransposeDom() {};
    TransposeDom(BlackBoxDomain_t * D) : _BB_domain(D) {}
    
        //-- Cstor of recopy: compiler's generated
    TransposeDom(const Self_t& M) : _BB_domain(M._BB_domain) {}

        //-- Usefull to use the same domain to perform other operations
    Domain_t getdomain() const { return _BB_domain->getdomain(); }
    BlackBoxDomain_t * getBBdomain() const { return _BB_domain; }

        // Size
    long n_row() const { return _BB_domain->n_col();};
    long n_col() const { return _BB_domain->n_row();};


    template<class OutMatrix, class InMatrix>
    OutMatrix& Apply(OutMatrix& res, const InMatrix& vect ) {
        return _BB_domain->ApplyTrans(res, vect);
    }
    
    template<class OutMatrix, class InMatrix>
    OutMatrix& ApplyTrans(OutMatrix& res, const InMatrix& vect ) {
        return _BB_domain->Apply(res, vect);
    }
    
};




#endif __B_B_DOM_TRANSPOSE_H_
