#ifndef __B_B_DOM_AAT_H
#define __B_B_DOM_AAT_H

// ======================================================
// LinBox Project 1999
// ======================================================

template <class BlackBoxDomain, class Inter = typename BlackBoxDomain::PreferredInMatrix_t>
class AAT {
public:
    typedef          BlackBoxDomain                             BlackBoxDomain_t;
    typedef typename BlackBoxDomain::PreferredInMatrix_t        PreferredInMatrix_t;
    typedef typename BlackBoxDomain::PreferredOutMatrix_t       PreferredOutMatrix_t;
    typedef typename BlackBoxDomain::Domain_t                   Domain_t;
    typedef typename BlackBoxDomain::Type_t                     Type_t;
    typedef          AAT< BlackBoxDomain >             Self_t;
private:
    BlackBoxDomain_t * _BB_domain;
    Inter _inter;

public:
        // Transpose information
//     typedef BB_with_transpose BB_info;
   
        //-- Default cstors:
    AAT() {};
    AAT(BlackBoxDomain_t * D) : _BB_domain(D), _inter(D->n_col()) {}
    
        //-- Cstor of recopy: compiler's generated
    AAT(const Self_t& M) : _BB_domain(M._BB_domain), _inter(M._inter) {}

        //-- Usefull to use the same domain to perform other operations
    Domain_t getdomain() const { return _BB_domain->getdomain(); }
    BlackBoxDomain_t * getBBdomain() const { return _BB_domain; }

        // Size
    long n_row() const { return _BB_domain->n_row();};
    long n_col() const { return _BB_domain->n_row();};


    template<class OutMatrix, class InMatrix>
    OutMatrix& Apply(OutMatrix& res, const InMatrix& vect) {
        return _BB_domain->Apply(res, _BB_domain->ApplyTrans( _inter, vect) );
    }
    
    template<class OutMatrix, class InMatrix>
    OutMatrix& ApplyTrans(OutMatrix& res, const InMatrix& vect ) {
        return _BB_domain->Apply(res, _BB_domain->ApplyTrans( _inter, vect) );
    }
    
};




#endif // __B_B_DOM_AAT_H
