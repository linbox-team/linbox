// ================================================================
// LinBox Project 1999
// Symmetrizing iterator (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is u^t u, (A u)^t (A u) = u^t (A^t A) u, 
// (A^t (A u))^t (A^t (A u)) = u^t (A^t A)^2 u , etc.
// Time-stamp: <06 Mar 00 19:24:27 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __THUNK_SYMMETRIZE_H__
#define __THUNK_SYMMETRIZE_H__
#include <lin_dom_thunk_bbase.h>

template<class BlackBoxDomain, class Vecteur = BlackBoxDomain::PreferredInMatrix_t>
class BB_Symmetrize_Iterator : public Base_BB_Iterator< BlackBoxDomain> {
public:
        //-- Constructors
    BB_Symmetrize_Iterator() {} 

    BB_Symmetrize_Iterator(BlackBoxDomain_t * BD) 
            : Base_BB_Iterator<BlackBoxDomain>(BD) { init(); }
    BB_Symmetrize_Iterator(BlackBoxDomain_t * BD, const Vecteur& u0) 
            : Base_BB_Iterator<BlackBoxDomain>(BD) { init(u0); }

private:
    bool even;
    Vecteur u, v;
    
    Type_t& init() {
        even = 1;
        u.resize(_BB_domain->n_col());
        for(long i=u.size();i--;)
            _domain.random(u[i]);
        v.resize(_BB_domain->n_row());
        return DOTPROD(_value, u, u);
    }        

    Type_t& init(const Vecteur& uu) {
        even = 1;
        u = uu;
        v.resize(_BB_domain->n_row());
        return DOTPROD(_value, u, u);
    }
    
    void _launch () {
        if (even) {
            even = 0;
            _BB_domain->Apply(v, u);
            DOTPROD(_value,v,v); 
        } else {
            even = 1;
            _BB_domain->ApplyTrans( u, v); 
            DOTPROD(_value,u,u);
        }  
    }
    
    void _wait () {}
};

template<class BlackBoxDomain, class Vecteur = BlackBoxDomain::PreferredInMatrix_t>
class BB_Symmetrize_Container : public Base_BB_Container< BlackBoxDomain> {
public:
    typedef BB_Symmetrize_Iterator< BlackBoxDomain, Vecteur > const_iterator;
};

#endif
