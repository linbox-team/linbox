// ================================================================
// LinBox Project 1999
// Symmetrizing iterator (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is u^t u, (A u)^t (A u) = u^t (A^t A) u, 
// (A^t (A u))^t (A^t (A u)) = u^t (A^t A)^2 u , etc.
// Time-stamp: <10 Mar 00 18:49:00 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __THUNK_SYMMETRIZE_H__
#define __THUNK_SYMMETRIZE_H__
#include <lin_rand.h>
#include <lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t>
class BB_Symmetrize_Iterator : public Base_BB_Iterator< BlackBoxDomain > {
public:
        //-- Constructors
    BB_Symmetrize_Iterator() {} 

    BB_Symmetrize_Iterator(BlackBoxDomain_t * BD, const Vecteur& u0)
            : Base_BB_Iterator<BlackBoxDomain>(BD), even(1), u(u0) {
        v.resize(_BB_domain->n_row());
        DOTPROD(_value, u, u);
    }

private:
    bool even;
    Vecteur u, v;
    
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

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t, class RandIter = Random>
class BB_Symmetrize_Container : public Base_BB_Container< BlackBoxDomain> {
    Vecteur u;
public:
    BB_Symmetrize_Container() {} 

    BB_Symmetrize_Container(BlackBoxDomain_t * D, const Vecteur& u0) 
            : Base_BB_Container< BlackBoxDomain>(D), u(u0) {}
    
    BB_Symmetrize_Container(BlackBoxDomain_t * D, RandIter& g ) 
            : Base_BB_Container< BlackBoxDomain>(D) {
        u.resize(_BB_domain->n_col());
        for(long i=u.size();i--;)
            _domain.random(g, u[i]);
    }
    
    typedef BB_Symmetrize_Iterator< BlackBoxDomain, Vecteur > const_iterator;
    const_iterator& begin() const { return *(new const_iterator(_BB_domain,u)); }
    const_iterator& end() const { return *(new const_iterator()); }

};

#endif
