// ================================================================
// LinBox Project 1999
// Black Box iterator and container
// Classical one
// the sequence is u^t v, u^t A v, ...,  u^t A^n v,  
// Time-stamp: <15 Mar 00 17:46:32 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __THUNK_DEFAULT_H__
#define __THUNK_DEFAULT_H__

#include <lin_rand.h>
#include <lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t>
class BB_Iterator : public Base_BB_Iterator< BlackBoxDomain> {
public:
        //-- Constructors
    BB_Iterator() {} 

    BB_Iterator(BlackBoxDomain_t * BD, const Vecteur& u0) 
            : Base_BB_Iterator<BlackBoxDomain>(BD), even(1) { init(u0,u0); }
    BB_Iterator(BlackBoxDomain_t * BD, const Vecteur& u0, const Vecteur& v0) 
            : Base_BB_Iterator<BlackBoxDomain>(BD), even(1) { init(u0,v0); }

private:

        /// User Left and Right vectors 
    Type_t& init(const Vecteur& uu, const Vecteur& vv) {
        u = uu;
        w = vv;
        v.resize(_BB_domain->n_row());
        return DOTPROD(_value, u, u);
    }
    
private:
    bool even;
    Vecteur u, v, w;
    
    void _launch () {
        if (even) {
            _BB_domain->Apply( v, w);  // GV
            DOTPROD(_value,u,v);       // GV 
            even = 0;
        } else {
            _BB_domain->Apply( w, v);  // GV
            DOTPROD(_value,u,w);       // GV
            even = 1;
        }  
    }
    
    void _wait () {}
};

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t, class RandIter = Random>
class BB_Container : public Base_BB_Container< BlackBoxDomain> {
public:
    BB_Container() {} 

    BB_Containers(BlackBoxDomain_t * D, const Vecteur& u0) 
            : Base_BB_Container< BlackBoxDomain>(D) {}
    
    BB_Containers(BlackBoxDomain_t * D, const Vecteur& u0, const Vecteur& v0) 
            : Base_BB_Container< BlackBoxDomain>(D) {}
    
    BB_Containers(BlackBoxDomain_t * D, RandIter& g) 
            : Base_BB_Container< BlackBoxDomain>(D) {}
    
protected:
        /// User Left and Right vectors 
    Type_t& init(const Vecteur& uu, const Vecteur& vv) {
        u = uu;
        w = vv;
        v.resize(_BB_domain->n_row());
        return DOTPROD(_value, u, u);
    }
    
    bool even;
    Vecteur u, v, w;
    
    void _launch () {
        if (even) {
            _BB_domain->Apply( v, w);  // GV
            DOTPROD(_value,u,v);       // GV 
            even = 0;
        } else {
            _BB_domain->Apply( w, v);  // GV
            DOTPROD(_value,u,w);       // GV
            even = 1;
        }  
    }
    
    void _wait () {}
};


#endif
