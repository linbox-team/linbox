// ================================================================
// LinBox Project 1999
// Black Box iterator and container
// Classical one
// the sequence is u^t v, u^t A v, ...,  u^t A^n v.  
// Time-stamp: <15 Mar 00 18:20:44 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __BBContainer_Default_H__
#define __BBContainer_Default_H__
#include <LinBox/lin_rand.h>
#include <LinBox/lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t, class RandIter = Random>
class BB_Container : public Base_BB_Container< BlackBoxDomain, Vecteur> {
public:
    BB_Container() {} 

    BB_Containers(BlackBoxDomain_t * D, const Vecteur& u0) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(u0,u0); w=u; }
    
    BB_Containers(BlackBoxDomain_t * D, const Vecteur& u0, const Vecteur& v0) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(u0,v0); w=u;}
    
    BB_Containers(BlackBoxDomain_t * D, RandIter& g) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(g); w=u; }
    
protected:
    Vecteur w;
    
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


#endif // __BBContainer_Default_H__
