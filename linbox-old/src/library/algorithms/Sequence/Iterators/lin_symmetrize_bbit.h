// ================================================================
// LinBox Project 1999
// Symmetrizing iterator (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is u^t u, (A u)^t (A u) = u^t (A^t A) u, 
// (A^t (A u))^t (A^t (A u)) = u^t (A^t A)^2 u , etc.
// Time-stamp: <15 May 00 15:21:57 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __BBContainer_SYMMETRIZE_H__
#define __BBContainer_SYMMETRIZE_H__
#include <lin_rand.h>
#include <lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t, class RandIter = Random>
class BB_Symmetrize_Container : public Base_BB_Container< BlackBoxDomain, Vecteur > {
public:
    BB_Symmetrize_Container() {} 

    BB_Symmetrize_Container(BlackBoxDomain_t * D, const Vecteur& u0) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(u0); }
    
    BB_Symmetrize_Container(BlackBoxDomain_t * D, RandIter& g ) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(g); }
    
private:
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

#endif // __BBContainer_SYMMETRIZE_H__
