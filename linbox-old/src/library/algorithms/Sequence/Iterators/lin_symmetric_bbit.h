// ================================================================
// LinBox Project 1999
// Black Box iterator and container 
// For symmetric matrix with same left and right vector
// the sequence is u^t v, u^t A v, ...,  u^t A^n v,  
// Time-stamp: <15 Mar 00 18:18:46 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __BBContainer_SYMMETRIC_H__
#define __BBContainer_SYMMETRIC_H__

#include <LinBox/lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t>
class BB_Symmetric_Container : public Base_BB_Container< BlackBoxDomain, Vecteur > {
public:
    BB_Symmetric_Container() {} 

    BB_Symmetric_Container(BlackBoxDomain_t * D, const Vecteur& u0) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(u0,u0); }
    
    BB_Symmetric_Container(BlackBoxDomain_t * D, RandIter& g ) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(g); }

protected:
    void _launch () {
        switch(even) {
            case 1:
                _BB_domain->Apply( v, u);  // v <- B(B^i u_0) = B^(i+1) u_0
                DOTPROD(_value,u,v);       // t <- u^t v = u_0^t B^(2i+1) u_0
                even = 2;
                break;
            case 2:
                DOTPROD(_value,v,v);       // t <- v^t v = u_0^t B^(2i+2) u_0
                even = 3;
                break;
            case 3:
                _BB_domain->Apply( u, v);  // u <- B(B^(i+1) u_0) = B^(i+2) u_0
                DOTPROD(_value,v,u);       // t <- v^t u = u_0^t B^(2i+3) u_0
                even = 0;
                break;
            default:
                DOTPROD(_value,u,u);       // t <- u^t u = u_0^t B^(2i+4) u_0
                even = 1;
                break;
        }
    }
    
    void _wait () {}
};


#endif // __BBContainer_SYMMETRIC_H__
