// ================================================================
// LinBox Project 1999
// Symmetric Thunk
// Same left and right vector, A is supposed to be symmetric
// Therefore only half the matrix-vector products are needed
// the sequence is u^t u, u^t (A u), (A u)^t (A u) = u^t A^2 u, etc.
// Time-stamp: <05 Jan 00 12:49:41 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __THUNK_SYMMETRIC_H__
#define __THUNK_SYMMETRIC_H__
#include <lin_dom_thunk_bbase.h>

template<class BlackBoxDomain>
class BBSymmetricThunk : public BBThunkDom< BlackBoxDomain > {
public:
        //-- Constructors
    BBSymmetricThunk() {} 

    BBSymmetricThunk(BlackBoxDomain_t * BD) 
            : BBThunkDom<BlackBoxDomain>(BD) { init(); }
    
    BBSymmetricThunk(BlackBoxDomain_t * BD, const Domain_t& D) 
            : BBThunkDom<BlackBoxDomain>(BD,D) { init(); }

    Type_t& init() {
        even = 1;
        u.resize(_BB_domain->n_col());
        for(long i=u.size();i--;)
            _domain.random(u[i]);
        v.resize(_BB_domain->n_row());
        DOTPROD(_value, u, u);
    }        

    Type_t& init(const Vecteur& uu) {
        even = 1;
        u = uu;
        v.resize(_BB_domain->n_row());
        return DOTPROD(_value, u, u);
   }
    
private:
    bool even;
    Vecteur u, v;
    
    void _launch () {
        if (even) {
            _BB_domain->Apply(v,u);
            DOTPROD(_value,v,u); 
            even = 0;
        } else {
            u = v;
            DOTPROD(_value,u,u);
            even = 1;
        }  
    }
    
    void _wait () {}
};


#endif
