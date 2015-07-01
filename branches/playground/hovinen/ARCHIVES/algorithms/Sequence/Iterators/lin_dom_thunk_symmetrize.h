// ================================================================
// LinBox Project 1999
// Symmetrizing Thunk (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is u^t u, (A u)^t (A u) = u^t (A^t A) u, 
// (A^t (A u))^t (A^t (A u)) = u^t (A^t A)^2 u , etc.
// Time-stamp: <05 Jan 00 12:53:27 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __THUNK_SYMMETRIZE_H__
#define __THUNK_SYMMETRIZE_H__
#include <lin_dom_thunk_bbase.h>

template<class BlackBoxDomain>
class BBSymmetrizeThunk : public BBThunkDom< BlackBoxDomain > {
public:
        //-- Constructors
    BBSymmetrizeThunk() {} 

    BBSymmetrizeThunk(BlackBoxDomain_t * BD) 
            : BBThunkDom<BlackBoxDomain>(BD) { init(); }
    
    BBSymmetrizeThunk(BlackBoxDomain_t * BD, const Domain_t& D) 
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
            _BB_domain->Apply(v, u);
            DOTPROD(_value,v,v); 
            even = 0;
        } else {
            _BB_domain->ApplyTrans( u, v); 
            DOTPROD(_value,u,u);
            even = 1;
        }  
    }
    
    void _wait () {}
};


#endif
