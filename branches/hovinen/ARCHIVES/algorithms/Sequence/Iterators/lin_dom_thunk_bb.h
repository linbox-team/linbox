// ================================================================
// LinBox Project 1999
// Black Box Thunk
// Classical one, with the ssame left and right vector
// the sequence is u^t u, u^t A u, ...,  u^t A^n u,  
// Time-stamp: <05 Jan 00 12:53:33 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __THUNK_DEFAULT_H__
#define __THUNK_DEFAULT_H__

#include <lin_dom_thunk_bbase.h>

template<class BlackBoxDomain>
class BBThunk : public BBThunkDom< BlackBoxDomain > {
public:
        //-- Constructors
    BBThunk() {} 

    BBThunk(BlackBoxDomain_t * BD) 
            : BBThunkDom<BlackBoxDomain>(BD) { init(); }
    
    BBThunk(BlackBoxDomain_t * BD, const Domain_t& D) 
            : BBThunkDom<BlackBoxDomain>(BD,D) { init(); }

    Type_t& init() {
        even = 1;
        u.resize(_BB_domain->n_col());
        for(long i=u.size();i--;)
            _domain.random(u[i]);
        w = u;
        v.resize(_BB_domain->n_row());
        DOTPROD(_value, u, u);
    }        

    Type_t& init(const Vecteur& uu) {
        even = 1;
        u = uu;
        w = uu;
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


#endif
