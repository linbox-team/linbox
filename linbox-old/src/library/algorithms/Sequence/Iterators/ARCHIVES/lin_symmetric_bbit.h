// ================================================================
// LinBox Project 1999
// Black Box iterator and container 
// For symmetric matrix with same left and right vector
// the sequence is u^t v, u^t A v, ...,  u^t A^n v,  
// Time-stamp: <07 Mar 00 18:46:56 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __THUNK_DEFAULT_H__
#define __THUNK_DEFAULT_H__

#include <lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t>
class BB_Symmetric_Iterator : public Base_BB_Iterator< BlackBoxDomain> {
public:
        //-- Constructors
    BB_Symmetric_Iterator() {} 

    BB_Symmetric_Iterator(BlackBoxDomain_t * BD) 
            : Base_BB_Iterator<BlackBoxDomain>(BD) { init(); }
    BB_Symmetric_Iterator(BlackBoxDomain_t * BD, const Vecteur& u0) 
            : Base_BB_Iterator<BlackBoxDomain>(BD) { init(u0); }

private:

        /// Random Left vector, Right is the same
    Type_t& init() {
        even = 1;
        u.resize(_BB_domain->n_row());
        for(long i=u.size();i--;)
            _domain.random(u[i]);
        v.resize(_BB_domain->n_row());
        DOTPROD(_value, u, u);
    }        

        /// Same Left and Right user vectors 
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
        switch(even) {
            case 1:
                _BB_domain->Apply( v, u);  // GV
                DOTPROD(_value,u,v);       // GV 
                even = 2;
                break;
            case 2:
                DOTPROD(_value,v,v);       // GV
                even = 3;
                break;
            case 3:
                _BB_domain->Apply( u, v); 
                DOTPROD(_value,v,u);
                even = 0;
                break;
            default:
                DOTPROD(_value,u,u);
                even = 1;
                break;
        }
    }
    
    void _wait () {}
};

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t>
class BB_Symmetric_Container : public Base_BB_Container< BlackBoxDomain> {
public:
    BB_Symmetric_Container() {} 

    BB_Symmetric_Container(BlackBoxDomain_t * D) 
            : Base_BB_Container< BlackBoxDomain>(D) {}
    
    typedef BB_Symmetric_Iterator< BlackBoxDomain, Vecteur > const_iterator;
    const_iterator begin() const { return const_iterator(_BB_domain); }
    const_iterator end() const { return const_iterator(); }
};


#endif
