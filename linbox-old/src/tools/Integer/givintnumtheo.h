// =================================================================== //
// Givaro : Euler's phi function
//          Primitive roots.
//          RSA scheme.
// Time-stamp: <31 Aug 00 20:14:20 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_NUMTHEORY_
#define _GIVARO_NUMTHEORY_

#include <iostream.h>
#include "givinteger.h"
#include "givintprime.h"
#include "givintfactor.h"

// =================================================================== //
// Givaro : Theoreme Chinois des Restes
// =================================================================== //
// const Integer tcr(Container<Integer>& Lr, Container<Integer>& Lm);
// const Integer& tcr(Integer& res, Integer& pm, Container<Integer>& Lr, Container<Integer>& Lm);

template<class RandIter = Random>
class IntNumTheoDom : public IntFactorDom<RandIter> {
public:
    IntNumTheoDom(RandIter& g = *(new RandIter())) 
            :  IntFactorDom<RandIter>(g) {}
// =================================================================== //
// Euler's phi function
// =================================================================== //
    Rep& phi(Rep& r, const Rep& n) const ;
    template< template<class> class Container> Rep& phi(Rep& res, const Container<Rep>& Lf, const Rep& n) const ;
// =================================================================== //
// Primitive Root
// =================================================================== //
    Rep& prim_root(Rep&, const Rep&) const ;
    Rep& prim_root(Rep&, int&, const Rep&) const ;
    Rep& lowest_prim_root(Rep&, const Rep&) const ;
    int is_prim_root(const Rep&, const Rep&) const ;
    Rep& order(Rep&, const Rep&, const Rep&) const ;
    int isorder(const Rep&, const Rep&, const Rep&) const ;

    template< template<class> class Container> short mobius(const Container<unsigned long>& lpow) const ;
    short mobius(const Rep& a) const;
};


#include "givintnumtheo.inl"

#endif
