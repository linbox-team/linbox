#ifndef _LinBox_Map_H
#define _LinBox_Map_H
// ==========================================================================
// Authors: JGD
// Time-stamp: <08 Sep 08 14:23:02 Jean-Guillaume.Dumas@imag.fr> 
// ==========================================================================
#include "linbox/field/hom.h"

namespace LinBox {

template< class Source, class Target > 
struct Map {

    typedef Hom<Source,Target> Homomorphism;

    Map(const Homomorphism& H) : _hom(H) {}
    Map(const Source& S, const Target& T) : _hom(S,T) {}


    template< typename ContainerTargetElement, 
              typename ContainerSourceElement >
    ContainerTargetElement& operator() (
        ContainerTargetElement& tgt, 
        const ContainerSourceElement& src) {

        tgt.resize(src.size());
        typename ContainerTargetElement::iterator tgt_it(tgt.begin());
        typename ContainerSourceElement::const_iterator src_it(src.begin());
        
        for( ; src_it != src.end(); ++src_it, ++tgt_it)
            this->_hom.image (*tgt_it, *src_it);

        return tgt;
    }


private:
    Homomorphism _hom;
};


template< class Source, class Target > 
struct PreMap {

    typedef Hom<Source,Target> Homomorphism;

    PreMap(const Homomorphism& H) : _hom(H) {}
    PreMap(const Source& S, const Target& T) : _hom(S,T) {}

    template< typename ContainerSourceElement, 
              typename ContainerTargetElement >
    ContainerSourceElement& operator() (
        ContainerSourceElement& src, 
        ContainerTargetElement tgt) {

        src.resize(tgt.size());
        typename ContainerSourceElement::iterator src_it(src.begin());
        typename ContainerTargetElement::const_iterator tgt_it(tgt.begin());
        
        for( ; tgt_it != tgt.end(); ++src_it, ++tgt_it) {
            this->_hom.preimage (*src_it, *tgt_it);
        }

        return src;
    }


private:
    Homomorphism _hom;
};

 
}
#endif // _LinBox_Map_H

    
        
            
    
