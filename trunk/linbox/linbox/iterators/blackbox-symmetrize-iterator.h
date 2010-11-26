/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/blackbox-symmetrize-iterator.h
 * Copyright (C) 1999, 2001, 2003 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_bbcontainer_symmetrize_H
#define __LINBOX_bbcontainer_symmetrize_H

// ================================================================
// Symmetrizing iterator (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is u^t u, (A u)^t (A u) = u^t (A^t A) u, 
// (A^t (A u))^t (A^t (A u)) = u^t (A^t A)^2 u , etc.
// ================================================================


#include <linbox/algorithms/blackbox-container-base.h>

namespace LinBox 
{

template<class Field, class Vector>
class BlackboxSymmetrizeIterator : public BlackboxContainerBase< Field, Vector > {
public:
    BlackboxSymmetrizeIterator() {} 

    BlackboxSymmetrizeIterator(Blackbox * D, const Field& F, const Vector& u0) 
            : BlackboxContainerBase< Field, Vector >(D, F) { init(u0); }
    
    BlackboxSymmetrizeIterator(Blackbox * D, const Field& F) 
            : BlackboxContainerBase< Field, Vector >(D, F) { init( Field::RandIter(_field) ); }
    
private:
    void _launch () {
        if (casenumber) {
            casenumber = 0;
            _BB_domain->Apply(v, u);
            DOTPROD(_value,v,v); 
        } else {
            casenumber = 1;
            _BB_domain->ApplyTrans( u, v); 
            DOTPROD(_value,u,u);
        }  
    }
    
    void _wait () {}

};
 
};


#endif // __LINBOX_bbcontainer_symmetrize_H

