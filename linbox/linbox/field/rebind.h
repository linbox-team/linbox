/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/field/rebind.h
 * Copyright (C) 2005 Jean-Guillaume Dumas
 * Time-stamp: <10 Jun 05 11:52:32 Jean-Guillaume.Dumas@imag.fr> 
 */

template<class XXX, class U>
struct Rebind {
    typedef typename XXX::template rebind<U>::other other;
};


template<class T, template<class> class Container>
template<class U>
struct Rebind< Container<T>, U > {
    typedef Container<U> other;
};

    
