/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ======================================================================= //
// Time-stamp: <09 Mar 07 20:20:44 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_RATIONAL_CRA_H
#define __LINBOX_RATIONAL_CRA_H

#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/cra-full-multip.h"

namespace LinBox {

#include "linbox/algorithms/cra-early-multip.h"

namespace LinBox {

    template<class Domain_Type>
    struct EarlyMultipRatCRA : public virtual EarlyMultipCRA<Domain-Type>, public virtual FullMultipRatCRA<Domain_Type> {
        typedef Domain_Type				Domain;
        typedef typename Father_t::DomainElement 	DomainElement;
        typedef EarlyMultipRatCRA<Domain>		Self_t;
   public:

 
        EarlyMultipRatCRA(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD) 
                : Father_t(EARLY) {}

        
        template<template<class T> class Vect>
        void initialize (const Domain& D, const Vect<DomainElement>& e) {
		    EarlyMultipCRA<Domain>::initialize(D, e);
        }

        template<template<class T> class Vect>
        Vect<Integer>& result(Vect<Integer>& num, Integer& den) {
            return FullMultipRatCRA<Domain>::result(num, den);
        }

    };
}

#endif
