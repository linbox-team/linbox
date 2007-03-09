/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ======================================================================= //
// Time-stamp: <09 Mar 07 20:10:52 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_RATIONAL_CRA_H
#define __LINBOX_RATIONAL_CRA_H

#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/cra-early-single.h"

namespace LinBox {

    template<class Domain_Type>
    struct EarlySingleRatCRA : public EarlySingleCRA<Domain_Type> {
        typedef Domain_Type				Domain;
        typedef EarlySingleCRA<Domain> 			Father_t;
        typedef typename Father_t::DomainElement 	DomainElement;
        typedef EarlySingleRatCRA<Domain>		Self_t;
        PID_integer _ZZ;

        Integer					Numer0;
        Integer					Denom0;

    public:

	    EarlySingleRatCRA(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD) 
			    : Father_t(EARLY) {}
	    
	    void progress (const Domain& D, const DomainElement& e) {
		    DomainElement u0, m0;
		    
		    fieldreconstruct(residue_, D, e, D.init(u0,residue_), D.init(m0,primeProd_), Integer(residue_), primeProd_);
		    D.characteristic( nextm_ );
		    primeProd_ *= nextm_;
		    Integer a, b;
		    _ZZ.reconstructRational(a, b, residue_, primeProd_);
		    if ((a == Numer0) && (b == Denom0))
			    ++occurency;
		    else {
			    occurency = 1;
			    Numer0 = a;
			    Denom0 = b;
		    } 
	    }
	    
	    void initialize (const Domain& D, const DomainElement& e) {
		    Father_t::initialize(D, e);
		    _ZZ.reconstructRational(Numer0, Denom0, residue_, primeProd_);
	    }
	    
    };
}

#endif
