/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ======================================================================= //
// Time-stamp: <12 Mar 07 19:34:19 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_RATIONAL_EARLY_SINGLE_CRA_H
#define __LINBOX_RATIONAL_EARLY_SINGLE_CRA_H

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
		    
		    fieldreconstruct(this->residue_, D, e, D.init(u0,this->residue_), D.init(m0,this->primeProd_), Integer(this->residue_), this->primeProd_);
		    D.characteristic( this->nextM_ );
		    this->primeProd_ *= this->nextM_;
		    Integer a, b;
		    _ZZ.reconstructRational(a, b, this->residue_, this->primeProd_);
		    if ((a == Numer0) && (b == Denom0))
			    ++this->occurency_;
		    else {
			    this->occurency_ = 1;
			    Numer0 = a;
			    Denom0 = b;
		    } 
	    }
	    
	    void initialize (const Domain& D, const DomainElement& e) {
		    Father_t::initialize(D, e);
		    _ZZ.reconstructRational(Numer0, Denom0, this->residue_, this->primeProd_);
	    }
	    
    protected:
	    
        Integer& fieldreconstruct(Integer& res, const Domain& D1, const DomainElement& u1, DomainElement& u0, DomainElement& m0, const Integer& r0, const Integer& P0) {
                // u0 and m0 are modified
            D1.negin(u0);   	// u0 <-- -u0
            D1.addin(u0,u1);	// u0 <-- u1-u0
            D1.invin(m0);	// m0 <-- m0^{-1} mod m1
            D1.mulin(u0, m0);   // u0 <-- (u1-u0)( m0^{-1} mod m1 )
            D1.convert(res, u0);// res <-- (u1-u0)( m0^{-1} mod m1 )         and res <  m1
            res *= P0;      	// res <-- (u1-u0)( m0^{-1} mod m1 ) m0      and res <= (m0m1-m0)
            return res += r0;   // res <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0 and res <  m0m1
        }
    };
}

#endif