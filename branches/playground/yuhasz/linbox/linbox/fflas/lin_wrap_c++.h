/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#ifndef _LINBOX_OperatorWrapper_Dom_H_
#define _LINBOX_OperatorWrapper_Dom_H_
#include <iostream>
#include "linbox/integer.h"
using namespace LinBox;
// ==========================================================================
// Time-stamp: <11 Apr 03 17:03:23 Jean-Guillaume.Dumas@imag.fr>
// (c) Givaro Team
// date: 1999
// version: 
// Description: Wraps class with +,*,-,\ into a domain
// ==========================================================================

  // ------------------------------------------------- class OperatorWrapper

template<class TT> class OperatorWrapper {
    typedef TT Rep;
public:
    typedef Rep element;
    typedef Rep Element;

        // ----- Representation of vector of the element
    typedef Rep* Array;
    typedef const Rep* constArray;
        // ----- Constantes 
    const Rep zero;
    const Rep one;

    size_t size() const { return 0; }
    size_t cardinality() const { return 0; }
    
    OperatorWrapper() : zero(0), one(1) {};
    ~OperatorWrapper() {};

        // Initialization of elements
    Rep& init( Rep& a ) const;
    Rep& init( Rep& r, const double& a ) const { r = (Rep)a; }
    Rep& init( Rep& r, const long int& a ) const { r = (Rep)a; }
    Rep& init( Rep& r, const int& a ) const { r = (Rep)a; }
    Rep& init( Rep& r, const unsigned long int& a ) const { r = (Rep)a; }
    Rep& init( Rep& r, const unsigned int& a ) const { r = (Rep)a; }
    Rep& init( Rep& r, const integer& a ) const { r = (Rep)a; }

        // Assignment of elements
    Rep& assign(Rep& r, const Rep&a) const;
    double& convert( double& r, const Rep& a ) const { return r = double(a);}
    int& convert( int& r, const Rep& a ) const { return r = int(a);}

        // Test operators
    inline int operator== (const OperatorWrapper<TT>& a) const ;
    inline int operator!= (const OperatorWrapper<TT>& a) const ;

        // Miscellaneous functions
    bool iszero( const Rep& ) const;
    bool isnzero( const Rep& ) const;
    bool isone ( const Rep& ) const;
    bool isnone ( const Rep& ) const;
    bool isequal( const Rep&, const Rep&) const;
    bool isnequal( const Rep&, const Rep&) const;
    bool islt( const Rep&, const Rep&) const;
    bool isgt( const Rep&, const Rep&) const;
    bool areEqual( const Rep a, const Rep b ) const { return isequal(a,b);  } ;
    bool isZero( const Rep a ) const { return iszero(a); }
    bool isOne ( const Rep a ) const { return isone(a); }

    
    size_t length ( const Rep& ) const;

        // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
    Rep& mul (Rep& r, const Rep&a, const Rep&b) const;
    Rep& div (Rep& r, const Rep&a, const Rep&b) const;
    Rep& add (Rep& r, const Rep&a, const Rep&b) const;
    Rep& sub (Rep& r, const Rep&a, const Rep&b) const;
    Rep& neg (Rep& r, const Rep&a) const;
    Rep& inv (Rep& r, const Rep&a) const;
    Rep& sqrt (Rep& r, const Rep&a) const;

    Rep& mulin (Rep& r, const Rep&a) const;
    Rep& divin (Rep& r, const Rep&a) const;
    Rep& addin (Rep& r, const Rep&a) const;
    Rep& subin (Rep& r, const Rep&a) const;
    Rep& negin (Rep& r) const;
    Rep& invin (Rep& r) const;
    Rep& sqrtin (Rep& r) const;

    Rep& axpy (Rep& r, const Rep&a, const Rep&b, const Rep&c) const;
    Rep& axmy (Rep& r, const Rep&a, const Rep&b, const Rep&c) const;
    Rep& axpyin (Rep& r, const Rep&a, const Rep&b) const;
    Rep& axmyin (Rep& r, const Rep&a, const Rep&b) const;


        // --- IO methods
    std::istream& read ( std::istream& s );
    std::ostream& write( std::ostream& s ) const;

    TT write (const Rep& a) const { return a; }
    Rep& read (Rep&, const TT ) const;

    std::istream& read ( std::istream& s, Rep& a ) const;
    std::ostream& write( std::ostream& s, const Rep& a ) const;

        // ----- random generators
    template<class RandIter> Rep& random(RandIter& g, Rep& r) const ;
    template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const ;
};


#include "lin_wrap_c++.inl"



#endif
