/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#ifndef _LINBOX_DOM_TT_H_
#define _LINBOX_DOM_TT_H_

// ==========================================================================
// Time-stamp: <03 Jul 00 12:18:15 Jean-Guillaume.Dumas@imag.fr>
// (c) Givaro Team
// date: 1999
// version: 
// Description: Wraps class with +,*,-,\ into a domain
// ==========================================================================

  // ------------------------------------------------- class TTDom

template<class TT> class TTDom {
    typedef TT Rep;
public:
    typedef TT element;
    typedef element Element;

        // ----- Representation of vector of the element
    typedef Rep* Array;
    typedef const Rep* constArray;
        // ----- Constantes 
    const Rep zero;
    const Rep one;

    size_t size() const { return 0; }
    size_t cardinality() const { return 0; }
    
    TTDom() : zero(0), one(1) {};
    ~TTDom() {};

        // Initialization of elements
    Rep& init( Rep& a ) const;
    Rep& init( Rep& r, const double& a ) const { r = (Rep)a; };
    Rep& init( Rep& r, const long int& a ) const { r = (Rep)a; };
    Rep& init( Rep& r, const int& a ) const;
        // Assignment of elements
    Rep& assign(Rep& r, const Rep&a) const;
    double& convert( double& r, const Rep& a ) const { return r = double(a);}
    int& convert( int& r, const Rep& a ) const { return r = int(a);}

        // Test operators
    inline int operator== (const TTDom<TT>& a) const ;
    inline int operator!= (const TTDom<TT>& a) const ;

        // Miscellaneous functions
    bool iszero( const Rep& ) const;
    bool isnzero( const Rep& ) const;
    bool isone ( const Rep& ) const;
    bool isnone ( const Rep& ) const;
    bool isequal( const Rep&, const Rep&) const;
    bool isnequal( const Rep&, const Rep&) const;
    bool islt( const Rep&, const Rep&) const;
    bool isgt( const Rep&, const Rep&) const;
    
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
    istream& read ( istream& s );
    ostream& write( ostream& s ) const;

    TT write (const Rep& a) const { return a; }
    Rep& read (Rep&, const TT ) const;

    istream& read ( istream& s, Rep& a ) const;
    ostream& write( ostream& s, const Rep& a ) const;

        // ----- random generators
    template<class RandIter> Rep& random(RandIter& g, Rep& r) const ;
    template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const ;
};


#include "lin_wrap_c++.inl"



#endif
