#ifndef _LINBOX_DOM_TT_H_
#define _LINBOX_DOM_TT_H_

// ==========================================================================
// Time-stamp: <03 Jul 00 12:18:15 Jean-Guillaume.Dumas@imag.fr>
// (c) Givaro Team
// date: 1999
// version: 
// Description: Wraps class with +,*,-,\ into a domain
// ==========================================================================

  // ------------------------------------------------- class OperatorWrapper

template<class TT> class OperatorWrapper {
    typedef TT Rep;
public:
    typedef TT element;
    typedef OperatorWrapper<TT> Self_t;

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
    template<typename XXX> Rep& init( Rep& r, const XXX& a ) const { return r = TT(a); }
        // Assignment of elements
    Rep& assign(Rep& r, const Rep&a) const;

        // Test operators
    inline int operator== (const Self_t& a) const ;
    inline int operator!= (const Self_t& a) const ;

        // Miscellaneous functions
    short iszero( const Rep& ) const;
    short isZero( const Rep& a ) const { return iszero(a); }
    short isnzero( const Rep& ) const;
    short isone ( const Rep& ) const;
    short isnone ( const Rep& ) const;
    short isequal( const Rep&, const Rep&) const;
    short isnequal( const Rep&, const Rep&) const;
    short islt( const Rep&, const Rep&) const;
    short isgt( const Rep&, const Rep&) const;
    
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

    template <class XXX> XXX& convert (XXX& x, const Rep& a) const { return x=XXX(a); }
    TT write (const Rep& a) const { return a; }
    Rep& read (Rep&, const TT ) const;

    istream& read ( istream& s, Rep& a ) const;
    ostream& write( ostream& s, const Rep& a ) const;

        // ----- random generators
    template<class RandIter> Rep& random(RandIter& g, Rep& r) const ;
    template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const ;
};


#include "LinBox/lin_wrap_c++.inl"



#endif _LINBOX_DOM_TT_H_
