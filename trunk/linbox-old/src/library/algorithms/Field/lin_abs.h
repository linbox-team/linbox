#ifndef _LINBOX_DOM_ABS_H_
#define _LINBOX_DOM_ABS_H_

// ==========================================================================
// (c) Givaro Team
// date: 1999
// version: 
// Description: Absolute value in an ordered domain
// Time-stamp: <05 Apr 00 12:54:45 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

#ifndef DOMABS
#define DOMABS(a) ( islt((a), zero) ? Domain::negin(a) : (a) )
#endif

// ------------------------------------------------- class AbsDom

template<class Domain> class AbsDom : public Domain {
public:
    typedef typename Domain::Rep Rep;

        // ----- Representation of vector of the element
    typedef Rep* Array;
    typedef const Rep* constArray;
        // ----- Constantes 

    AbsDom() : Domain( ) {};
    AbsDom(Domain& D) : Domain(D) { };
    ~AbsDom() {};

    Rep& init(Rep& a) const { Domain::init(a); return DOMABS(a); }  
    Rep& init(Rep& r, const Rep& a) const { Domain::init(r,a); return DOMABS(r); }  
    Rep& read(Rep& a, const int i) const { Domain::read(a,i); return DOMABS(a); }  
    Rep& read(Rep& a, const unsigned int i) const { Domain::read(a,i); return DOMABS(a); }  
    Rep& read(Rep& a, const long i) const { Domain::read(a,i); return DOMABS(a); }  
    Rep& read(Rep& a, const unsigned long i) const { Domain::read(a,i); return DOMABS(a); }  

    Rep& assign(Rep& r, const Rep a) const { Domain::assign(r,a); return DOMABS(r); }
    istream& read ( istream& s, Rep& a ) const { Domain::read(s,a); DOMABS(a); return s; }

    Rep& sub (Rep& r, const Rep&a, const Rep&b) const { Domain::sub(r,a,b); return DOMABS(r); }
    Rep& subin (Rep& r, const Rep&a) const { Domain::subin(r,a); return DOMABS(r); }
    Rep& neg (Rep& r, const Rep&a) const { return Domain::assign(r,a) ; }
    Rep& negin (Rep& r) const { return r; }

    Rep& axmy (Rep& r, const Rep&a, const Rep&b, const Rep&c) const { Domain::axmy(r,a,b,c); return DOMABS(r); }
    Rep& amxyin (Rep& r, const Rep&a, const Rep&b) const { Domain::axmyin(r,a,b); return DOMABS(r); }

    Rep& max(Rep& a, const Rep& b, const Rep& c) { return Domain::assign(a, islt(b,c) ? c : b); }
    Rep& maxin(Rep& a, const Rep& b) { if ( islt(a,b) ) Domain::assign(a, b); return a; }
    Rep& min(Rep& a, const Rep& b, const Rep& c) { return Domain::assign(a, islt(b,c) ? b : c); }
    Rep& minin(Rep& a, const Rep& b) { if ( islt(b,a) ) Domain::assign(a, b); return a; }
    
        // ----- random generators
    template<class RandIter>
    Rep& nonzerorandom(RandIter& g, Rep& a) { Domain::nonzerorandom(g,a); return DOMABS(a); }
    
        // ----- random generators
    template<class RandIter>
    Rep& random(RandIter& g, Rep& a) { Domain::random(g,a); return DOMABS(a); }
    
          
};


#endif _LINBOX_DOM_ABS_H_
