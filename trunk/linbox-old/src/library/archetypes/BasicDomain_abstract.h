// -bds 7/00
#ifndef __BasicDomain_abstract__
#define __BasicDomain_abstract__
#include <iostream.h>
#include <BasicObject_abstract.h>
#include <lin_integer.h>

namespace linbox{
/*
BasicDomain_abstract is also known as Set_abstract.
The BasicDomain members are the BasicObject functionality plus:
class element, 
functions areEqual(), areNotEqual(), random(), cardinality().
*/

class BasicDomain_abstract : public BasicObject_abstract
{ protected:
    typedef BasicObject_abstract baseElement;
  public:
    typedef baseElement element;

    // should be able to use a.init(), I think.
    //virtual element& init( element& a ) const 
    //= 0;
    //virtual element& init( element& a, const element& b ) const 
    //= 0;
    //virtual element& assign ( element& r, const element& a ) const 
    //= 0;
    //copy now archaic?
    // should be able to use a.read(s), I think.
    //virtual istream& read( istream& s, element& a ) const 
    //= 0;
    //virtual ostream& write( ostream& s, const element& a ) const 
    //= 0;

    virtual bool areEqual( const element& a, const element& b ) const 
    = 0;
    virtual bool areNotEqual( const element& a, const element& b ) const 
    = 0;
    virtual const integer& cardinality(integer& n) const 
    = 0;
    // random() to be replaced by random iterator construct.
    virtual element& random ( element& r, const integer& n ) const 
    = 0;
}; // BasicDomain_abstract

/// constants for use with cardinality
enum{ INF = -1, UNDEFINED = -2, MINF = -3};

// Set is alias for BasicDomain
typedef BasicDomain_abstract Set_abstract;

} // namespace linbox
#endif

