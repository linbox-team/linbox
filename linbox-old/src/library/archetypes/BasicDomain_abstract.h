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
    typedef BasicObject_abstract eltbase;
  public:
    typedef eltbase element;

    virtual bool areEqual( const eltbase& a, const eltbase& b ) const 
    = 0;
    virtual bool areNotEqual( const eltbase& a, const eltbase& b ) const 
    = 0;
    virtual const integer& cardinality(integer& n) const 
    = 0;
    // random() to be replaced by random iterator construct.
    virtual eltbase& random ( eltbase& r, const integer& n ) const 
    = 0;
}; // BasicDomain_abstract

/// constants for use with cardinality
enum{ INF = -1, UNDEFINED = -2, MINF = -3};

// Set is alias for BasicDomain
typedef BasicDomain_abstract Set_abstract;

} // namespace linbox
#endif

