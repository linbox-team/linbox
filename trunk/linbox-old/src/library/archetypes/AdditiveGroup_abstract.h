// -bds 7/00
#ifndef __AdditiveGroup_abstract__
#define __AdditiveGroup_abstract__
#include <iostream.h>
#include <BasicDomain_abstract.h>
#include <lin_integer.h>

namespace linbox{

class AdditiveGroup_abstract : public BasicDomain_abstract
{ protected:
  public:
    typedef BasicDomain_abstract::element element;
    virtual element& add(element& r, const element& a, const element& b)const = 0;
}; // AdditiveGroup_abstract

} //namespace linbox
#endif
