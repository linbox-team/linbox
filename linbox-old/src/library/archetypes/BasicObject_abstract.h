// -bds 7/00
#ifndef __BasicObject_abstract__
#define __BasicObject_abstract__
#include <iostream.h>

namespace linbox{
class BasicObject_abstract
{ protected:
    typedef BasicObject_abstract self;
  public:
    // virtual self() = 0; i.e. you must have a public default constructor
    // virtual self(const self&) = 0; i.e. you must have a public copy constructor
    virtual self& init() = 0; // default cstor
    virtual self& init(const self& b) = 0; // copy cstor
    //virtual ~BasicObject_abstract() = 0;
    virtual self& operator=(const self& b) = 0;
    virtual istream& read(istream& instr) = 0;
    virtual ostream& write(ostream& outstr) const = 0;
    // copy, isequal are now archaic?
}; // BasicObject_abstract

template<class BO> 
// On class BO we assume 
// default constructor, copy constructor, destructor, op=, op<<, op>>.
class BasicObject_envelope : public BasicObject_abstract
{ protected: 
    typedef BasicObject_envelope<BO> myself;
  public:
    BasicObject_envelope<BO>(BO rep)
    : _rep(rep) { }
    myself& operator=(const myself& b)
    { _rep = b._rep;  return *this; }
    BO _rep;

    BasicObject_envelope<BO>()
    : _rep() { }
    BasicObject_envelope<BO>(const self& b) 
    { init(b); }
    virtual self& init() // a kind of clone
    { /* need a call to _rep constructor here */ 
      /////////return *this; } // default cstor
      return static_cast<self&> (*new myself); 
    }
    virtual self& init(const self& b)
    { _rep = static_cast<const myself&>(b)._rep; return *this; } // copy cstor
    virtual ~BasicObject_envelope()
    { _rep.~BO(); }
    virtual self& operator=(const self& b)
    { return *this = static_cast<const myself&>(b); }
    virtual istream& read(istream& instr)
    { return instr >> _rep; }
    virtual ostream& write(ostream& outstr) const 
    { return outstr << _rep; }

}; // BasicObject_envelope

} // namespace linbox
#endif
