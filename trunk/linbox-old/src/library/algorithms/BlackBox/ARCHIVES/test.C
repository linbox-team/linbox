// ==========================================================================
// $Source$
// Copyright(c)'99 by LINBOX Team
// Authors: 
// $Id$
// ==========================================================================

#include <iostream.h>
#include "abstract_bb.h"

// --
// - Concrete class
// --
class My_BB : public BB_Base {
public:
  // - Concrete entry object for My_BB
  class Entry : public BB_Base::Entry {
  public:
    Entry( double d =0.0 ) : _d (d) {}

    // - make a copy of the object
    Entry* clone() const { return new Entry(_d); }

    // - delete a copy make by clone: (delete this ???)
    void killclone() { delete this; } 

    ostream& print( ostream& o ) const
    { return o << _d; }

    friend ostream& operator<<( ostream& o, const Entry& e ) 
    { return o << e._d; }

  protected:
    double _d; // a data
    friend class My_BB;
  };

  // - Ctor
  My_BB() {}

  // - Dstor
  ~My_BB() {}

  // - Make a copy of the My_BB object 
  My_BB* clone() const { return new My_BB(); }

  // - Delete a copy make by clone: (delete this ???)
  void killclone() { delete this; }

  // - Apply method: warning with the type casting:
  void apply( BB_Base::Entry*& r, const BB_Base::Entry* const& s ) 
  {
    if (s ==0) { cerr << "invalid argument" << endl; return; }
    if (r ==0) r = new Entry();
    // - here res is a My_BB::Entry and src also. 
    // They could differs.
    My_BB::Entry* res = (My_BB::Entry*) r; 
    const My_BB::Entry* src = (const My_BB::Entry*) s; 
  
    // - Direct access to internal state of My_BB::Entry objects
    res->_d = 10 * src->_d;

    // -
    cout << "In apply of My_BB: in: " << *src << ", out: " << *res << endl;
  } 
};


// -- An algorithm that compute  B o B where B is an abstract BB
void un_petit_algorithme( 
   BB_Abstract& B, 
   BB_Abstract::Entry& e_out, 
   const BB_Abstract::Entry& e_in 
)
{
  BB_Abstract::Entry tmp;
  B.apply( tmp, e_in );
  B.apply( e_out, tmp );
}

int main()
{
    BB_Abstract BB( new My_BB );
    BB_Abstract::Entry six( new My_BB::Entry(6) );
    BB_Abstract::Entry tentimestentimessix;
    un_petit_algorithme( BB, tentimestentimessix, six);
    cout << "in Main: Out: " << tentimestentimessix << endl;
} 
