/*
wrapint.cc   
compare intdom (element type int), no inheritance
             vs 
Intdom (element type Int), where Int (derived from objectbase)
                 wraps int and Intdom is derived from domainbase.                    

Which is faster: wrapint 1 or wrapint 2 ?

----------  wrapint.cc --------------
*/
#include "base.h"
#include <iostream.h>

class Int : public objectbase
{private:
  int _n;
  friend class Intdom;
 public:

  void cstor(const objectbase& b)
  {((Int*) this)->_n = ((Int&)b)._n; };

  objectbase& operator=(const objectbase& b)
  {((Int*) this)->_n = ((Int&)b)._n; return *this; }

  Int& operator=(const int b)
  { _n = b; return *this; }

  friend ostream& operator<<(ostream& out, const Int& b)
  ;
};
  ostream& operator<<(ostream& out, const Int& b)
  { return out << b._n;}

class Intdom : public domainbase
{private: 
  int _p;  // domains can have their data
 public:

  typedef Int element;

  void cstor() 
  { _p = 0; }

  void cstor(const objectbase& b) 
  { ((Intdom*)this)->_p = ((Intdom&)b)._p; }

  objectbase& operator=(const objectbase& b)
  { ((Intdom*)this)->_p = ((Intdom&)b)._p; return *this; }

  domainbase::element& mul(domainbase::element& r, const domainbase::element& a, const domainbase::element& b) const
  { ((Intdom::element&)r)._n = ((Intdom::element&)a)._n + ((Intdom::element&)b)._n; 
    return r;
  }
};
   
class intdom //: public domainbase
{protected: 
  int _p;
 public:

  //typedef something_derived_from_domainbase element;
  typedef int element;

  void cstor() 
  {_p = 0; }

  void cstor(const objectbase& b) 
  { _p = ((intdom&)b)._p; }

  //objectbase& operator=(const objectbase& b)
  //{ _p = ((intdom&)b)._p; return *this; }

  // a fake just to not be an abstract class.
  domainbase::element& mul(domainbase::element& r, const domainbase::element& a, const domainbase::element& b) const
  { (intdom::element&)r = (intdom::element&)a + (intdom::element&)b; 
    return r;
  }

  // the real mul
  int& mul(int& r, const int& a, const int& b) const
  { r = a + b;
    return r;
  }
};

//////////  application /////////////////////
#include <iostream.h>

objectbase& genericalg(const domainbase& D, domainbase::element& b, const domainbase::element& a)
;

int random(){static int n = 1; return n = n- 1; }

int main(int ac, char* av[])
{
  if (av[1][0] == '1') // measure time on raw ints
  { intdom A;
    int a, b;
    a = 0;
    b = 1;
    //a = b - 1;
    cout << a << endl;
    for (int i = 0; i < 100000000; ++i)
    { A.mul(a, a, b);}  
    cout << a << endl;
  } 
  else  // measure time on wrapped ints
  { Intdom A;
    Intdom::element a, b;
    cout << sizeof(a) << endl;
    a = 0;
    b = 1;
     cout << a << endl;
    for (int i = 0; i < 100000000; ++i)
    { A.mul(a, a, b);}  
    //{ b = random(); A.mul(a, a, b);}  
    cout << a << endl;
  }
}                                         
