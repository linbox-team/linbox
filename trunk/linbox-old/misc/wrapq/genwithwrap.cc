/*
Which is faster:

/usr/proc/bin/ptime wrapint 1
or
/usr/proc/bin/ptime wrapint 2
?

----------  wrapint.cc --------------
*/
#include "base.h"
#include <iostream.h>

class Int : public objectbase
{private:
  int _n;
  friend class Intdom;
 public:
  void cstor(){};
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
  { ((Intdom::element&)r)._n = ((Intdom::element&)a)._n * ((Intdom::element&)b)._n; 
    return r;
  }
};
   
class Dbl : public objectbase
{//protected:
 public:
  double _n;
 public:
  void cstor(){};
  void cstor(const objectbase& b)
  {((Dbl*) this)->_n = ((Dbl&)b)._n; };
  objectbase& operator=(const objectbase& b)
  {((Dbl*) this)->_n = ((Dbl&)b)._n; return *this; }
  Dbl& operator=(const double b)
  { _n = b; return *this; }
  friend ostream& operator<<(ostream& out, const Dbl& b)
  ;
};
  ostream& operator<<(ostream& out, const Dbl& b)
  { return out << b._n;}

class Dbldom : public domainbase
{private:
  int _p;
 public:
  typedef Dbl element;
  void cstor() 
  { _p = 0; }
  void cstor(const objectbase& b) 
  { ((Dbldom*)this)->_p = ((Dbldom&)b)._p; }
  objectbase& operator=(const objectbase& b)
  { ((Dbldom*)this)->_p = ((Dbldom&)b)._p; return *this; }
  domainbase::element& mul(domainbase::element& r, const domainbase::element& a, const domainbase::element& b) const
  { ((Dbldom::element&)r)._n = ((Dbldom::element&)a)._n * ((Dbldom::element&)b)._n; 
    return r;
  }
};
   
class intdom : public domainbase
{protected: 
  int _p;
 public:
  typedef int element;
  void cstor() 
  {_p = 0; }
  void cstor(const objectbase& b) 
  { _p = ((intdom&)b)._p; }
  objectbase& operator=(const objectbase& b)
  { _p = ((intdom&)b)._p; return *this; }

  domainbase::element& mul(domainbase::element& r, const domainbase::element& a, const domainbase::element& b) const
  { return r; } // a fake
  int& mul(int& r, const int& a, const int& b) const
  { r = a * b;
    return r;
  }
};

//////////  application /////////////////////
#include <iostream.h>

objectbase& genericalg(const domainbase& D, domainbase::element& b, const domainbase::element& a)
;

int random(){static int n = 2; return n = n*2; }
double rand(){static double n = 2; return n = n*2; }

int main(int ac, char* av[])
{
  if (av[1][0] == '0') // illustrate separately compiled generic fn - no bloat.
  {
  Intdom A;
  Intdom::element a, b, c;
  a = random(); b = random(); c = random(); 
  cout << a << endl;
  cout << b << endl;
  cout << c << endl;
  c = genericalg(A, a, b);
  cout << a << endl;
  cout << b << endl;
  cout << c << endl;
  cout << endl;
  Dbldom B;
  Dbldom::element x, y, z;
  x = rand(); y = rand(); z = rand(); 
  cout << x << endl;
  cout << y << endl;
  cout << z << endl;
  z = genericalg(B, x, y);
  cout << x << endl;
  cout << y << endl;
  cout << z << endl;
  cout << endl;
  z = genericalg(B, a, y);
  cout << a << endl;
  cout << z << endl;
  cout << endl;
  a = genericalg(B, x, y);
  cout << a << endl;
  cout << x << endl;
  cout << endl;
  a = genericalg(A, x, y);
  cout << a << endl;
  cout << x << endl;
  } else
  if (av[1][0] == '1') // measure time on raw ints
  { intdom A;
    intdom::element a, b, c;
    for (int i = 0; i < 100000000; ++i)
    { b = random(); c = random(); A.mul(a, b, c);}  
  } 
  else  // measure time on wrapped ints
  { Intdom A;
    Intdom::element a, b, c;
    for (int i = 0; i < 100000000; ++i)
    { b = random(); c = random(); A.mul(a, b, c);}  
  }
}                                         
