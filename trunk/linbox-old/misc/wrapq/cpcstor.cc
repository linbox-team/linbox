class foo
{private:
  foo(const foo& b); // no copy constructors please.
 public:
  foo(){}; // no arg constructor does the mem allocation, and possibly initialization.
  virtual foo& operator=(const foo& b) =0; // assignment op.
};

class bar : public foo
{//protected:
 public:
  int _p;
 public:
  bar() {_p = 7; }
 private:
  bar(const bar& b){}
 public:
  foo& operator=(const foo& b) 
  { _p = ((const bar&)b)._p; return *this; }
};

#include <iostream.h>
main()
{
  bar A;
  cout << A._p << endl;
  bar B;
  B = A; 
  cout << B._p << endl;
//  bar C = A;  can't use this syntax when copy const is private
//  cout << C._p << endl;
}
