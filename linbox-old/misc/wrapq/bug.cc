class foo
{private:
  foo(const foo& b); // no copy constructors please.
 public:
  // no arg constructor does the mem allocation.
  virtual foo& operator=(const foo& b) =0; // assignment op.
};

class bar : public foo
{
  int _p;
  bar() _p(3) {}
  foo& operator=(const foo& b) 
  { _p = ((bar&)b)._p; return *this; }
};

main()
{
  bar A;
  bar B;
  B = A; 
  bar C = A;
}

warning: `class foo' only defines private constructors
and has no friends
cpcstor.cc:12: parse error in method specification before `('
cpcstor.cc:15: parse error before `('
cpcstor.cc:24: parse error at end of input
cpcstor.cc:24: Internal error #19990916.
cpcstor.cc:24: Internal compiler error in poplevel, at ../gcc/cp/decl.c
:1315
Please submit a full bug report.
See <URL:http://www.gnu.org/software/gcc/bugs.html> for instructions.                    
