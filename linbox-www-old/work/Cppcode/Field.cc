#include <iostream>





// This is a specialized field. Use it to build 
template<class ELEMENT>
class GFp {

private:

  ELEMENT characteristic;

public:

  typedef ELEMENT GFpELEMENT;

  GFp();
  GFp( ELEMENT& elt) { characteristic = elt; }
  GFp( const ELEMENT& elt) { characteristic = elt; }
  void mul (ELEMENT& a, ELEMENT& b, ELEMENT& c) { a = (b * c) % characteristic;}
  void add (ELEMENT& a, ELEMENT& b, ELEMENT& c) { a = (b + c) % characteristic;}
  void sub (ELEMENT& a, ELEMENT& b, ELEMENT& c) { a = (b - c+ characteristic) % characteristic;}
};


template <>
class FieldTrait< GFp > {
public:
  GFp *Fptr;

  FieldTrait( GFp&  Fp) {
    Fptr = &Fp;
  }
  typedef   GFp::GFpELEMENT  FTELEMENT;

  void mul (FTELEMENT& a, FTELEMENT& b, FTELEMENT& c) { a = (b * c) % characteristic;}
  void add (FTELEMENT& a, FTELEMENT& b, FTELEMENT& c) { a = (b + c) % characteristic;}
  void sub (FTELEMENT& a, FTELEMENT& b, FTELEMENT& c) { a = (b - c+ characteristic) % characteristic;}



  
};


main(int ac,
     char **av)
{
  int x(2), y(9), z(5);
  GFp<int> F(7), H(11);
  F.add( x, y, z);
  cout << x << endl;
 F.mul( x, y, z);
  cout << x << endl;
 H.sub( x, x, z);
  cout << x << endl;
}
