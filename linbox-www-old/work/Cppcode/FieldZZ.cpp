#include "ZZ.h"
#include "givtimer.h"
#include <iostream.h>


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

  void mul (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    { a = (b * c) % characteristic;}

  void add (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    { a = (b + c) % characteristic;}

  void sub (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    { a = (b - c+ characteristic) % characteristic;}
};



template < class KK >
class FieldTrait {  };

template < class ELEMENT>
class FieldTrait<  GFp<ELEMENT>  >  {
public:
  typedef GFp<ELEMENT> Field;

    Field *Fptr;

      FieldTrait( Field&  Fp) 
      { Fptr = &Fp; }
      
      void mul (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
      {  Fptr->mul(a,b,c); }
      
      
      void add (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
      {  Fptr->add(a,b,c); }
      
      void sub (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
      {  Fptr->sub(a,b,c); }


};




typedef ZZ ETYPE;
main(int ac,
     char **av)
{
  Timer T;
  T.start();
  ETYPE x(2), y(9), z(5);
  GFp<ETYPE> F(7), H(11);
  FieldTrait<GFp<ETYPE> >  AbsField(F);
  FieldTrait<GFp<ETYPE> >  HAbsField(H);
  AbsField.add( x, y, z);
  cout << x << endl;

 AbsField.mul( x, y, z);
  cout << x << endl;

  ETYPE p(5), q(6);
  int i;
  int n = atoi(av[1]);
  int iter = atoi(av[2]);
  ETYPE* A = new ETYPE [n]; 
  ETYPE* B = new ETYPE [n]; 
  for (i=0; i<n; ++i) {
    A[i] << i; B[i] << i;
  }
    ETYPE dotprod, tmp;
    dotprod << 0;
    ETYPE    qq ; qq << 7;
  T.stop();
  cout << "Time for setup\n" << T << endl;

  T.start();
  for (i=0; i<iter; ++i) {
    dotprod << 0;
    for (int j=0; j<n; ++j) {
      //          rem(tmp, A[j] * B[j], qq);
      //          rem(dotprod, tmp+dotprod, qq);
          dotprod = (dotprod + (A[j] * B[j]) % qq) %qq;
    }
  }

  T.stop();
  cout << "Time for some kind of bogus case \n" << T << endl;
 cout << dotprod << endl << endl;
  
  T.start();

  for (i=0; i<iter; ++i) {
    dotprod<< 0;
    for (int j=0; j<n; ++j) {
      AbsField.mul(tmp, A[j], B[j]);
      AbsField.add(dotprod, dotprod, tmp);
    }
  }
  T.stop();
  cout << "Time for NTL ZZ used in some bogus way\n" << T << endl;
  cout << dotprod << endl << endl;

}



