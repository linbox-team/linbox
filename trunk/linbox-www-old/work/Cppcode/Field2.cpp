#include "givtimer.h"


// This is a specialized field. Use it to build a trait
template<class ELEMENT>
class GFp {

private:

  ELEMENT characteristic;

public:
  typedef ELEMENT GFpELEMENT;

  GFp();
  GFp( ELEMENT& elt) { characteristic = elt; }
  GFp( const ELEMENT& elt) { characteristic = elt; }

  inline void mul (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    { a = (b * c) % characteristic;}

  inline void add (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    { a = (b + c) % characteristic;}

  inline void sub (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    { a = (b - c+ characteristic) % characteristic;}
};



template < class KK >
class FieldTrait {  };

template < class ELEMENT>
class FieldTrait<  GFp<ELEMENT>  > : public GFp<ELEMENT> {
public:
  typedef GFp<ELEMENT> Field;
  //  Field *Fptr;

  FieldTrait( const ELEMENT& p ) : GFp<ELEMENT>(p){;}
  
  FieldTrait( Field& F) : GFp<ELEMENT>(F) {
    ;
  }
  /*
  void mul (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    {  mul(a,b,c); }
  
  
  void add (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    {  add(a,b,c); }
  
  void sub (ELEMENT& a, ELEMENT& b, ELEMENT& c) 
    {  sub(a,b,c); }
  */

};





main(int ac,
     char **av)
{
  Timer T;
  T.start();
  int x(2), y(9), z(5);
  GFp<int> F(7), H(11);
  FieldTrait<GFp<int> >  AbsField(F);
  FieldTrait<GFp<int> >  HAbsField(H);
  AbsField.add( x, y, z);
  cout << x << endl;

 AbsField.mul( x, y, z);
  cout << x << endl;

  int p(5), q(6);
  int i;
  int n = atoi(av[1]);
  int iter = atoi(av[2]);
  int* A = new int [n]; 
  int* B = new int [n]; 
  for (i=0; i<n; ++i) {
    A[i] = i; B[i] = i+1;
  }

  int dotprod =0, tmp;
  T.stop();
  cout << "Time setup\n" << T << endl;

  T.start();
  for (i=0; i<iter; ++i) {
    dotprod = 0;
    for (int j=0; j<n; ++j) {
      AbsField.mul(tmp, A[j], B[j]);
      AbsField.add(dotprod, dotprod, tmp);
    }
  }

  T.stop();
  cout << "Time AbsField \n" << T << endl;
  cout << dotprod << endl << endl;
  int chr= 7;
  T.start();
  for (i=0; i<iter; ++i) {
    dotprod =0;
    
    for (int j=0; j<n; ++j) {
      dotprod = (dotprod + (A[j] * B[j]) % chr) %chr;
    }
  }
  T.stop();
  cout << "Time C array direct \n" << T << endl;
  cout << dotprod << endl << endl;

}
