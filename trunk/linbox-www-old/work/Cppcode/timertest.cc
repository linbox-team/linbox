#include "givtimer.h"
#include <vector>

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





main(int ac,
     char **av)
{
  Timer T;

  T.start();

  const int qq=101;
  int x(2), y(9), z(5);

  GFp<int> F(qq), H(11);

  FieldTrait<GFp<int> >  KKF(F);
  FieldTrait<GFp<int> >  KKH(H);

  int i;
  int n    = atoi(av[1]);
  int iter = atoi(av[2]);

  vector<int> A(n);
  vector<int> B(n); 

  for (i=0; i< n; ++i) 
    {
      A[i] = i; B[i] = i+1;
    }

  T.stop();
  cout << "Time stl vector setup\n" << T << endl;
  T.start();

  // Inner product using generic field mul and add on stl vectors

  int dotprod =0, tmp;
  for (i=0; i<iter; ++i) {
    dotprod = 0;
    for (int j=0; j< n; ++j) 
      {
	KKF.mul(tmp, A[j], B[j]);
	KKF.add(dotprod, dotprod, tmp);
      }
  }

  T.stop();
  cout << "Time GFp\n" << T << endl;
  cout << dotprod << endl << endl;
  T.start();
  B[1] = 3;

  // Inner product using generic field mul and add on stl vectors
  for (i=0; i<iter/2; ++i) {
    dotprod = 0;
    for (int j=0; j< n; ++j) 
      {
	KKF.mul(tmp, A[j], B[j]);
	KKF.add(dotprod, dotprod, tmp);
      }
  }
  T.stop();
  cout << "Time half as many GFp\n" << T << endl;
  cout << dotprod << endl << endl;
  T.start();
  A[1] = 3;

  // Inner product using generic field mul and add on stl vectors
  for (i=0; i<iter; ++i) {
    dotprod = 0;
    for (int j=0; j< n; ++j) 
      {
	KKF.mul(tmp, A[j], B[j]);
	KKF.add(dotprod, dotprod, tmp);
      }
  }
  T.stop();
  cout << "Time full GFp again \n" << T << endl;
  cout << dotprod << endl << endl;
  
  // Inner product using direct mul and add on stl vectors

  T.start();
  for (i=0; i<iter; ++i) 
    {
      dotprod =0;
      for (int j=0; j<n; ++j) 
	{
	  dotprod = (dotprod + (A[j] * B[j]) % qq) % qq;
	}
    }
  T.stop();
  cout << "Time stl direct \n" << T << endl;
  cout << dotprod << endl;
  
  T.start();
  // Inner product using direct mul and add on C arrays.
  int* Ap = new int[n];
  int* Bp = new int[n];
  for (i=0; i< n; ++i) 
    {
      Ap[i] = i; Bp[i] = i+1;
    }


  T.stop();
  cout << "Time array setup \n" << T << endl;

  T.start();
  for (i=0; i<iter; ++i) 
    {
      dotprod =0;
      for (int j=0; j<n; ++j) 
	{
	  dotprod = (dotprod + (Ap[j] * Bp[j]) % qq) % qq;
	}
    }
  T.stop();
  cout << "Time arrays direct \n" << T << endl;
  cout << dotprod << endl;

}
