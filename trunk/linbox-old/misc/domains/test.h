/* CR1_test.h   file to test some of your rings  -bds 3/00 */
/* --------------------------------------------------*/

#include <fstream.h>

template<class Set>
int testSet(Set S, int reps, Integer rb)
// functions tested are areEqual(), areEq(), random(), cardinality().
// SHOULD test for reflexivity, symmetry, transitivity of equality ops.
{ 
  Integers Z;
  typedef typename Set::element elt;

  Integer n; S.cardinality(n);
  //  cerr << ".  Cardinality is "; Z.write(cerr, n); cerr << endl;

  int e = 0; // error count

  Set R = S; // copy
  Integer m;  R.cardinality(m);
  if (! Z.areEqual(m,n))
    cerr << ++e << ". cardinality or set-copy bug\n"; 

  elt a, b, c;
  for (int i = 0; i < reps; i++)
  {
    R.random(a, rb);

    if (! R.areEq(a, a)) 
      cerr << ++e << ". areEq(a,a) is false\n"; 
    b = a;
// //R.write(cout, a); cout << " "; R.write(cout, b); cout << endl;
    if (R.areEq(a, b)) 
      cerr << ++e << ". areEq(a,b) is wrongly true\n"; 
    if (! R.areEqual(a, b)) 
      cerr << ++e << ". areEqual or element-assign bug\n"; 
  
  // open a file, write to it, read back, check equality
    ofstream os; os.open("junk");
    R.write(os, a);
    os.flush(); os.close(); 
    ifstream is; is.open("junk");
    R.read(is, b);
    if (! R.areEqual(a, b)) 
      cerr << ++e << ". areEqual or read/write bug\n"; 
    is.close();
  }

  // random check  
  /* What would be a good check (but assuming the underlying
      pseudo-random num gen is good)? */
  R.random(a, rb); R.random(b, rb); Integer p = rb; 
  while (R.areEqual(a, b) && p < 1000) R.random(b,rb);
  if (R.areEqual(a,b) && p >= 1000)
  { cerr << ++e << ". likely random bug\n"; 
    cerr << "some random values: ";
    for (int i =0; i < 20; i++) 
    { R.random(a, rb); R.write(cerr, a); cerr << " "; }
    cerr << endl;
  }

  // no generic check for second and third forms of random()
  // or for class element, for constructors,  

  if ( e > 0 ) cerr << "  Set test FAIL\n";
  else         cerr << "  Set test OK\n";
  return e;
}

template<class AdditiveGroup>
int testAdditiveGroup(AdditiveGroup G, int reps, Integer rb)
{ Integers Z;
  typedef typename AdditiveGroup::element elt;

  Integer n; G.cardinality(n);

  int e = testSet(G, reps, rb); // error count
  if (e > 0) reps = 1;

// functions tested are 
// zero(), isZero(), add(), neg(), sub(), addin(), negin(), subin(),
// Zprod, characteristic

  elt a, b, c, r, s, t;

  a = G.zero();
  if (! G.isZero(a))  
    cerr << ++e << ". zero() or isZero() bug\n";

  int i;
  for (i = 0; i < 10; i++) { G.random(a,rb); if (! G.isZero(a) ) break; }
  if ( i == 10 ) cerr << ++e << ". likely isZero or random bug\n";
  for (i = 0; i < 10; i++) { G.random(b,rb); if (! G.isZero(b) ) break; }
  if ( i == 10 ) cerr << ++e << ". likely isZero or random bug\n";
  for (i = 0; i < 10; i++) { G.random(c,rb); if (! G.isZero(c) ) break; }
  if ( i == 10 ) cerr << ++e << ". likely isZero or random bug\n";
  // commutativity
  G.add(r, a, b);
  G.add(s, b, a);
  G.sub(t, r, s);
  if (! G.isZero(t)) cerr << ++e << ". commutativity bug\n";
  // neg and sub
  G.sub(t, a, a);
  if (! G.isZero(t)) cerr << ++e << ". sub bug\n";
  G.neg(r, a);
  G.add(t, a, r);
  if (! G.isZero(t)) cerr << ++e << ". neg bug\n";
  G.sub(s, b, c);
  G.neg(r, c); G.add(t, b, r);
  if (! G.areEqual(s, t)) cerr << ++e << ". sub vs add-neg bug\n";
  // associativity
  G.add(r, a, b); G.add(s, r, c);
  G.add(r, b, c); G.add(t, a, r);
  if (! G.areEqual(s, t)) cerr << ++e << ". assoc bug\n";
  
  // inplace forms
  G.add(s, a, b); 
  t = a; G.addin(t, b); 
  if (! G.areEqual(s, t)) cerr << ++e << ". addin bug\n";
  G.sub(s, a, b); 
  t = a; G.subin(t, b); 
  if (! G.areEqual(s, t)) cerr << ++e << ". subin bug\n";
  G.neg(s, a); 
  t = a; G.negin(t); 
  if (! G.areEqual(s, t)) cerr << ++e << ". negin bug\n";

  Integer m, j;
  m = reps*(reps+1)/2;
  G.Zprod(r, m, a);
  s = G.zero();
  for(i = 1; i <= reps; i++)
  { j = i;
    G.Zprod(t, j, a);
    G.addin(s, t);
  }
  if (! G.areEqual(r, s)) cerr << ++e << ". Zprod/add bug\n";

  if ( e > 0 ) cerr << "  AdditiveGroup test FAIL\n";
  else         cerr << "  AdditiveGroup test OK\n";
  return e;
}

template<class Ring>
int testRing(Ring R, int reps, Integer rb)
{ Integers Z;
  typedef typename Ring::element elt;

  Integer n; R.cardinality(n);

  int e = testAdditiveGroup(R, reps, rb); // error count
  if (e > 0) reps = 1;

// functions tested are 
// mul(), mulin(), axpy(), axpyinx(), axpyiny()

  elt a, b, c, r, s, t;
  for (int i = 0; i < reps; i++)
  {
    R.random(a,rb); R.random(b,rb); R.random(c,rb);
  // mul assoc
    R.mul(r, a, b); R.mul(s, r, c);
    R.mul(r, b, c); R.mul(t, a, r);
    if (! R.areEqual(s, t)) cerr << ++e << ". mul assoc bug\n";

  // mulin, etc valid
    R.mul(s, a, b); 
    t = a; R.mulin(t, b);
    if (! R.areEqual(s, t)) cerr << ++e << ". mulin bug\n";
 
    R.axpy(s, a, b, c); 
    R.mul(r, a, b); R.add(t, r, c);
    if (! R.areEqual(s, t)) cerr << ++e << ". axpy bug\n";
    t = b; R.axpyinx(a, t, c);
    if (! R.areEqual(s, t)) cerr << ++e << ". axpyinx bug\n";
    t = c; R.axpyiny(a, b, t);
    if (! R.areEqual(s, t)) cerr << ++e << ". axpyiny bug\n";
 
    R.axmy(s, a, b, c); 
    R.mul(r, a, b); R.sub(t, r, c);
    if (! R.areEqual(s, t)) cerr << ++e << ". axmy bug\n";
    t = b; R.axmyinx(a, t, c);
    if (! R.areEqual(s, t)) cerr << ++e << ". axmyinx bug\n";
    t = c; R.axmyiny(a, b, t);
    if (! R.areEqual(s, t)) cerr << ++e << ". axmyiny bug\n";
  }

  if ( e > 0 ) cerr << "  Ring test FAIL\n";
  else         cerr << "  Ring test OK\n";
  return e;
}

template<class Ring1>
int testRing1(Ring1 R, int reps, Integer rb)
{ Integers Z;
  typedef typename Ring1::element elt;

  Integer n; R.cardinality(n);

  int e = testRing(R, reps, rb); // error count
  if (e > 0) reps = 1;

// functions tested are 
// one(), isOne()

  elt a, b, c, r, s, t;
  for (int i = 0; i < reps; i++)
  {
    R.random(a,rb); R.random(b,rb); 
  // oneness
    R.mul(r, a, R.one()); 
    if (! R.areEqual(r, a)) cerr << ++e << ". one bug\n";
    R.mul(r, R.one(), b); 
    if (! R.areEqual(r, b)) cerr << ++e << ". one bug\n";
    if (! R.isOne(R.one())) cerr << ++e << ". isOne bug\n";
  }
  R.random(a,rb); R.random(b,rb); 
  if ( R.isOne(a) && R.isOne(b)) cerr << ++e << ". likely isOne bug\n";

  // 1 - a^n = (1-a) * sum^{n-1) a^i
  R.sub(r, R.one(), a); // r = 1 - a
  s = R.one();
  t = R.zero();
  for(int i = 0; i < reps; i++) { R.addin(t, s); R.mulin(s, a); }
  // t = sum^{n-1} a^i, s = a^n  (n is reps)
  R.subin(s, R.one()); // s = 1 - a^n
  R.mulin(t, r);
  if (! R.areEqual(s, t) ) cerr << ++e << ". sum-powers bug\n";

  if ( e > 0 ) cerr << "  Ring1 test FAIL\n";
  else         cerr << "  Ring1 test OK\n";
  return e;
}

template<class CR1>
int testCR1(CR1 R, int reps, Integer rb)
{ Integers Z;
  typedef typename CR1::element elt;

  Integer n; R.cardinality(n);

  int e = testRing1(R, reps, rb); // error count
  if (e > 0) reps = 1;

// functions tested are 
// no new functions, just mul-comm

  elt a, b, c, r, s, t;
  for (int i = 0; i < reps; i++)
  {
    R.random(a,rb); R.random(b,rb); 
    // mul comm
    R.mul(r, a, b); R.mul(s, b, a);
    if (! R.areEqual(r, s)) cerr << ++e << ". mul comm bug\n";
  }

  // 1 - (ab)^n = (1-a) * sum^{n-1) (ab)^i
  s = R.one();
  t = R.zero();
  for(int i = 0; i < reps; i++) 
  { R.addin(t, s); 
    R.mul(r, b, s);
    R.mul(s, r, a); }
  // t = sum^{n-1} b^i a^i, s = b^n a^n  (n is reps)
  R.subin(s, R.one()); // s = b^n a^n - 1
  R.mul(r, a, b);
  R.subin(r, R.one()); // r = ab - 1
  R.mulin(t, r);
  if (! R.areEqual(s, t) ) cerr << ++e << ". sum-powers bug\n";

  if ( e > 0 ) cerr << "  CR1 test FAIL\n";
  else         cerr << "  CR1 test OK\n";
  return e;
}

template<class Field>
int testField(Field R, int reps, Integer rb)
{ Integers Z;
  typedef typename Field::element elt;

  Integer n; R.cardinality(n);

  int e = testCR1(R, reps, rb); // error count
  if (e > 0) reps = 1;

// functions tested are 
// div(), inv(), divin(), invin().

  elt a, b, c, r, s, t;
  for (int i = 0; i < reps; i++)
  {
    R.random(a,rb); R.random(b,rb); 

  // inv prop
    if (! R.isZero(a)) 
    { R.inv(r, a); R.mul(s, r, a);
      if (! R.isOne(s)) cerr << ++e << ". " << a << " inv bug" << b << "\n";

      R.inv(s, a); 
      t = a; R.invin(t);
      if (! R.areEqual(s, t)) cerr << ++e << ". invin bug\n";
    }

  // div, negin, divin according to def

    if (! R.isZero(b)) 
    {
      R.div(r, a, b); 
      R.inv(s, b); R.mul(t, a, s);
      if (! R.areEqual(r, t)) cerr << ++e << ". div bug\n";

      R.div(s, a, b); 
      t = a; R.divin(t, b);
      if (! R.areEqual(s, t)) cerr << ++e << ". divin bug\n";
    }
  }

  // 1 - (ab)^n/(1-a) = sum^{n-1) (ab)^i
  while (R.isZero(a)) R.random(a, rb);
  R.inv(b, a); R.addin(b, R.one());
  s = R.one();
  t = R.zero();
  for(int i = 0; i < reps; i++) 
  { R.addin(t, s); 
    R.mul(r, b, s);
    R.mul(s, r, a); }
  // t = sum^{n-1} b^i a^i, s = b^n a^n  (n is reps)
  R.subin(s, R.one()); // s = b^n a^n - 1
  //R.mul(r, a, b);
  //R.subin(r, R.one()); // r = ab - 1
  R.divin(s, a);
  if (! R.areEqual(s, t) ) cerr << ++e << ". sum-powers bug\n";


  if ( e > 0 ) cerr << "  Field test FAIL\n";
  else         cerr << "  Field test OK\n";
  return e;
}

/********

class SetTrait{};
class AdditiveGroupTrait{};
class RingTrait{};
class Ring1Trait{};
class CR1Trait{};
class FieldTrait{};

template<class trait>
int test()
{}

...

  Integer p; R.characteristic(p);
  cout << endl << "begin: with p = " << p << " and n = " << n << endl;
  if (p != n ) 
    cout << "note that sums are not based on the characteristic\n";

  elt a,b;
  int i;

  a = G.zero();

  for (int i = 0; i < 10; i++) { b = G.random(); if (! G.isZero(b) ) break; }
  if ( i == 10 ) cerr << ++e << ". isZero or random bug\n;


  for (i = 1; i < n; i++) 
  {
    R.addin(a, R.one());
    R.addin(s, a);
  }
  cout << "sum i from 1 to n-1 (should be (n-1)n/2 = 0 mod p) ";
  R.write(cout, s);
  cout << endl;

  a = R.zero();
  s = R.zero();
  for (i = 1; i < n-1; i++) 
  {
    R.addin(a, R.one());
    R.addin(s, a);
  }
  cout << "sum i from 1 to n-2 (should be (n-2)(n-1)/2 = 1 mod p) ";
  R.write(cout, s);
  cout << endl;

  a = R.zero();
  s = R.zero();
  for (i = 1; i < n-2; i++) 
  {
    R.addin(a, R.one());
    R.addin(s, a);
  }
  cout << "sum i from 1 to n-3 (should be (n-3)(n-2)/2 = 3 mod p) ";
  R.write(cout, s);
  cout << endl;

  a = R.zero();
  s = R.zero();
  for (i = 1; i < (n+1)/2; i++) 
  {
    R.addin(a, R.one());
    R.addin(s, a);
  }
  a = R.one(); R.addin(a, a); R.mulin(a, a); R.addin(a, a); // a is 8.
  R.mulin(s, a);
  cout << "sum i from 1 to (n-1)/2 multiplied by ";
  R.write(cout, a); cout << " (should be (n-1)(n+1) = -1 mod p) ";
  R.write(cout, s);
  cout << endl;

  a = R.zero();
  s = R.zero();
  for (i = 1; i < (n+1)/2; i++) 
  {
    R.addin(a, R.one());
    R.addin(s, R.mul(a,a));
  }
    //R.write(cout, a);  cout << " (n-1)/2 " << endl;
  a = R.one(); 
    //R.write(cout, a);  cout << " 1 " << endl;
  R.addin(a, a); 
    //R.write(cout, a);  cout << " 2 " << endl;
  R.mulin(a, a); 
    //R.write(cout, a);  cout << " 4 " << endl;
    //R.write(cout, R.add(R.one(), R.one())); cout << " adding 2 ones " << endl;
  R.addin(a, R.add(R.one(), R.one())); // a is 6.
    //R.write(cout, a);  cout << " 6 " << endl;
  R.mulin(s, a);
  cout << "sum i^2 from 1 to (n-1)/2 multiplied by ";
  R.write(cout, a); cout << " (should be (n-1)(n+1)n/4 = 0 mod p) ";
  R.write(cout, s);
  cout << endl;

  const elt a0 = R.random();
  if (R.isZero(a0)) cout << "random produced a zero\n";
  a = R.one();
  s = R.zero();
  for (i = 0; i < n; i++) 
  {
    R.addin(s, a);
    R.mulin(a, a0);
  }
  cout << "sum a^i for i from 1 to n-1), for random a (";
  R.write(cout, a0); cout << " in this case), \n";
  cout << "  (should be (a^n - 1)/(a-1) = 1 if n = char p) ";
  R.write(cout, s);
  cout << endl;

  cout << "end:  with p = " << p << " and n = " << n << endl;

}
*/
