/* interfaces.h              Linbox            -bds 8/99, 4/00       */
/*-------------------------------------------------------------------*/
#ifndef LB_interfaces_h
#define LB_interfaces_h

/* Domain interfaces:

Use of these interfaces will allow for succinct 
documentation of the properties expected of domain objects in linbox code.  
They are akin to "concepts" in the Standard Template Library.  

For example, if an algorithm is declared as 

template <class Field>
void foobar(Field& F, const Field::element& a, Field::element& b, )
// sets b such that proposition P(a, b) holds.
{...}

then it is understood (because of the use of the name "Field" for the 
template parameter) that this proposition holds:
  if F satisfies the "Field_a" interface conditions,
  then foobar promises property P of the result b on input a.
Thus the algorithm writer has suscinctly stated some requirements on 
its input just by using the name "Field" for the domain.

It may be that the algorithm also works correctly when the domain class,
the template class parameter, does not instantiate a complete implementation 
of the Field interface.  The algorithm may be documented for such more 
specific, less restrictive conditions if wished.

The goal of this suite of interfaces is to capture the most widely
used domain types in linbox.

The interfaces are described as abstract C++ classes (member functions 
are pure virtual).  They can be thought of as concrete categories, 
in the mathematical sense of the term.
Classes which implement the interface may or may not inherit
from the interface class.
For example, if Field_a is the interface, then the assertion that a class
implements Field means that the class works in any template declared
as "template <class Field>" both in the sense that the code compiles 
and in the sense that the semantics asserted for Field functionality 
is implemented.
*/

/* Proposed interfaces for Linbox:

Integers/Integer, Set, AdditiveGroup, Ring/Ring1/CR1, Field, Module/FreeModule
...
( future:
ED PID UFD RichCR1 RichRing PIR LocalRing
)
*/

/* Naming conventions:

A domain usually has several variant implementations of a mathematical 
operation, "op", such as muliplication.  In this document we reserve the 
terms function and operator to speak of a C++ entity, and use the term 
"operation" when we speak of the mathematical "function", "mapping",  
or "operator" being implemented.  

1.  On functions, usually use a noun form (rank, minPoly, dotProd, etc.).
However, verb forms are appropriate for predicates (isZero) 
and we use verb forms by established custom for certain operations 
(add, sub, mul, div, apply).

2. Use suffix -in when a parameter is both input and output.
For example, consider mul and mulin forms.

    F.mul(r, a, b) does r <- a*b (and returns a reference to r). 
mul may use already allocated memory held by r, but must check
that it is sufficient.  mul may allocate new memory for the result.
(and is responsible to free any storage at r now unused - uh, what 
issues here if r contains pointers, possibly shared, into heap?
--- we can address this issue iff it comes up in domains we care about.)

    F.mulin(a, b) does a <- a*b (and returns a reference to a).
mulin may use already allocated memory held by a, but must check
that it is sufficient.  mulin knows to read before writing in a
as necessary and may be more efficient that mul.

    F.mul(a, b) returns a reference to a*b in new heap storage.
On second thought, this is an invitation for memory leaks.  Lets try:
    F.mul(a, b) returns a stack value a*b, to be copied...
and deprecate it's use.  It can be implemented entirely automatically
with an adaptor, and thus this form can be excised from the basic
domain requirements.  See DomainAllocatingForms adaptor.
Type variables that use this adaptor should have the suffix _al
as in "template<class Field_al>"

[[[[[[[[Design question:  We are dealing with two issues here.  One is whether
the user wishes to handle allocation (and reuse memory) for performance,
or wishes the field to handle allocation.  The other is to give the field
a chance to minimimize memory use, and possibly save time, by special
treatement of the mulin form.  However for other fields this is just a burden...
]]]]]]]]]

{\em Optionally} Domains may have operator forms of the mathematical ops.
This is "syntactic sugar".  I don't know how this can/should be done 
because many different domains may use the same element type (which
determines the operator signature).
When a template is depending on the operator forms, it should use the 
suffix _ss on the name of the template parameter.  
template<class Field_ss> ...  c = a*b ...;
*/
/* -------- end prolog -------------------------------------------------- */

//#include "Integers.h" 
  // defines Integers domain and it's element type Integer 
//#include "Set_a.h" 
  // defines set interface applying to all domains
//#include "AdditiveGroup.h" 
  // defines additive abelian group interface
//#include "Ring.h" 
  // defines Ring_a Ring1_a CR1_a interface
//#include "Field.h" 
  // defines Field_a interface
//#include "Module.h" 
  // defines Module_a, FreeModule_a interface

#endif
/* -------- cut -------------------------------------------------- */
/* Set.h                 Linbox                   -bds 3/00       */
/* ---------------------------------------------------------- */
#ifndef LB_Set_h
#define LB_Set_h

#include <iostream.h>
#include "Integers.h" 
// for Integers domain and it's element type: Integer
/** 
Set_a - interface definition for the category Set.

The Set members 
class element, constructors,  
functions are areEqual(), areEq(), random(), cardinality().
*/

template<class Rep>
class Set_a 
{public: 
 /** 
 The type of the underlying set of elements
 */
  typedef Rep element;

  /** 
  Constructors
  Domains must have null constructor, copy constructor, assignment op.
  If other constructors are defined, these must be also be.
  */
  //Set_a(); Set_a(Set_a S); 

  /** 
  Label briefly identifying the implementation.
  */
  virtual char* label()const=0;

// set element equality. 
  /** 
  areEqual(a,b) is true iff a and b represent same set element.
  */
  virtual bool areEqual(const element& a, const element& b)const=0; 

  /// areEq(a,b) is true iff a and b are identical in storage.
  virtual bool areEq(const element& a, const element& b)const=0;
  
// IO 
  /** 
  read(is, a) Reads a char stream format, allocates and assigns to a. 
  */
  virtual istream& read(istream& is, element& a)const=0;

  /** 
  write(os, a) puts out a char string in format compatible with read. 
  */
  virtual ostream& write(ostream& os, const element& a)const=0;

  ///A Set_ss also defines ==, <<, >>

// size
  /** 
  cardinality() is the set size.  It is -1 if the set is infinite.
  */
  virtual Integer& cardinality(Integer& r)const=0;

// random
  /** 
  Assigns to r and returns a reference to a random element whose probability 
  is <= 1/n.   Argument n greater than set size is an error.  
  */ 
  virtual element& random(element& r, const Integer& n)const=0;

  /** 
  Assigns to r and returns a random value.
  Uses a default distribution.  
  For finite sets the default is the uniform distribution. 
  For infinite sets I don't know what the heck this is.
  */
  virtual element& random(element& r)const=0; 

  /** 
  Uses key to determine the distribution (nice idea of J-G).
  For instance if element is integers, distribution used might be 
  uniform over 0 .. key-1 if key is pos, 
  uniform over key/2 ..  -key/2, if key is neg.
  If element is univariate poly, degree of key could determine degree
  of result.  Coeff of key determine coeff range, etc.
  */
  virtual element& random(element& r, const element& key)const=0;

};

#endif
/* -------- cut -------------------------------------------------- */
/* AdditiveGroup.h                 Linbox                   -bds 3/00       */
/* ---------------------------------------------------------- */
#ifndef LB_AdditiveGroup_h
#define LB_AdditiveGroup_h

/** 
AdditiveGroup_a interface - definition for the category AbelianGroup
with add() as the basic binary operation.

zero(), isZero(), 
add(), neg(), sub(),
addin(), negin(), subin(),
and the Set functionality.
*/
template<class Rep>
class AdditiveGroup_a : virtual public Set_a<Rep>
{public: 

// Start with Set_a interface and add these:
  /** 
  zero() is the identity element: zero() + a == a, zero()*a == zero().
  */
  virtual const element& zero()const=0;   

  /** 
  isZero(a) // a == 0
  */
  virtual bool isZero(const element& a)const=0;

  /** 
  add(r,a,b) // r = a + b, commutative, associative.
  */
  virtual element& add(element& r, const element& a, const element& b)const=0; 

  /** 
  neg(r,a) // r = -a.  i.e. isZero(add(a,neg(a))).
  */
  virtual element& neg(element& r, const element& a)const=0; 

  /** 
  sub(r,a,b) // r = a + -b.
  */
  virtual element& sub(element& r, const element& a, const element& b)const=0; 

// in place forms
  /**
  addin(r,b) // r += b.
  */
  virtual element& addin(element& r, const element& b)const=0;  /// r += b;

  /**
  negin(r) // r = -r.
  */
  virtual element& negin(element& r)const=0;  

  /**
  subin(r,b) // r -= b.
  */
  virtual element& subin(element& r, const element& b)const=0;  

// scalar mul by Integer, (abelian gp is Z_Module)
  /**
  Zprod(r, n, a) // r <- n*a, scalar product
  */
  virtual element& Zprod(element& r, const Integer n, const element& a)const=0;

  /**
  Zprodin(n, a) // a <- n*a, scalar product
  */
  virtual element& Zprodin(const Integer n, element& a)const=0;

  /**
  n = characteristic() if n is least > 0 such that Zprod(n,a) = 0, for all a.
  0 = characteristic() if no such positive n.
  */
  virtual Integer& characteristic(Integer& r)const=0;

};

#endif
/* -------- cut -------------------------------------------------- */
/* Ring.h                 Linbox                   -bds 3/00       */
/* ---------------------------------------------------------- */
#ifndef LB_Ring_h
#define LB_Ring_h

/** 
Ring_a - interface definition for the category Ring.

A multiplicative identity is not required, see Ring1.
Multiplicative commutativity is not required, see CR1.

Ring functions are 

mul(), axpy(), axmy()
mulin(), axpyinx(), axpyiny(), axmyinx(), axmyiny(),
and AdditiveGroup functions.
*/
template<class Rep>
class Ring_a : virtual public AdditiveGroup_a<Rep>
{public: 

// Start with AdditiveGroup_a interface and add these:

// multiplicative semiring
  /**
  mul(r,a,b) // r = a * b, commutative, associative, distributes over +.
  */
  virtual element& mul(element& r, const element& a, const element& b)const=0; 

// mixed forms
  /**
  axpy(r,a,x,y) // r = a*x + y.
  */
  virtual element& axpy
    (element& r, const element& a, const element& x, const element& y)const=0;

  /**
  axmy(r,a,x,y) // r = a*x - y.
  */
  virtual element& axmy
    (element& r, const element& a, const element& x, const element& y)const=0;

// in place forms
  /**
  mulin(r,b) // r *= b.
  */
  virtual element& mulin(element& r, const element& b)const=0; 

  /**
  axpyinx(a,x,y) // x = a*x + y.
  */
  virtual element& axpyinx(const element& a, element& x, const element& y)const=0; 

  /**
  axpyiny(a,x,y) // y = a*x + y.
  */
  virtual element& axpyiny(const element& a, const element& x, element& y)const=0; 

  /**
  axmyinx(a,x,y) // x = a*x - y.
  */
  virtual element& axmyinx(const element& a, element& x, const element& y)const=0; 

  /**
  axmyiny(a,x,y) // y = a*x - y.
  */
  virtual element& axmyiny(const element& a, const element& x, element& y)const=0; 

/*
// vector forms
  typedef Vector<Self> Vect; /// should be iterator pair instead
  element& reduceadd(element& r, const Vect& v)const=0; 
  // r <-- sum of entries in v, possibly efficient.
  // Returns zero if v.length is 0.
*/

};

/** 
Ring1_a interface - definition for the category Ring1
of rings with multiplicative unit.

Multiplicative commutativity is not required, see CR1.

Ring1 functions are 

one(), isOne(), 
and the Ring functions.

Basic matrix operations that don't depend on multiplicative commutativity
should be declared over Ring1, i.e. using "template <class Ring1>".  
For example, when the entries are matrices (blocks)...
Also, matrix rings themselves are Ring1's.
*/
template<class Rep>
class Ring1_a 
: virtual public Ring_a<Rep>
{public: 

// start with Ring_a interface and add these:
  /**
  one()*a == a.
  */
  virtual const element& one()const=0;

  /**
  isOne(a) same as areEqual(a, one()).
  */
  virtual bool isOne(const element& a)const=0;

/*
  element& reducemul(element& r, const Vect& v)const=0; 
  // r <-- product of entries in v, possibly efficient.
  // Returns one if v.length is 0.

  dotprod(element& r, const Vect& u, const Vect& v)const=0;
  // r <-- inner product of u and v
  // i.e. r <-- reduceadd( M.mul(u,v) ), where M is the appropriate module.
*/
};

/** 
CR1_a interface - definition for the category CR1
of commutative rings with one.

CR1 functions are just
the Ring1 functions, but mul is commutative.

Some of our basic matrix operations are valid over CR1s, but not R1s.
They should use "template <class CR1>".
*/
template<class Rep>
class CR1_a : virtual public Ring1_a<Rep> 
{ // mul is commutative
  class Commutative{};
  // use a trait trick to prevent Ring1 matching a CR1 interface requirement
};

#endif
/* -------- cut -------------------------------------------------- */
/* Field.h                 Linbox                   -bds 3/00       */
/* ---------------------------------------------------------- */
#ifndef LB_Field_h
#define LB_Field_h

/** 
Field_a interface - definition for the category Field

Field functions are 
inv(), div(), invin(), divin(),
and the CR1 functions

The non-zero elements form an abelian group under mul(). inv(), one().
*/
template<class Rep>
class Field_a : virtual public CR1_a<Rep>
// should be refined thru ID (integral domain).
{ 

// start with CR1_a interface and add these:
  /**
  inv(r,a) // r <- 1/a, i.e. a * inv(a) = one().  Error if a is zero.
  */
  virtual element& inv(element& r, const element& a)const=0; 

  /**
  div(r,a,b) // r <- a / b.
  */
  virtual element& div(element& r, const element& a, const element& b)const=0; 

// inplace forms 
  /**
  invin(r) // r = 1/r.
  */
  virtual element& invin(element& r)const=0; 

  /**
  divin(r,a) // r /= a.
  */
  virtual element& divin(element& r, const element& a)const=0; 

};

#endif
/* -------- cut -------------------------------------------------- */
/* Module.h                 Linbox                   -bds 3/00       */
/* ---------------------------------------------------------- */
#ifndef LB_Module_h
#define LB_Module_h
namespace stl
{
#include <vector.h>
}

/** 
Module_a interface - definition for the category Module

Module functions are 
scprod(), scprodin(), (scalar multiplications)
and the AdditiveGroup functions

*/
template<class Ring>
class Module_a : virtual public AdditiveGroup_a<stl::vector<Ring> >
{public:

  typedef typename Ring::element scalar;
  /**
  Sprod(r, s, m) // r <- s*m, scalar product
  */
  virtual element& Sprod(element& r, const scalar& s, const element& m)const=0;

  /**
  Sprodin(s, m) // m <- s*m, scalar product
  */
  virtual element& Sprodin(const scalar& s, element& m)const=0;

};

/** 
FreeModule_a interface - definition for the category FreeModule

FreeModule functions are 
dotprod(), coefficient(), 
and the Module functions.
*/
template<class Ring>
class FreeModule_a : virtual public Module_a<Ring>
{public:
  // Module_a members plus these:
  class basis_element;

  // add
  virtual element& add(element& r, const element& a, const basis_element b)const=0;
  virtual element& add(element& r, const basis_element& a, const element b)const=0;
  virtual element& add(element& r, const basis_element& a, const basis_element b)const=0;

  // neg
  virtual element& neg(element& r, const basis_element b)const=0;

  // sub
  virtual element& sub(element& r, const element& a, const basis_element b)const=0;
  virtual element& sub(element& r, const basis_element& a, const element b)const=0;
  virtual element& sub(element& r, const basis_element& a, const basis_element b)const=0;

  // mul (scalar product)
  virtual element& Sprod(element& r, const scalar& s, const basis_element e)const=0;

  /**
  dotprod(r, a, b) r <- a.b (dot product aka innerproduct)
  note that dot with ith basis vector is ith coefficient access.
  */
  virtual scalar& dotprod(scalar& r, const element& a, const element b)const=0;
  virtual scalar& dotprod(scalar& r, const basis_element& a, const element b)const=0;
  virtual scalar& dotprod(scalar& r, const element& a, const basis_element b)const=0;
  virtual scalar& dotprod(scalar& r, const basis_element& a, const basis_element b)const=0;

  // need function integer to ith basis element
};
/* *******************
// matrices and vectors

template<class R>
class Vector
{
iterators
set and get
natural length
}

///
template<class R>
class BBMatrix  // base matrix class
{
  ///
  Vector<R> apply(Vector<R> v, const Vector<R> u)const=0;
  /// v = Au;
  Vector<R> apply_transpose(Vector<R> v, const Vector<R> u)const=0;
  /// v = uA;
  Integer nrows;
  /// natural number of rows.  In particular, higher rows are zero.
  Integer ncols; 
  // natural number of columns.  In particular, higher columns are zero.
};
*/
#endif
