// =========================================================
// (C) The Linbox Group 1999
// Examples for using fields 
// Thu Sep 27 11:33:58 MEST 2001 / Gilles Villard
// =========================================================


// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------

#include "LinBox/field_archetype.h"
#include "LinBox/gmp-rational-field.C"
#include "LinBox/abstract_double.h"

using namespace LinBox;
 
// ---------------------------------------------

template <class Field> 
int in_a_field(const Field&  K) {
 
  typedef typename Field::element K_elt;

  K_elt a,b,r; 

  K.read(cin,a);
  K.read(cin,b);

  K.div(r,a,b);

  K.write(cout,r);

}

// ---------------------------------------------

int main() {

  //GMP_Rational_Field  K;

  abstract_double K;

  in_a_field<  abstract_double > (K);

  return 0;
};
