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

#include "LinBox/gmp-rational-field.C"
#include "LinBox/unparam_field.h"
#include "LinBox/param_modular.h"

using namespace LinBox;
 
// ---------------------------------------------

template <class Field> 
int fct(const Field&  K) {
 
  typedef typename Field::element K_elt;

  K_elt a,b,r; 

  K.read(cin,a);
  K.read(cin,b);

  K.div(r,a,b);

  K.write(cout,r) << "\n";

}

// ---------------------------------------------

int main() {

  GMP_Rational_Field  K;

  //unparam_field<double> K;

  //param_modular K(4);

  fct< GMP_Rational_Field > (K);

  return 0;
};
