// =========================================================
// (C) The Linbox Group 1999   Examples for using fields
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

int main() {

  /*  Use of LinBox basic field capabilities.
   *  A concrete field will be a domain (an object)
   *  with a given type that one may choose using: 
   */ 

  typedef GMP_Rational_Field  Field;
  //typedef unparam_field<double> Field;

  //typedef param_modular Field; 

  /* Once this type "Field" is chosen, a domain K */
    
  //Field K(7);  /* Modulo 7 */
  Field K;

  /* The element of "K" are addressed using the type "Field::element" */ 

  Field::element a,b,r;

  /* "a" and "b" elements of K, initialized, read from 
   * std input, divided, written on std output 
   */ 

  K.init(a); K.init(b); K.init(r);
  K.read(cin,a);  K.read(cin,b);

  K.div(r,a,b);
  K.write(cout,r) << "\n";

  return 0;
};
