// =========================================================
// (C) The Linbox Group 1999   Examples for using fields
// Fri Feb  8 14:10:35 MET 2002 / Gilles Villard 
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
   *  with a given type that one may choose between 
      LinBox kernel capabilities : 
   */ 

  typedef GMP_Rational_Field  Field;     // Rational numbers based on Gmp
  //typedef unparam_field<double> Field;   // Or the doubles via the LinBox interface 

  //typedef param_modular Field;             // Or the % on the integers for a mod p field  


  /* Once this type "Field" is chosen, a domain K is declared */
    
  Field K;
  //Field K(7);  /* If choice of the modular field, modulo 7 */

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
