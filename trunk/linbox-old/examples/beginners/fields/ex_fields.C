// =========================================================
// (C) The Linbox Group 1999   Examples for using fields
// Thu Sep 27 11:33:58 MEST 2001 / Gilles Villard
// Fri Feb  8 14:24:31 MET 2002
// =========================================================


// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------
#include "LinBox/param_modular.h"

using namespace LinBox;
 
/*  The template function "fct" reads two elements "a" and "b" of the 
 *  field "K" from the standard input and writes a/b on the standard output */ 

template <class Field> 
int fct(const Field&  K) {
 
  /* "K" is a field domain (a C++ object) of type "Field" (here the template 
   *  parameter). The type of the elements of "K" is obtained through 
   * "Field" by typename Field::element */

  typedef typename Field::element K_elt;

  K_elt a,b,r; 

  K.init(a); K.init(b); K.init(r);
  K.read(cin,a);  K.read(cin,b);
  K.div(r,a,b);
  K.write(cout,r) << "\n";
}

// ---------------------------------------------

int main() {

  /* Using the parameterized domain capabilities, several domains 
   * representing integers modulo may be used simultaneously. */

  param_modular D(3), K(7);

  fct(K);  fct(D);
  return 0;
};
