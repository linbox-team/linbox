// =========================================================
// (C) The Linbox Group 1999
// Examples for using fields 
// Fri Feb  8 16:25:22 MET 2002 / Gilles Villard
// =========================================================

#define _REENTRANT   

// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------

#include "LinBox/field_archetype.h"

#include "LinBox/gmp-rational-field.C"
#include "LinBox/gmp-rational-random.h"

#include "LinBox/abstract_double.h"
#include "LinBox/unparam_field.h"
#include "LinBox/param_modular.h"

#include "LinBox/givaro_zpz.h"
#include "LinBox/ntl.h"

using namespace LinBox;
 
// ---------------------------------------------
/*  The template function "fct" reads two elements "a" and "b" of the 
 *  field "K" from the standard input and writes a/b on the standard output */ 

template <class Field> 
int fct(const Field& K) {
 
  typedef typename Field::element K_elt;

  K_elt a,b,r; 

  K.init(a); K.init(b); K.init(r);
  K.read(cin,a);
  K.read(cin,b);
  K.div(r,a,b);
  K.write(cout,r) << "\n";
  }

// ---------------------------------------------

int main() {

  /* The field objects "K_o" and "Q_o" are constructed as in previous examples 
   */ 

  //abstract_double K_o;
  //unparam_field<double> K_o;
  //param_modular K_o(4);
  ZpzDom<Std16> K_o(4);

  //GMP_Rational_Field  Q_o;

  //unparam_field<NTL::RR> K_o;
  //NTL::RR::SetPrecision(400);
  //NTL::RR::SetOutputPrecision(50);

  /* These field objects "K_o" and "Q_o" of different types can be converted to 
   * objects Q and K of a unique type "Field_archetype" for instance using 
   * a constructor: */ 

  //Field_archetype Q( & Q_o );
  Field_archetype K( & K_o );

  /* The template function "fct" is called with two different fields but the 
   * template is instantiated only once since it is called with a unique 
   * template parameter "Field_archetype" */

  //fct(Q);
  fct(K);

  return 0;
};
