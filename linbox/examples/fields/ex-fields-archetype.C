// =========================================================
// (C) The Linbox Group 1999
// Examples for using fields 
// Fri Feb  8 16:25:22 MET 2002 
// Wed Apr 17 17:37:12 MEST 2002/ Gilles Villard
// =========================================================

// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
// ---------------------------------------------


//#include "LinBox/gmp-rational-field.C"
//#include "LinBox/gmp-rational-random.h"

//#include "LinBox/abstract_double.h"
#include "linbox/field/unparametric.h"
//#include "LinBox/param_modular.h"

//#include "LinBox/givaro_zpz.h"
#include "linbox/field//ntl.h"

#include "linbox/field/archetype.h"

using namespace LinBox;
 
// ---------------------------------------------
/*  The template function "fct" reads two elements "a" and "b" of the 
 *  field "K" from the standard input and writes a/b on the standard output */ 

template <class Field> 
int fct(const Field& K) {
 
  typedef typename Field::Element K_elt;

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
  //ZpzDom<Std16> K_o(4);
  //GMP_Rational_Field  Q_o;

  //UnparametricField<NTL::RR> K_o;
  //NTL::RR::SetPrecision(400);
  //NTL::RR::SetOutputPrecision(50);

  UnparametricField<NTL::zz_p> K_o;   
  NTL::zz_p::init(553);

  /* These field objects "K_o" and "Q_o" of different types can be converted to 
   * objects Q and K of a unique type "Field_archetype" for instance using 
   * a constructor: */ 

  //Field_archetype Q( & Q_o );
  FieldArchetype K( & K_o );

  /* The template function "fct" is called with two different fields but the 
   * template is instantiated only once since it is called with a unique 
   * template parameter "Field_archetype" */

  //fct(Q);
  fct(K);

  return 0;
};
