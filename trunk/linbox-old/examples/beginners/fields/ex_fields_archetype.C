// =========================================================
// (C) The Linbox Group 1999
// Examples for using fields 
// Thu Sep 27 11:33:58 MEST 2001 / Gilles Villard
// =========================================================

#define _REENTRANT  // pour le rand_r ?? 
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

#include "LinBox/lin_zpz_giv.h"
#include "LinBox/ntl.h"



using namespace LinBox;
 
// ---------------------------------------------

template <class Field> 
int fct(const Field& K) {
 
  typedef typename Field::element K_elt;

  K_elt a,b,r; 

  K.read(cin,a);
  K.read(cin,b);

  K.div(r,a,b);

  K.write(cout,r) << "\n";

  }

// ---------------------------------------------

int main() {

  //GMP_Rational_Field  K;

  //abstract_double K;

  //unparam_field<double> K;

  //param_modular K(4);

  //ZpzDom<Std16> K(4);

  unparam_field<NTL::RR> K;
  NTL::RR::SetPrecision(500);
  NTL::RR::SetOutputPrecision(50);

  //fct<  unparam_field<NTL::RR>  > (K);

  Field_archetype AK( & K );

  fct< Field_archetype >(AK);

  return 0;
};
