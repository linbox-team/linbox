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
/* LinBox is connected to external libraries through 
 * wrappers. Here, wrappers for Givaro and NTL are included */

#include "LinBox/lin_zpz_giv.h"
#include "LinBox/ntl.h"

using namespace LinBox;
 // ---------------------------------------------

template <class Field> 
int fct(const Field K) {
 
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


  // Givaro parameterized modulo domain  
  //ZpzDom<Std16> K(4);

  // NTL arbitrary precision real field
  // (Could be parameterized by the precision)  
  unparam_field<NTL::RR> K;
  NTL::RR::SetPrecision(500);
  NTL::RR::SetOutputPrecision(50);

  fct(K);

  return 0;
};
