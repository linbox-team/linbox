// =========================================================
// (C) The Linbox Group 1999
// Examples for using fields 
// Fri Feb  8 14:28:08 MET 2002 / Gilles Villard
// =========================================================

// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------
/* LinBox is connected to external libraries through 
 * wrappers. Here, wrappers for Givaro and NTL are included */
/* Need to have the corresponding prefexes at installation ! */


#include "LinBox/givaro_gfq.h"
//#include "LinBox/ntl.h"


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


  //Givaro parameterized modulo domain  
  //givaro_zpz<Std16> K(4);
  givaro_gfq K(7,3);

  // NTL arbitrary precision real field
  // (Could be parameterized by the precision)  
  //unparam_field<NTL::RR> K;
  //NTL::RR::SetPrecision(500);
  //NTL::RR::SetOutputPrecision(50);

  // NTL modulo p field 
  //unparam_field<NTL::zz_p> K;   
  //NTL::zz_p::init(553);

  fct(K);

  return 0;
};
