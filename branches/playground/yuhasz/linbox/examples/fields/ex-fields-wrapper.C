/** @name examples/fields/ex-fields-wrapper.C
 * @author Gilles Villard
 * @memo wrapping a field implementation of another system.
 * @doc
 * LinBox is connected to external libraries through 
 * wrappers. Here, a wrapper for  NTL is included */
/* Need to have the corresponding prefixes at installation ! */

 */
//@{
// =========================================================
// (C) The Linbox Group 1999
// Examples for using fields 
// Fri Feb  8 14:28:08 MET 2002 
// Wed Apr 17 17:31:24 MEST 2002/ Gilles Villard
// =========================================================

// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
// ---------------------------------------------

#include "linbox/field/ntl.h"

using namespace LinBox;

// ---------------------------------------------

template <class Field> 
int fct(const Field K) {
 
  typedef typename Field::Element K_elt;
  K_elt a,b,r; 
  K.init(a); K.init(b); K.init(r);
  K.read(cin,a);
  K.read(cin,b);
  K.div(r,a,b);
  K.write(cout,r) << "\n";
}

// ---------------------------------------------

/// no command line args
int main() {

  // NTL arbitrary precision real field
  // (Could be parameterized by the precision)  

  UnparametricField<NTL::RR> K;
  NTL::RR::SetPrecision(500);
  NTL::RR::SetOutputPrecision(50);

  // NTL modulo p field 

  //UnparametricField<NTL::zz_p> K;   
  //NTL::zz_p::init(553);

  fct(K);

  return 0;
};
//@}
