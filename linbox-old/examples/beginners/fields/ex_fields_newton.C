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
#include "LinBox/unparam_field.h"
#include "LinBox/gmp-rational-field.C"

using namespace LinBox;
 // ---------------------------------------------
/* Newton iterations to invert a number */

template <class Field> 
int fct(const Field K) {
 
  typedef typename Field::element K_elt;

  K_elt a,w,r,two; 
  K.init(a); K.init(w); K.init(r); K.init(two,2);

  long nb_iter=0,i;

  cout << "Number to invert > "; 
  K.read(cin,a);
  cout << "\n";

  cout << "Initialization point > "; 
  K.read(cin,w);
  cout << "\n";

  cout << "Number of iterations > "; 
  cin >> nb_iter;
  cout << "\n";

  for (i=0; i<nb_iter; i++){
     K.mul(r,a,w);
     K.sub(r, two, r);
     K.mulin(w,r);
  } 


  K.write(cout,w) << "\n";
}

// ---------------------------------------------

int main() {

  // Rational numbers from Gmp 
  typedef GMP_Rational_Field  Field;

  // Or double numbers wrapped to LinBox 
  //typedef unparam_field<double> Field;

  // NTL arbitrary precision real field
  // (Could be parameterized by the precision)  

  //typedef unparam_field<NTL::RR> Field;
  //NTL::RR::SetPrecision(500);
  //NTL::RR::SetOutputPrecision(50);

  Field K; 

  fct(K);

  return 0;
};
