/** \file examples/fields/ex-fields.C
 * \author Gilles Villard
 * \brief Using a template function with two distinct fields. 
 */

// =========================================================
// (C) The Linbox Group 1999   Examples for using fields
// Thu Sep 27 11:33:58 MEST 2001 / Gilles Villard
// Fri Feb  8 14:24:31 MET 2002
// =========================================================

// ---------------------------------------------

#include "linbox/linbox-config.h"

#include <iostream>

// ---------------------------------------------
#include "linbox/field/modular.h"
#include "linbox/field/ntl.h"

using namespace LinBox;
using namespace std;
 
/**  The template function "fct" reads two elements "a" and "b" of the 
 *  field "K" from the standard input and writes a/b on the standard output */ 

template <class Field> 
void divide_ex(const Field&  K) {
 
  /* "K" is a field domain (a C++ object) of type "Field" (here the template 
   *  parameter). The type of the elements of "K" is obtained through 
   * "Field" by typename Field::Element */

  typedef typename Field::Element K_elt;

  K_elt a,b,r; 

  K.init(a); K.init(b); K.init(r);
  cout << "division example: enter two numbers: ";
  K.read(cin,a);  K.read(cin,b);
  K.div(r,a,b);
  K.write( K.write(cout<< "the quotient in ") << " is ", 
           r) << endl;

}

// ---------------------------------------------
/// no command line args
int main() {

  /* Using the parameterized domain capabilities, several domains 
   * representing integers modulo may be used simultaneously. */

  Modular<uint32> D(3), K(7);

  divide_ex(D);  divide_ex(K);

  // NTL arbitrary precision real field
  // (Could be parameterized by the precision)  

  UnparametricField<NTL::RR> K2;
  NTL::RR::SetPrecision(500);
  NTL::RR::SetOutputPrecision(50);

  // NTL modulo p field 

  //UnparametricField<NTL::zz_p> K2;   
  //NTL::zz_p::init(553);

  divide_ex(K2);

  return 0;
};

