/** \file examples/integer-mul.C
 * \author Gilles Villard
 * \brief The LinBox arbitrary precision integer type illustrated.
\ingroup examples
 *
 * class `integer' is a wrapper of GMP integers.
 */
// ============================================================
// (C) The Linbox Group 1999   Examples for using long integers
// Fri Feb  8 14:00:35 MET 2002 Gilles Villard 
// Wed Apr 17 15:41:35 MEST 2002
// Sun Feb  9 23:15:15 MET 2003
// ============================================================

// ---------------------------------------------
#include <iostream>
#include <fstream> 
// ---------------------------------------------

#include "linbox-config.h"

// Use of Gmp based LinBox integers 
#include "linbox/integer.h"

using namespace LinBox;
using namespace std;
 
// ---------------------------------------------

/// no command line args.  Prompts for two integers.
int main() {

  integer a,b;

  cout << "1st integer > ";
  cin >> a;
  cout << "2nd integer > ";
  cin >> b;

  cout << "The product " << a*b << "\n";

  return 0;
};
