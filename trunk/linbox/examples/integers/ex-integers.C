// ============================================================
// (C) The Linbox Group 1999   Examples for using long integers
// Fri Feb  8 14:00:35 MET 2002 / Gilles Villard 
// ============================================================

// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------

// Use of Gmp based LinBox integers 
#include "LinBox/integer.h"

using namespace LinBox;
 
// ---------------------------------------------

int main() {

  integer a,b;

  cout << "1rst integer > ";
  cin >> a;
  cout << "2nd integer > ";
  cin >> b;

  cout << "The product " << a*b << "\n";

  return 0;
};
