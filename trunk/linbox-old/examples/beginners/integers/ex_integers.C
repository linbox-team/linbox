// ============================================================
// (C) The Linbox Group 1999   Examples for using long integers
// Thu Sep 27 11:33:58 MEST 2001 / Gilles Villard
// ============================================================

// ---------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------

#include "LinBox/integer.h"

using namespace LinBox;
 
// ---------------------------------------------

int main() {

  Integer a,b;

  cout << "1rst integer > ";
  cin >> a;
  cout << "2nd integer > ";
  cin >> b;

  cout << "The product " << a*b << "\n";

  return 0;
};
