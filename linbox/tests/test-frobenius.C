/* -*- mode: C++; tab-width: 2; indent-tabs-mode: t; c-basic-offset: 2 -*- */

/* tests/test-direct-sum.C
 * Written by Austin Lobo, David Saunders
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/frobenius.h"

#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
  bool pass = true;
  
  static size_t n = 10;
  static integer q = 101;
  static int iterations1 = 100;
  static int iterations2 = 1;
  
  static Argument args[] = {
    { 'n', "-n N", "Set dimension of test matrices to NxN (default 10)", TYPE_INT,     &n },
    { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)", TYPE_INTEGER, &q },
    { 'i', "-i I", "Perform each test for I iterations (default 100)",   TYPE_INT,     &iterations1 },
    { 'j', "-j J", "Apply test matrix to J vectors (default 1)",         TYPE_INT,     &iterations2 },
  };
  
  typedef Modular<uint32> Field;
  typedef vector<Field::Element> Vector;
  typedef Vector Polynomial;
  typedef vector<Polynomial> Plist;
  
  parseArguments (argc, argv, args);
  Field F (q);
  Plist plist(3);
  
  Field::RandIter r(F);
  
  size_t  pdeg = 10;
  plist[0].resize(pdeg+1);
  for ( size_t ideg=0; ideg < pdeg; ++ideg) r.random(plist[0][ideg]);
  F.init(plist[0][pdeg],1);

  pdeg = 6;
  plist[1].resize(pdeg+1);
  for ( size_t ideg=0; ideg < pdeg; ++ideg) r.random(plist[1][ideg]);
  F.init(plist[1][pdeg],1);

  pdeg = 4;
  plist[2].resize(pdeg+1);
  for ( size_t ideg=0; ideg < pdeg; ++ideg) r.random(plist[2][ideg]);
  F.init(plist[2][pdeg],1);

  cout << endl << "black box frobenius-form test suite" << endl;
  Frobenius<Field>  A(F, plist.begin(), plist.end());

  pass = pass && testBlackbox(F, A);
  //pass = pass && testSmallBlackbox(F, A);
  
  return pass ? 0 : -1;
}
