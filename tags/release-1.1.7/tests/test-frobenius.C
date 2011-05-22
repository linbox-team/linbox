/* tests/test-direct-sum.C
 * Copyright (C) LinBox
 * Written by Austin Lobo, David Saunders
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

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
    { 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
    { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
    { 'i', "-i I", "Perform each test for I iterations.",   TYPE_INT,     &iterations1 },
    { 'j', "-j J", "Apply test matrix to J vectors.",         TYPE_INT,     &iterations2 },
	{ '\0' }
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

	commentator.start("Frobenius form black box test suite", "frobenius");
  Frobenius<Field>  A(F, plist.begin(), plist.end());

  pass = pass && testBlackbox(A);
  
	commentator.stop("Frobenius form black box test suite");
  return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
