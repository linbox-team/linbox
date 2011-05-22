/* Copyright (C) LinBox
 *
 * Copyright (C) 2003 Austin Lobo, B. David Saunders
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#include <iostream>
#include <fstream>
#include <vector>
#include "linbox/field/ntl.h"

#include <linbox/field/ntl-ZZ_p.h>
#include <linbox/integer.h>
#include <linbox/blackbox/ntl-hankel.h>


#include "test-generic.h"


using namespace std;


int main(int argc, char* argv[])
{
  LinBox::commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);
  ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
  
  bool pass = true;
  
  static size_t n = 1000;
  static long q = 2147483647;


  static int iterations = 1;
  
  static Argument args[] = {
    { 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
    { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
    { 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
	{ '\0' }
  };
  
  parseArguments (argc, argv, args);
  

  
  //------ Read q and construct F(q)
  NTL::ZZ modulus; 	// prime modulus
  modulus = q;
  //std::cout << std::endl << "Enter a prime number for the modulus of the field: ";
  //std::cin >> modulus;
  report <<  "The modulus is " << modulus << std::endl;
  NTL::ZZ_p::init(modulus); // NOTE: This is essential for using NTL
  

	commentator.start("Hankel black box test test suite", "Hankel");
  report << "\tn= " <<  n << " \tq= " << q <<   endl ;

  typedef LinBox::UnparametricField<NTL::ZZ_p> Field;
  typedef Field::Element element;
  typedef std::vector<element> Vector;  
  
  // Now we are using the NTL wrapper as the field, call the instance F
  Field F;
  element zero;
  F.init(zero, 0);
  
  // Use the default constructor to create a matrix
  LinBox::Hankel<Field> T;
  
  // Use a special constructor to construct a matrix of dim TSIZE
  int TSIZE = 2*n-1;
  Vector tdata(TSIZE); 
  report << "The random vector is:" << std::endl;
  for (unsigned int i=0; i < tdata.size(); i++) {
    tdata[i] = NTL::random_ZZ_p() ;
    report << tdata[i] << " ";
  }
  report << std::endl;
  
  LinBox::Hankel<Field> TT(F,tdata);
  report << "The matrix is: " << std::endl;
  TT.print(report);
  
  // Create an interesting input vector called idata
  Vector idata((TSIZE+2)/2), odata((TSIZE+2)/2);
  report << "A random col vector:\t" << std::endl;
  for (unsigned int i=0; i < idata.size(); i++) {
    idata[i] = NTL::random_ZZ_p() ;
    report << idata[i] << " ";
  }
  report << std::endl;
  
  // Apply the matrix to the vector just created
  // Testing the apply function when both input and output are over ZZ_p
  TT.applyTranspose(odata, idata);
  report << "\tTesting apply Transpose:----------------- \nResult is[";
  for (unsigned int i = 0; i < odata.size(); i++)
    report << odata[i] << " ";
  report << "]\n";
  
  
  
  TT.apply(odata, idata);
  report << "\n\nTesting  apply :--------------------- \nResult is[";
  for (unsigned int i = 0; i < odata.size(); i++)
    report << odata[i] << " ";
  report << "]\n";
  
   report << "Setting the matrix to UniModular Lower Triangular";
   TT.setToUniModLT();
   TT.print(report);
   report << "\nSetting the matrix to UniModular Upper Triangular";
      TT.setToUniModUT();
      TT.print(report);
  
  pass = testBlackbox(TT);

	commentator.stop("Hankel black box test test suite");
  return pass ? 0 : -1;
  
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
