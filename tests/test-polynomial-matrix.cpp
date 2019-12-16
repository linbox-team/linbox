/*
 * Copyright (C) 2013  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */




#include <iostream>
using namespace std;

#include <linbox/ring/modular.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/randiter/random-fftprime.h>
#include <givaro/zring.h>
#include <recint/rint.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/polynomial-matrix.h>
//#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>
#include <linbox/algorithms/polynomial-matrix/matpoly-mult-naive.h>

using namespace LinBox;


int main(int argc, char** argv){
  static size_t  m = 10; // matrix dimension
  static size_t  n = 15; // matrix dimension
  static size_t  d = 50; // polynomial size
  static long    seed = time(NULL);

  static Argument args[] = {
			    { 'm', "-n M", "Set row dimension of test matrices to M.", TYPE_INT,     &m },
			    { 'n', "-n N", "Set col dimension of test matrices to N.", TYPE_INT,     &n },
			    { 'd', "-d D", "Set degree of test matrices to D.", TYPE_INT,     &d },
			    { 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
			    END_OF_ARGUMENTS
  };
  parseArguments (argc, argv, args);


  commentator().start ("Testing polynomial matrix", "testMatpoly", 1);


  typedef Givaro::Modular<double> Field;

  Field F(65537);
  typename Field::RandIter G(F,seed);
  typedef PolynomialMatrix<Field, PMType::polfirst> MatrixP;
  //typedef PolynomialMatrix<Field, PMType::matfirst> PMatrix;


  PolynomialMatrixNaiveMulDomain<Field> PMMD(F);
  
  
  ostream& report = LinBox::commentator().report();
  report<<"Polynomial matrix (polfirst) testing over ";F.write(report)<<std::endl;


  MatrixP A1(F,m,n,d),B1(F,n,n,d), C1(F,m,n,2*d-1);
  MatrixP A2(F,m,n,d),B2(F,n,n,d), C2(F,m,n,2*d-1);

  PMMD.mul(C1,A1,B1);
  
  bool pass= true;   
  commentator().stop(MSG_STATUS(pass),(const char *) 0,"testMatpoly");
    
  return (pass? 0: -1);
} 



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
