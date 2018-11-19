/* Copyright (C) 2018 The LinBox group
 * Written by Hongguang Zhu <zhuhongguang2014@gmail.com>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file tests/test-solveCRA.C
 * @ingroup benchmarks
 * @brief Testing the OMP parallel/serial rational solver
 */
#include <cassert>


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "linbox/linbox-config.h"
#include "givaro/modular.h"
#include "givaro/zring.h"
#include "linbox/matrix/sparse-matrix.h"

#include "linbox/algorithms/cra-domain-omp.h" //<---Only compile without MPI

#include "linbox/solutions/methods.h"
#include "linbox/solutions/solve.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;
using namespace std;


template <class Field, class Matrix, class Vector>
static bool checkResult (const Field  &ZZ,
			 Matrix &A,
			 Vector &B,
			 Vector &X,
			 Integer &d)
{
  Vector B2(ZZ, A.coldim());
  Vector B3(ZZ, A.coldim());
  A.apply(B2,X);
  for (size_t j = 0 ; j < B.size() ; ++j) B3.setEntry(j,d*B.getEntry(j));
  
  for (size_t j = 0 ; j < A.coldim() ; ++j)
    if(!ZZ.areEqual(B2[j],B3[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      return false;
    }
  return true;
}


template <class Field, class Matrix>
void genData (const Field  &ZZ, Matrix  &Mat, size_t bits, int q){
  typedef typename Field::RandIter RandIter;  
  if(q<0){  
  RandIter RI(ZZ, bits) ;
  LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
  RDM.randomFullRank(Mat);
  }else{
  RandIter RI(ZZ, bits, q);
  LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
  RDM.randomFullRank(Mat);
  }
  

}


template <class Field>
void genData (const Field  &ZZ, BlasVector<Field>  &Vec, size_t bits, int q){
  typedef typename Field::RandIter RandIter; 
  if(q<0){       
  RandIter RI(ZZ, bits);
  Vec.random(RI);
  }else{
  RandIter RI(ZZ, bits, q);
  Vec.random(RI);
  }
}



template <class Field, class Matrix, class Vector>
bool test_set(const Field  &ZZ, Vector &X2,
	      Matrix &A, Vector &B
	      ){
  bool tag = false;
//  Givaro::ZRing<Integer> ZZ;
  Integer d;//  Givaro::ZRing<Integer>::Element d;
  std::cout<<"Computation is done over Q"<<std::endl;
  std::cout << "OMP solveCRA" << std::endl;

  
  /***********************
    Results verification 
  ***********************/
  RingCategories::IntegerTag tg;
  

  Timer chrono;
  chrono.start();

  solveCRA (X2, d, A, B, tg, 
	    Method::BlasElimination()
	    //Method::Hybrid(*Cptr)
	    );	
  
  chrono.stop();

  //  DenseVector B2(ZZ, A.coldim());
  

    std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl; 
    
    tag=checkResult (ZZ, A, B, X2, d);

  return tag;
}



int main(int argc, char ** argv)
{
  

  int bits,bitsize, niter,ni,nj,nt,n, q;
  bool peak = false, loop=false;
  bits=10, niter=1, ni=3,nj=3,nt=1, n=3, q=-1; 
  
  static Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &n },
    { 'b', "-b B", "Set the bitsize for input value.", TYPE_INT,     &bitsize },
    { 'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT,     &niter },
    { 't', "-t T", "Set the number of threads to run unit tests.", TYPE_INT,     &nt },
    { 'q', "-q Q", "Set negative value for random input or positive value to always generate the same input.", TYPE_INT,     &q },
    { 'l', "-l L", "Set if the infinte testing loop should be applied.", TYPE_BOOL,     &loop },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); 
  
  
  omp_set_num_threads(nt);

  nj=n;
  ni=n;  

  Givaro::ZRing<Integer> ZZ;  
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;




  for(int j=0;loop || j<niter;j++){  
            

            std::cout << " Test with dimension: " << ni << " x " << ni << std::endl;
            std::cout << " Test with bitsize: " << bits << std::endl;
        {
          DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);
          DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
          genData (ZZ, A, bits, q);
          genData (ZZ, B, bits, q);

/*
	std::cerr << ">>>>Compute with B: " << std::endl;      
	for(long j=0;j<(long)nj;j++) std::cerr << B.getEntry(j) << std::endl; 
	
	A.write(std::cout << ">>>>Compute with A: " << A.rowdim() << " by " << A.coldim() << "\n"<< "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
*/
    
        if(!test_set(ZZ, X2, A, B )) break;
        }
            ni = rand() % n + 1;
            if (ni < n / 2 && ni % 2 == 0 && !peak) ni = 1;

            bits = rand() % bitsize + 1;
            if (bits < bitsize / 2 && bitsize % 2 == 0 && !peak) bits = 1;

            nj=ni;
            peak = !peak;

  }

    return 0;

  
}


