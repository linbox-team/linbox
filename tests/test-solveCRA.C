/* Copyright (C) 2018 The LinBox group
 *
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
static bool checkResult (const Field  &F,
			 Matrix &A,
			 Vector &B,
			 Vector &X,
			 typename Field::Element &d
			 )
{
  Vector B2(F, A.coldim());
  Vector B3(F, A.coldim());
  A.apply(B2,X);
  for (size_t j = 0 ; j < B.size() ; ++j) B3.setEntry(j,d*B.getEntry(j));
  
  for (size_t j = 0 ; j < A.coldim() ; ++j)
    if(!F.areEqual(B2[j],B3[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      return false;
    }
  return true;
}


template <class Field, class Matrix>
void genData (const Field  &F, Matrix  &Mat, size_t bits, int seed=0){
  typedef typename Field::RandIter RandIter;  
std::cout << " Test with seed: " << seed << std::endl;
  RandIter RI(F, bits, seed);
  LinBox::RandomDenseMatrix<RandIter,Field>  RDM(F,RI);
  RDM.randomFullRank(Mat);

}


template <class Field>
void genData (const Field  &F, BlasVector<Field>  &Vec, size_t bits, int seed=0){
  typedef typename Field::RandIter RandIter; 

  RandIter RI(F, bits, seed);
  Vec.random(RI);

}


template <class Field, class Matrix, class Vector>
bool test_set(const Field  &F, Vector &X2,
	      Matrix &A, Vector &B
	      ){
  bool tag = false;
  typename Field::Element d;
  std::cout<<"Computation is done over Q"<<std::endl;
  std::cout << "OMP solveCRA" << std::endl;

  
  /***********************
    Results verification 
  ***********************/
 typename ClassifyRing<Field>::categoryTag tg;//RingCategories::IntegerTag tg;
  

  Timer chrono;
  chrono.start();
  
  solveCRA (X2, d, A, B, tg, 
	    Method::BlasElimination()
	    //Method::Hybrid(*Cptr)
	    );	
  
  chrono.stop();
  
  //  DenseVector B2(F, A.coldim());
  
  
  std::cout << "Real time (seconds): " << chrono.realtime() << std::endl; 
  
  tag=checkResult (F, A, B, X2, d);
  
  return tag;
}

uint64_t getSeed(){
	struct timeval tp;
	gettimeofday(&tp, 0) ;
    return static_cast<uint64_t> (tp.tv_usec + tp.tv_sec*1000000);
}

int main(int argc, char ** argv)
{
  
    int seed=1; // This value should be other than 0 if not the generated input will always be random
    size_t nt = 1;
    size_t n = 100;
    size_t bitsize = 10;
    size_t niter = 3;

    int q = -1;
    bool peak = false, loop=false;
  
  static Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &n },
    { 'b', "-b B", "Set the bitsize for input value.", TYPE_INT,     &bitsize },
    { 'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT,     &niter },
    { 't', "-t T", "Set the number of threads to run unit tests.", TYPE_INT,     &nt },
    { 'q', "-q Q", "Set negative value for random input or positive value to always generate the same input.", TYPE_INT,     &q },
    { 's', "-s S", "Set the seed to always generate the same input.", TYPE_INT,     &seed },
    { 'l', "-l L", "Set if the infinte testing loop should be applied.", TYPE_BOOL,     &loop },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); 
    
  omp_set_num_threads(nt);
 
  Givaro::ZRing<Integer> F;  
  
  typedef Givaro::ZRing<Integer> TF;
  
  
  typedef BlasVector<TF> DenseVector;
  DenseMatrix< TF > A (F,n,n);
  DenseVector X2(F, A.coldim()),  B(F, A.coldim());
  size_t ni=n;
  size_t bits=bitsize;

  for(size_t j=0;loop || j<niter;j++){  
    
    
    std::cout << " Test with dimension: " << ni << " x " << ni << std::endl;
    std::cout << " Test with bitsize: " << bits << std::endl;

    A.resize(ni,ni);
    B.resize(ni,ni);
    X2.resize(ni);
    if(q<0){
        genData (F, A, bits, getSeed());
        genData (F, B, bits, getSeed());

    }else{
        genData (F, A, bits, seed);
        genData (F, B, bits, seed);

    }
   /*
	std::cerr << ">>>>Compute with B: " << std::endl;      
	for(long j=0;j<(long)ni;j++) std::cerr << B.getEntry(j) << std::endl; 
	
	A.write(std::cout << ">>>>Compute with A: " << A.rowdim() << " by " << A.coldim() << "\n"<< "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
   */
omp_set_num_threads(nt);    
   if(!test_set(F, X2, A, B )) break;

    if(q<0){
        ni = rand() % n + 1;
        if (ni < n / 2 && ni % 2 == 0 && !peak) ni = 1;
       
        bits = rand() % bitsize + 1;
        if (bits < bitsize / 2 && bitsize % 2 == 0 && !peak) bits = 1;
    } 

    peak = !peak;

  }
  
  return 0;
  
}


