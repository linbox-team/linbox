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
 * @ingroup test
 * @brief Testing the MPI parallel/serial rational solver
 */
//#define __Detailed_Time_Measurement
#define __LINBOX_HAVE_MPI




#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "linbox/linbox-config.h"
#include "givaro/modular.h"
#include "givaro/zring.h"
#include "linbox/matrix/sparse-matrix.h"

#ifdef __LINBOX_HAVE_MPI
#include <mpi.h>
#include "linbox/util/mpicpp.h"	
#endif
#include "linbox/solutions/methods.h"
#include "linbox/solutions/solve.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;
using namespace std;

template <class Field, class Matrix>
static bool checkResult (const Field  &F,
			 Matrix &A,	
			 BlasVector<Field> &B,
			 BlasVector<Field> &X,
			 typename Field::Element &d)
{

  BlasVector<Field> B2(F, A.coldim());
  BlasVector<Field> B3(F, A.coldim());
  A.apply(B2,X);
  
  for (size_t j = 0 ; j < B.size() ; ++j){
    B3.setEntry(j,d*B.getEntry	(j));
  }
  for (size_t j = 0 ; j < A.coldim() ; ++j){
    if(!F.areEqual(B2[j],B3[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cerr << " B2["<<j<<"] := "<< B2[j] << std::endl;
      std::cerr << " B3["<<j<<"] := "<< B3[j] << std::endl;
      std::cerr << " d*B["<<j<<"] := "<< d*B.getEntry(j) << std::endl;
      
      return false;
    }
  }
  return true;
}

/////////////////////////////////////////////////////
template <class Field, class Matrix>
void genData (const Field  &F, Matrix  &Mat, size_t bits, int seed){
  typedef typename Field::RandIter RandIter;  
  
  RandIter RI(F, bits, seed) ;
  LinBox::RandomDenseMatrix<RandIter,Field>  RDM(F,RI);
  RDM.randomFullRank(Mat);
  
}
/////////////////////////////////////////////////////
template <class Field>
void genData (const Field  &F, BlasVector<Field>  &Vec, size_t bits, int seed){
  typedef typename Field::RandIter RandIter; 
  
  RandIter RI(F, bits, seed);
  Vec.random(RI);
  
}
/////////////////////////////////////////////////////
template <class Field>
bool test_with_field(BlasVector<Field> &X2,
		     BlasMatrix<Field> &A,
		     BlasVector<Field> &B
		     ,int bits

#ifdef __LINBOX_HAVE_MPI
		     , Communicator *Cptr
#endif
		     
		     ){
  bool tag = false;
  Givaro::ZRing<Integer> F;
  Givaro::ZRing<Integer>::Element d;
  std::cerr<<"Computation is done over Q"<<std::endl;
  
#ifdef __LINBOX_HAVE_MPI
  std::cerr << "MPI solveCRA" << std::endl;
#endif 

  
  
  

  /***********************
    Results verification 
  ***********************/
  RingCategories::IntegerTag tg;
  
#ifdef __LINBOX_HAVE_MPI
  double starttime, endtime;
  starttime = MPI_Wtime(); 
#else
  //Timer chrono;
  double start;
  double end;
  
  start = omp_get_wtime();
#endif
  solveCRA (X2, d, A, B, tg, 
	    Method::BlasElimination()
	    //Method::Hybrid(*Cptr)
#ifdef __LINBOX_HAVE_MPI
	    ,Cptr
#endif
	    );	
  
#ifdef __LINBOX_HAVE_MPI
  endtime   = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(0 == Cptr->rank())  std::cout << "Total CPU time (seconds): " << endtime-starttime << std::endl;
#else
  end = omp_get_wtime();
  std::cout << "Total CPU time (seconds): " << end-start << std::endl;
#endif

#ifdef __LINBOX_HAVE_MPI
  if(0 == Cptr->rank()) {
#endif
  tag=checkResult (F, A, B, X2, d);
#ifdef __LINBOX_HAVE_MPI
  }
#endif
  
#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&tag, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif

  return tag;	      
}
/////////////////////////////////////////////////////
//subroutine for defaut seed if no user provided seed as parameter
uint64_t getSeed(){
  struct timeval tp;
  gettimeofday(&tp, 0) ;
  return static_cast<uint64_t> (tp.tv_usec + tp.tv_sec*1000000);
}


void rand_param_vary(size_t &n, size_t &ni, size_t &bits, size_t &bitsize, bool &peak){
  
  n = ni;
  
  n = rand() % ni + 1;
  if (n < ni / 2 && n % 2 == 0 && !peak) n = 1;
  
  bits = rand() % bitsize + 1;
  if (bits < bitsize / 2 && bitsize % 2 == 0 && !peak) bits = 1;
  
  peak = !peak;
  
}
/////////////////////////////////////////////////////
template <class Field>
void prepare_data_with_field(size_t bits, int seed,
			     BlasVector<Field> &X2,
			     BlasMatrix<Field> &A,
			     BlasVector<Field> &B
			     
#ifdef __LINBOX_HAVE_MPI
			     , Communicator *Cptr
#endif
			     ){
  Field F;
  
#ifdef __LINBOX_HAVE_MPI  
  if(0==Cptr->rank()){
#endif

    std::cerr << " Test with seed: "<<seed<<std::endl;
    
    genData (F, A, bits, seed);
    genData (F, B, bits, seed);
    
#ifdef __LINBOX_HAVE_MPI 	
  }//End of BLock for process(0)
#endif   
  
  //	B.write(std::cout << " Proc("<<Cptr->rank()<<") >>>> Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
  //	A.write(std::cout << " Proc("<<Cptr->rank()<<") >>>> Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;
  

#ifdef __LINBOX_HAVE_MPI
  //distribute big integer compatible data
  {
    
#ifdef __Detailed_Time_Measurement
    double starttime, endtime; 
    MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();
#endif
    
    //MPI data distribution for Integer type value
    Cptr->bcast(A,0);
    Cptr->bcast(B,0);
    
#ifdef __Detailed_Time_Measurement
    MPI_Barrier(MPI_COMM_WORLD);
    endtime   = MPI_Wtime(); 
    std::cout<<"In Proc("<<Cptr->rank()<<") MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
#endif
    
  }
#endif
  
  //Check if data are correctly distributed to all processes
  //	B.write(std::cout << " Proc("<<Cptr->rank()<<") <<<< Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
  //	A.write(std::cout << " Proc("<<Cptr->rank()<<") <<<< Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;
  
  
}
/////////////////////////////////////////////////////
template <class Field>
void update_input(BlasMatrix<Field>& A, BlasVector<Field>& B, BlasVector<Field>& X2,
		  int& seed, int q, size_t& n, size_t ni, size_t& bits, size_t bitsize, bool& peak
#ifdef __LINBOX_HAVE_MPI
		  , Communicator *Cptr
#endif
		  ){
#ifdef __LINBOX_HAVE_MPI
  if(0==Cptr->rank()){
#endif
    
    if(q<0){   
      seed = getSeed();
      rand_param_vary(n, ni, bits, bitsize, peak);
    }else{
      bits = bitsize;  
      n = ni;
    }
  std::cout << " Test with dimension: " << n << " x " << n << std::endl;
  std::cout << " Test with bitsize: " << bits << std::endl;
#ifdef __LINBOX_HAVE_MPI 	
  }
#endif    


#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bits, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  A.resize(n,n);  X2.resize( A.coldim()),  B.resize(A.coldim());
        
  prepare_data_with_field(bits, seed, X2, A, B	     
#ifdef __LINBOX_HAVE_MPI
			  , Cptr
#endif
			  );
}
/////////////////////////////////////////////////////

void get_input_param_ready(int& seed, int& q, size_t& n, size_t& ni, size_t& bits, size_t& bitsize, bool& peak, bool& loop, size_t& nt
#ifdef __LINBOX_HAVE_MPI
			   , Communicator *Cptr
#endif
			   ){
#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&loop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&q, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  
#ifdef __LINBOX_HAVE_MPI
  if(0==Cptr->rank()){
#endif
    
    if(q<0){   
      seed = getSeed();
      rand_param_vary(n, ni, bits, bitsize, peak);
    }else{
      bits = bitsize;  
      n = ni;
    }
    
#ifdef __LINBOX_HAVE_MPI 	
  }
#endif    
  
  
#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bits, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}
/////////////////////////////////////////////////////


int main(int argc, char ** argv)
{

#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv, MPI_THREAD_SERIALIZED);//Cptr = new Communicator(&argc, &argv);
#endif

  size_t bits=10;
  size_t bitsize=10;
  size_t niter=1;
  size_t ni=1;
  size_t nt=1;
  
  size_t n=1;
  int q=-1;
  bool peak = false;
  bool loop=false;
  int seed = 0; 

  
  Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
    { 'b', "-b B", "Set the bitsize for input value.", TYPE_INT,     &bitsize },
    { 'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT,     &niter },
    { 'q', "-q Q", "Set the randomness of test (<0 for random and >0 for derterministic).",    TYPE_INT , &q },
    { 'l', "-loop Y/N", "Set if run the test in an infinite loop.", TYPE_BOOL , &loop },
    { 's', "-s S", "Set the seed to fill the input matrices.", TYPE_INT,     &seed },
    { 't', "-t T", "Set the number of threads per process.", TYPE_INT,     &nt },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); 
  
  
  Givaro::ZRing<Integer> F;
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;
  
  
  get_input_param_ready(seed,  q,  n,  ni,  bits,  bitsize,  peak, loop, nt
#ifdef __LINBOX_HAVE_MPI
			,  Cptr
#endif
			);
  
  
  DenseMatrix<Givaro::ZRing<Integer> > A (F,n,n);
  DenseVector X(F, A.coldim()), X2(F, A.coldim()),  B(F, A.coldim()); 



  
  for(size_t j=0;loop || j<niter;j++){  
    
    update_input(A, B, X2, seed, q, n, ni, bits, bitsize, peak
#ifdef __LINBOX_HAVE_MPI
		 , Cptr
#endif
		 );
    
    omp_set_num_threads(nt);
    if(!test_with_field<Givaro::ZRing<Integer>>(X2, A, B, bits

#ifdef __LINBOX_HAVE_MPI
						, Cptr
#endif
						
						)){
      break;
    }

    
  }
  
  
#ifdef __LINBOX_HAVE_MPI
  delete Cptr;  
#endif

  return 0;
  
}

