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
#include "linbox/util/mpicpp.h"	//#include "linbox/util/mpi-gmp.inl"
#else
//#include "linbox/algorithms/cra-domain-omp.h" //<---Only compile without MPI
#endif
#include "linbox/solutions/methods.h"
#include "linbox/solutions/solve.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;
using namespace std;


template <class Field, class Matrix>
static bool checkResult (const Field  &ZZ,
			 Matrix &A,
			 BlasVector<Field> &B,
			 BlasVector<Field> &X,
			 Integer &d)
{
  BlasVector<Field> B2(ZZ, A.coldim());
  BlasVector<Field> B3(ZZ, A.coldim());
  A.apply(B2,X);

Integer tmp;
  for (size_t j = 0 ; j < B.size() ; ++j){
    B3.setEntry(j,d*B.getEntry(j));
  }
  for (size_t j = 0 ; j < A.coldim() ; ++j){
    if(!ZZ.areEqual(B2[j],B3[j])){
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




bool test_set(BlasVector<Givaro::ZRing<Integer> > &X2,
	      BlasMatrix<Givaro::ZRing<Integer> > &A,
	      BlasVector<Givaro::ZRing<Integer> > &B
#ifdef __LINBOX_HAVE_MPI
	      , Communicator *Cptr
#endif
	      ){
  bool tag = false;
  Givaro::ZRing<Integer> ZZ;
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
#else
  end = omp_get_wtime();
#endif
  //  DenseVector B2(ZZ, A.coldim());

  
#ifdef __LINBOX_HAVE_MPI
  if(0 == Cptr->rank()){  

    std::cout << "Total CPU time (seconds): " << endtime-starttime << std::endl;
#else
    std::cout << "Total CPU time (seconds): " << end-start << std::endl;
#endif

    tag=checkResult (ZZ, A, B, X2, d);
    if(!tag){
        B.write(std::cout << " >>>> Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
	    A.write(std::cout << " >>>> Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;
	    std::cout << " >>>> Found d:="<<d<< std::endl;
	    std::cout << " >>>> Found X:=";
	    for(long i=0; i<X2.size();i++)std::cout << X2[i] << std::endl;
    }
#ifdef __LINBOX_HAVE_MPI
  }
#endif 
  
#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&tag, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif
  return tag;
}



int main(int argc, char ** argv)
{
  
#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
#endif
  int bits,bitsize, niter,ni, nj, nt,n,q;
  bool peak = false, loop=false;
  bits=10, niter=1, ni=3, nj=3, nt=1, n=3, q=-1; 
  
  static Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
    { 'b', "-b B", "Set the bitsize for input value.", TYPE_INT,     &bitsize },
    { 'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT,     &niter },
    { 't', "-t T", "Set the number of threads to run unit tests.", TYPE_INT,     &nt },
    { 'q', "-q Q", "Set negative value for random input or positive value to always generate the same input.", TYPE_INT,     &q },
    { 'l', "-l L", "Set if the infinte testing loop should be applied.", TYPE_BOOL,     &loop },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); 
  

  

///////////////////////////////////////////////////////////////////////////////////////////////////////
#if 1

#ifdef __LINBOX_HAVE_MPI

  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

#else
  srand (time(NULL));
#endif
  nj=ni;
    
  Givaro::ZRing<Integer> ZZ;
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;

  DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);
  DenseVector X(ZZ, A.rowdim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
  
  for(long j=0;j<(long)niter;j++){  
    
#ifdef __LINBOX_HAVE_MPI  
    if(0==Cptr->rank()){
#endif

    bits = rand() % bitsize + 1;
    if (bits < bitsize / 2 && bitsize % 2 == 0 && !peak) bits = 1;
    
    peak = !peak;
    
    std::cout << " Test with dimension: " << ni << " x " << nj << std::endl;
    std::cout << " Test with bitsize: " << bits << std::endl;


      genData (ZZ, A, bits, q);//genData (A, bits);
      genData (ZZ, B, bits, q);//genData (B, bits);






#ifdef __LINBOX_HAVE_MPI 	
    }//End of BLock for process(0)
#endif

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
//	B.write(std::cout << " <<<< Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
//	A.write(std::cout << " <<<< Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl; 


	
    if(!test_set(X2, A, B
#ifdef __LINBOX_HAVE_MPI
		 , Cptr
#endif
		 )){

		 break;
		 }

  }
#else /////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&q, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&loop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif

    
  Givaro::ZRing<Integer> ZZ;
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;

  for(int j=0;loop || j<niter;j++){  


#ifdef __LINBOX_HAVE_MPI  
    if(0==Cptr->rank()){
#endif


    
    n = rand() % ni + 1;
    if (n < ni / 2 && n % 2 == 0 && !peak) n = 1;

    
    bits = rand() % bitsize + 1;
    if (bits < bitsize / 2 && bitsize % 2 == 0 && !peak) bits = 1;
    
    peak = !peak;
    
    std::cout << " Test with dimension: " << n << " x " << n << std::endl;
    std::cout << " Test with bitsize: " << bits << std::endl;
    
#ifdef __LINBOX_HAVE_MPI 	
    }//End of BLock for process(0)
#endif    



#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&bits, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&peak, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD); 
#endif


    DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,n,n);
    DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());


#ifdef __LINBOX_HAVE_MPI  
    if(0==Cptr->rank()){
#endif
      genData (ZZ, A, bits, q);
      genData (ZZ, B, bits, q);

#ifdef __LINBOX_HAVE_MPI 	
    }//End of BLock for process(0)
#endif   


//	B.write(std::cout << " >>>> Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
//	A.write(std::cout << " >>>> Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;




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
//	B.write(std::cout << " <<<< Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
//	A.write(std::cout << " <<<< Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl; 


	
    if(!test_set(X2, A, B
#ifdef __LINBOX_HAVE_MPI
		 , Cptr
#endif
		 )){
		 break;
		 }
		

    

  }
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __LINBOX_HAVE_MPI
  delete Cptr;  
#endif
    return 0;
}

