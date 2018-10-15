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
      return false;
    }
    }
  return true;
}




template <class Field>
bool genData (BlasMatrix<Field> &A, size_t bits){
  typename Field::Element ZZ;
  typedef typename Field::RandIter RandIter;    
  RandIter RI(ZZ,bits,5);//RandIter RI(ZZ) ;
  LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
  RDM.randomFullRank(A);
}

template <class Field>
bool genData (SparseMatrix<Field> &A, size_t bits){
  typename Field::Element ZZ;
  typedef typename Field::RandIter RandIter;    
  RandIter RI(ZZ,bits,5);//RandIter RI(ZZ) ;
  LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
  RDM.randomFullRank(A);
}

template <>
bool genData (BlasMatrix<Givaro::ZRing<Integer> >  &A, size_t bits){
  Givaro::ZRing<Integer> ZZ;
  typedef typename  Givaro::ZRing<Integer> ::RandIter RandIter;    
  RandIter RI(ZZ,bits,5);//RandIter RI(ZZ,bits) ;
  LinBox::RandomDenseMatrix<RandIter, Givaro::ZRing<Integer> >  RDM(ZZ,RI);
  RDM.randomFullRank(A);
}
template <>
bool genData (SparseMatrix<Givaro::ZRing<Integer> >  &A, size_t bits){
  Givaro::ZRing<Integer> ZZ;
  typedef typename  Givaro::ZRing<Integer> ::RandIter RandIter;    
  RandIter RI(ZZ,bits,5);//RandIter RI(ZZ,bits) ;
  LinBox::RandomDenseMatrix<RandIter, Givaro::ZRing<Integer> >  RDM(ZZ,RI);
  RDM.randomFullRank(A);
}

template <class Field>
bool genData (BlasVector<Field>  &B, size_t bits){
  typename Field::Element ZZ;
  typedef typename Field::RandIter RandIter;    
  RandIter RI(ZZ,bits,5);//RandIter RI(ZZ) ;
  B.random(RI);
}


template <>
bool genData (DenseVector<Givaro::ZRing<Integer> >  &B, size_t bits){
  Givaro::ZRing<Integer> ZZ;
  typedef typename  Givaro::ZRing<Integer> ::RandIter RandIter;    
  RandIter RI(ZZ,bits,5);//RandIter RI(ZZ,bits) ;
  B.random(RI);
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
  std::cout<<"Computation is done over Q"<<std::endl;
#ifdef __LINBOX_HAVE_MPI
  std::cout << "MPI solveCRA" << std::endl;
#else
  std::cout << "Sequential solveCRA" << std::endl;
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

    std::cout << "CPU time (seconds): " << endtime-starttime << std::endl;
#else
    std::cout << "CPU time (seconds): " << end-start << std::endl;
#endif
    
    tag=checkResult (ZZ, A, B, X2, d);

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
  size_t bits,niter,ni,nj;  
  
  bits=10, niter=1, ni=1,nj=1; 
  
  static Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
    { 'b', "-b B", "Set the mxaimum number of digits of integers to generate.", TYPE_INT,     &bits },
    { 'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT,     &niter },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); 
#ifdef __LINBOX_HAVE_MPI

  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);

#else
  srand (time(NULL));
#endif
  nj=ni;
    
  Givaro::ZRing<Integer> ZZ;
  DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);
  
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;
  DenseVector X(ZZ, A.rowdim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());

  for(long j=0;j<(long)niter;j++){  
    
#ifdef __LINBOX_HAVE_MPI  
    if(0==Cptr->rank()){
#endif

      genData (A, bits);
      genData (B, bits);


//	B.write(std::cout << " >>>> Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
//	A.write(std::cout << " >>>> Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;


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
		 )) break;
  }


#ifdef __LINBOX_HAVE_MPI
  delete Cptr;  
#endif
    return 0;
}

