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

/*! @file benchmarks/benchmark-solveCRA.C
 * @ingroup benchmarks
 * @brief Benchmarking the MPI parallel rational solver
 */

#define __LINBOX_HAVE_MPI

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <linbox/linbox-config.h>
#include <givaro/modular.h>
#include <givaro/zring.h>
#include <linbox/matrix/sparse-matrix.h>


#ifdef __LINBOX_HAVE_MPI
#include <mpi.h>
#include "linbox/util/mpicpp.h"
//#include "linbox/util/mpi-gmp.inl"
#endif
#include <linbox/solutions/methods.h>
#include <linbox/solutions/solve.h>
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;
using namespace std;

#include <iostream>
#include <fstream>
///////////////////////////////////////////////////////////////////////////////////////////////////////////
checkResult (ZZ, A, B, X2, d);
template <class Field>
static bool checkResult (const Field  &ZZ,
				BlasMatrix<Field> &A,
				BlasVector<Field> &B,
				BlasVector<Field> &X,
				Integer &d){
BlasMatrix<Field> B2(ZZ, A.coldim());
A.apply(B2,X);
for (size_t j = 0 ; j < A.coldim() ; ++j) B3.setEntry(j,d*B.getEntry(j));
for (size_t j = 0 ; j < A.coldim() ; ++j)
  for (size_t j = 0 ; j < A.coldim() ; ++j){
    if(!ZZ.areEqual(B[j],B2[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      break;
    }
  }

    return true;
}
template <class Field>
static bool checkResult (const Field  &ZZ,
				    BlasVector<Field> &B,
				    BlasVector<Field> &B2){

    for (size_t j = 0 ; j < B.size() ; ++j)
      if(!ZZ.areEqual(B[j],B2[j])){
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }

    return true;
}

template <class Field>
static bool genData (const Field  &ZZ,   BlasMatrix<Field> &A, size_t bits){
    typedef typename Field::RandIter RandIter;    
    RandIter RI(ZZ,bits) ;
    LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
    RDM.randomFullRank(A);
}
/*
template <class Field>
static bool genData (const Field  &ZZ,   DenseMatrix<Field> &A, size_t bits){
    typedef typename Field::RandIter RandIter;    
    RandIter RI(ZZ,bits) ;
    LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
    RDM.randomFullRank(A);
}
template <class Field>
static bool genData (const Field  &ZZ,   SparseMatrix<Field> &A, size_t bits){
    typedef typename Field::RandIter RandIter;    
    RandIter RI(ZZ,bits) ;
    LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
    RDM.randomFullRank(A);
}
*/
template <class Field>
static bool genData (const Field  &ZZ,   BlasVector<Field>  &B, size_t bits){
    typedef typename Field::RandIter RandIter;    
    RandIter RI(ZZ,bits) ;
    B.random(RI);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
  
#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
#endif
  size_t bits,ni,nj;  
  
  bits=10, ni=3,nj=3; 

  static Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
    { 'b', "-b M", "Set the number of bits of integers to generate.", TYPE_INT,     &bits },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); 
#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif
  nj=ni;
  Givaro::ZRing<Integer> ZZ;
  DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);
  
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;
  DenseVector X(ZZ, A.rowdim()), X2(ZZ, A.rowdim()),  B(ZZ, A.rowdim());
  Givaro::ZRing<Integer>::Element d;
  
#ifdef __LINBOX_HAVE_MPI  
  if(0==Cptr->rank()){
#endif

    genData (ZZ, A, bits);
    genData (ZZ, B, bits);


    //LinBox::rank (r, A); std::cout<<"The rank of generated matrix A is:"<<r<<std::endl;  
    
    /*
      std::cerr << ">>>>Compute with B: " << std::endl;      
      for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
      
      std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
      if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
    */  
#ifdef __LINBOX_HAVE_MPI 	
  }//End of BLock for process(0)
#endif
  
#ifdef __LINBOX_HAVE_MPI
  //distribute big integer compatible data
  {
    //double starttime, endtime; 
    //MPI_Barrier(MPI_COMM_WORLD);
    //starttime = MPI_Wtime();
    //MPI data distribution for Integer type value
    Cptr->bcast(A,0);   // MPIgmpBcast(A, ni, nj, 0, Cptr);
    Cptr->bcast(B,0);   // MPIgmpBcast(B, ni, 0, Cptr);
    //MPI_Barrier(MPI_COMM_WORLD);
    //endtime   = MPI_Wtime(); 
    //std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
  }
#endif
  
  
  //Check if data are correctly distributed to all processes
  /*
    if(0!=Cptr->rank()){std::cerr <<"process("<<Cptr->rank()<< ")<<<<Compute with B: " << std::endl;
    for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; }
    
    std::cerr <<"process("<<Cptr->rank()<< ")<<<<Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
    if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
  */
  
  
#ifdef __LINBOX_HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
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
  Timer chrono;
  chrono.start();
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
  chrono.stop();
#endif
//  DenseVector B2(ZZ, A.coldim());

#ifdef __LINBOX_HAVE_MPI
  if(0 == Cptr->rank()){  
#endif
    /*
      std::cerr << "Compute with B: " << std::endl;
      for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
      
      std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
      if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A::=",Tag::FileFormat::Maple) << ';' << std::endl;
    */
   
    // solveCRA
    /* 
       std::cout << "MPI CRA Solution is [";
       for(DenseVector::const_iterator it=X2.begin();it != X2.end(); ++it)
       ZZ.write(cout, *it) << " ";
       std::cout << "] / ";
       ZZ.write(std::cout, d) << std::endl;
    */
#ifdef __LINBOX_HAVE_MPI
    std::cout << "CPU time (seconds): " << endtime-starttime << std::endl;
#else
    std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;
#endif


checkResult (ZZ, A, B, X2, d);

#ifdef __LINBOX_HAVE_MPI    
    MPI_Finalize();
#endif
    return 0;
#ifdef __LINBOX_HAVE_MPI
  }
  MPI_Finalize();
#endif
  
}


















































