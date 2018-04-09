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

/*! @file tests/test-mpi-comm.C
 * @ingroup benchmarks
 * @brief Check MPI communicator interface
 */
#define __LINBOX_HAVE_MPI //! Necessary to compile correctly for the communicator !
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <linbox/linbox-config.h>
#include <givaro/modular.h>
#include <givaro/zring.h>
#include <linbox/matrix/sparse-matrix.h>


#include <mpi.h>
#include "linbox/util/mpicpp.h"
//#include "linbox/util/mpi-gmp.inl"


#include <linbox/solutions/methods.h>
#include <linbox/solutions/solve.h>
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;
using namespace std;
/*
template<class T>
T& myrand (T& r, long size)
{
  if (size < 0)
    return r = T( (lrand48() % (-size-size)) + size );
  else
    return r = T(  lrand48() % size ) ;
};

#include <gmp++/gmp++.h>
#include <string>

std::string gmp_rand ( size_t maxNdigits)
{ 
  std::string result, tmpStr;
  long tmp; 
  tmpStr = std::to_string(myrand(tmp, 10));
  while(result.size()+tmpStr.size()<maxNdigits){
    result+=tmpStr;
    tmpStr = std::to_string(myrand(tmp, 10));
  }
  return result;
}
*/
///////////////////////////////////////////////////////////////////////////////////////////
template <class Field>
static bool checkResult (const Field  &ZZ,
				    BlasMatrix<Field> &A,
				    BlasMatrix<Field> &A2){
  for (long i = 0 ; i < A.rowdim() ; ++i)  
    for (long j = 0 ; j < A.coldim() ; ++j){
      if(!ZZ.areEqual(A[i][j],A2[i][j])){
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }
    } 

    return true;
}
template <class Field>
static bool checkResult (const Field  &ZZ,
				    SparseMatrix<Field> &A,
				    SparseMatrix<Field> &A2){
  for (long i = 0 ; i < A.rowdim() ; ++i)  
    for (long j = 0 ; j < A.coldim() ; ++j){
      if(!ZZ.areEqual(A.getEntry(i,j),A2.getEntry(i,j))){
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class Field>
bool genData (BlasMatrix<Field> &A, size_t bits){
    typename Field::Element ZZ;
    typedef typename Field::RandIter RandIter;    
    RandIter RI(ZZ) ;
    LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
    RDM.randomFullRank(A);
}

template <class Field>
bool genData (SparseMatrix<Field> &A, size_t bits){
    typename Field::Element ZZ;
    typedef typename Field::RandIter RandIter;    
    RandIter RI(ZZ) ;
    LinBox::RandomDenseMatrix<RandIter,Field>  RDM(ZZ,RI);
    RDM.randomFullRank(A);
}

template <class Field>
bool genData (BlasVector<Field>  &B, size_t bits){
    typename Field::Element ZZ;
    typedef typename Field::RandIter RandIter;    
    RandIter RI(ZZ) ;
    B.random(RI);
}


template <>
bool genData (BlasMatrix<Givaro::ZRing<Integer> >  &A, size_t bits){
    Givaro::ZRing<Integer> ZZ;
    typedef typename  Givaro::ZRing<Integer> ::RandIter RandIter;    
    RandIter RI(ZZ,bits) ;
    LinBox::RandomDenseMatrix<RandIter, Givaro::ZRing<Integer> >  RDM(ZZ,RI);
    RDM.randomFullRank(A);
}
template <>
bool genData (SparseMatrix<Givaro::ZRing<Integer> >  &A, size_t bits){
    Givaro::ZRing<Integer> ZZ;
    typedef typename  Givaro::ZRing<Integer> ::RandIter RandIter;    
    RandIter RI(ZZ,bits) ;
    LinBox::RandomDenseMatrix<RandIter, Givaro::ZRing<Integer> >  RDM(ZZ,RI);
    RDM.randomFullRank(A);
}
template <>
bool genData (DenseVector<Givaro::ZRing<Integer> >  &B, size_t bits){
    Givaro::ZRing<Integer> ZZ;
    typedef typename  Givaro::ZRing<Integer> ::RandIter RandIter;    
    RandIter RI(ZZ,bits) ;
    B.random(RI);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void test_main(size_t bits, size_t ni, size_t nj, Communicator *Cptr)
{

  Givaro::ZRing<T> ZZ;
  DenseMatrix<Givaro::ZRing<T>> A (ZZ,ni,nj), A2(ZZ,ni,nj);

  typedef BlasVector<Givaro::ZRing<T> > DenseVector;
  DenseVector B2(ZZ, A.coldim()),  B(ZZ, A.coldim());


  if(0==Cptr->rank()){

    genData (A, bits);
    genData (B, bits);

B.write(std::cout << ">>>>Compute with B:\n",Tag::FileFormat::Maple) << ';' << std::endl;
A.write(std::cout << ">>>>Compute with A:\n",Tag::FileFormat::Maple) << ';' << std::endl;


  }//End of BLock for process(0)





if(0==Cptr->rank()){

//double starttime, endtime; 
//MPI_Barrier(MPI_COMM_WORLD);
//starttime = MPI_Wtime();
  //MPI data distribution for Integer type value
    Cptr->ssend(B,1);

//MPI_Barrier(MPI_COMM_WORLD);
//endtime   = MPI_Wtime(); 
//std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;

}else{

    Cptr->recv(B2,0);

}


if(0==Cptr->rank()){

//double starttime, endtime; 
//MPI_Barrier(MPI_COMM_WORLD);
//starttime = MPI_Wtime();
  //MPI data distribution for Integer type value

    Cptr->ssend(A,1); 
//MPI_Barrier(MPI_COMM_WORLD);
//endtime   = MPI_Wtime(); 
//std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;

}else{

    Cptr->recv(A2,0);
}



Cptr->bcast(B,0);
Cptr->bcast(A,0);

 if(0!=Cptr->rank()) checkResult (ZZ, A, A2);
 if(0!=Cptr->rank()) checkResult (ZZ, B, B2);

//Check if data are correctly distributed to all processes
  if(0!=Cptr->rank()){
B.write(std::cout << "process("<<Cptr->rank()<< ")<<<<Compute with B:\n",Tag::FileFormat::Maple) << ';' << std::endl;
A.write(std::cout << "process("<<Cptr->rank()<< ")<<<<Compute with A: \n",Tag::FileFormat::Maple) << ';' << std::endl; 
}

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char ** argv)
{
  
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
  size_t bits,ni,niter,nj;  

   bits=10,niter=3, ni=3,nj=3; 

 	static Argument args[] = {
		{ 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
		{ 'b', "-b B", "Set the mxaimum number of digits of integers to generate.", TYPE_INT,     &bits },
		{ 'i', "-i I", "Set the number of iteration over unit test sets.", TYPE_INT,     &niter },
		END_OF_ARGUMENTS
	};	
parseArguments (argc, argv, args); 

 MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); nj=ni;
 MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);

srand (time(NULL));

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
for(int j=0;j<niter;j++){
	test_main<float>(bits,ni,nj,Cptr);
	test_main<double>(bits,ni,nj,Cptr);
	test_main<int>(bits,ni,nj,Cptr);
	test_main<Integer>(bits,ni,nj,Cptr);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  MPI_Finalize();
  return 0;  
  
}


