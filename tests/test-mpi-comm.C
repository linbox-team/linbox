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
////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////
#include <limits.h>
#include <float.h>
double get_random_double() {
  //return DBL_MIN * ((double) rand() / (double) RAND_MAX) - DBL_MAX;
    double f = (double)rand() / RAND_MAX;
    return DBL_MIN + f * (DBL_MAX - DBL_MIN);
}
float get_random_float() {
//  return FLT_MIN * ((double) rand() / (double) RAND_MAX) - FLT_MAX;
    float f = (float)rand() / RAND_MAX;
    return FLT_MIN + f * (FLT_MAX - FLT_MIN);
}
int get_random_int() {
  return INT_MIN * (rand() / RAND_MAX) - INT_MAX;
}
///////////////////////////////////////////////////////////////////////////////////

//#include <map>
#include <iostream>
#include <fstream>


int main(int argc, char ** argv)
{
  
  
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
  size_t dg,ni,nj,max;  

   dg=10, ni=3,nj=3,max=1000; 

 	static Argument args[] = {
		{ 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
		{ 'm', "-m M", "Set the mxaimum for the range of integers to generate.", TYPE_INT,     &max },
		{ 'd', "-d M", "Set the mxaimum number of digits of integers to generate.", TYPE_INT,     &dg },
		END_OF_ARGUMENTS
	};	
parseArguments (argc, argv, args); 

 MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); nj=ni;

  Givaro::ZRing<double> ZZ;
  DenseMatrix<Givaro::ZRing<double>> A (ZZ,ni,nj), A2(ZZ,ni,nj);

  typedef BlasVector<Givaro::ZRing<double> > DenseVector;
  DenseVector B2(ZZ, A.coldim()),  B(ZZ, A.coldim());
  Givaro::ZRing<double>::Element d;

    


  if(0==Cptr->rank()){
////////////////////////////////////////////////////////////////////////////////////////////////
srand (time(NULL));
    for (long i = 0; i < ni; ++i)
      for (long j = 0; j < nj; ++j){
A.setEntry(i,j,get_random_float());
      }
 
    for (long j = 0; j < nj; ++j){
B.setEntry(j,get_random_float());
    }

//////////////////////////////////////////////////////////////////////////////////////////////////

//B.write(std::cout << ">>>>Compute with B:\n",Tag::FileFormat::Maple) << ';' << std::endl;

//A.write(std::cout << ">>>>Compute with A:\n",Tag::FileFormat::Maple) << ';' << std::endl;


  }//End of BLock for process(0)



////////////////////////////////////////////////////////////////////////////////////////////////////

if(0==Cptr->rank()){

//double starttime, endtime; 
//MPI_Barrier(MPI_COMM_WORLD);
//starttime = MPI_Wtime();
  //MPI data distribution for Integer type value

    Cptr->send(B,1);   // MPIgmpBcast(A, ni, nj, 0, Cptr);

//MPI_Barrier(MPI_COMM_WORLD);
//endtime   = MPI_Wtime(); 
//std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;

}else{

    Cptr->recv(B2,0);   // MPIgmpBcast(A, ni, nj, 0, Cptr);

}

if(0==Cptr->rank()){

//double starttime, endtime; 
//MPI_Barrier(MPI_COMM_WORLD);
//starttime = MPI_Wtime();
  //MPI data distribution for Integer type value

    Cptr->send(A,1);   // MPIgmpBcast(B, ni, 0, Cptr);
//MPI_Barrier(MPI_COMM_WORLD);
//endtime   = MPI_Wtime(); 
//std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;

}else{

    Cptr->recv(A2,0);   // MPIgmpBcast(B, ni, 0, Cptr
}



Cptr->bcast(B,0);
Cptr->bcast(A,0);

 if(0!=Cptr->rank()) checkResult (ZZ, A, A2);
 if(0!=Cptr->rank()) checkResult (ZZ, B, B2);

//Check if data are correctly distributed to all processes

  if(0!=Cptr->rank()){
/*
B.write(std::cout << "process("<<Cptr->rank()<< ")<<<<Compute with B:\n",Tag::FileFormat::Maple) << ';' << std::endl;
 A.write(std::cout << "process("<<Cptr->rank()<< ")<<<<Compute with A: \n",Tag::FileFormat::Maple) << ';' << std::endl;
*/
  
}
  MPI_Finalize();
  return 0;  
////////////////////////////////////////////////////////////////////////////////////////////////////


  
}


