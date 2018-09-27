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
#include <givaro/modular.h>
#include <givaro/zring.h>
#include <iostream>
#include <linbox/linbox-config.h>
#include <linbox/matrix/sparse-matrix.h>
#include <stdio.h>
#include <stdlib.h>

#include "linbox/util/mpicpp.h"
#include <mpi.h>
//#include "linbox/util/mpi-gmp.inl"

#include "linbox/blackbox/random-matrix.h"
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

template <class Field> static bool checkResult(const Field& ZZ, BlasMatrix<Field>& A, BlasMatrix<Field>& A2)
{
  // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  for (long i = 0; i < A.rowdim(); ++i)
    for (long j = 0; j < A.coldim(); ++j) {
      if (!ZZ.areEqual(A[i][j], A2[i][j])) {
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }
    }
  
  return true;
}


template <class Field>
static bool checkResult(const Field& ZZ, SparseMatrix<Givaro::ZRing<Integer>>& A, SparseMatrix<Givaro::ZRing<Integer>>& A2)
{
  // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  for (long i = 0; i < A.rowdim(); ++i)
    for (long j = 0; j < A.coldim(); ++j) {
      if (!ZZ.areEqual(A.getEntry(i, j), A2.getEntry(i, j))) {
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }
    }
  
  return true;
}

template <class Field, class T>
static bool checkResult(const Field& ZZ, SparseMatrix<Givaro::ZRing<T>>& A, SparseMatrix<Givaro::ZRing<T>>& A2)
{
  // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  
  for (long i = 0; i < A.rowdim(); ++i)
    for (long j = 0; j < A.coldim(); ++j) {
      if (!ZZ.areEqual(A.getEntry(i, j), A2.getEntry(i, j))) {
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }
    }
  
  return true;
}


template <class Field>
static bool checkResult(const Field& ZZ, SparseMatrix<Givaro::Modular<Integer>>& A, SparseMatrix<Givaro::Modular<Integer>>& A2)
{
  // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  for (long i = 0; i < A.rowdim(); ++i)
    for (long j = 0; j < A.coldim(); ++j) {
      if (!ZZ.areEqual(A.getEntry(i, j), A2.getEntry(i, j))) {
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }
    }
  
  return true;
}

template <class Field, class T>
static bool checkResult(const Field& ZZ, SparseMatrix<Givaro::Modular<T>>& A, SparseMatrix<Givaro::Modular<T>>& A2)
{
  // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  
  for (long i = 0; i < A.rowdim(); ++i)
    for (long j = 0; j < A.coldim(); ++j) {
      if (!ZZ.areEqual(A.getEntry(i, j), A2.getEntry(i, j))) {
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }
    }
  
  return true;
}

/*
  template <class Field>
  static bool checkResult(const Field& ZZ, SparseMatrix<Givaro::ModularBalanced<Integer>>& A, SparseMatrix<Givaro::ModularBalanced<Integer>>& A2)
  {
  // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  for (long i = 0; i < A.rowdim(); ++i)
  for (long j = 0; j < A.coldim(); ++j) {
  if (!ZZ.areEqual(A.getEntry(i, j), A2.getEntry(i, j))) {
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "               The data communicated is inconsistent                " << std::endl;
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  return false;
  }
  }
  
  return true;
  }
*/
template <class Field, class T>
static bool checkResult(const Field& ZZ, SparseMatrix<Givaro::ModularBalanced<T>>& A, SparseMatrix<Givaro::ModularBalanced<T>>& A2)
{
  // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  
  for (long i = 0; i < A.rowdim(); ++i)
    for (long j = 0; j < A.coldim(); ++j) {
      if (!ZZ.areEqual(A.getEntry(i, j), A2.getEntry(i, j))) {
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The data communicated is inconsistent                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	return false;
      }
    }
  
  return true;
}

template <class Field> static bool checkResult(const Field& ZZ, BlasVector<Field>& B, BlasVector<Field>& B2)
{
  // B.write(std::cout << " B: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  // B2.write(std::cout << " B2: \n",Tag::FileFormat::Maple) << ';' << std::endl;
  for (size_t j = 0; j < B.size(); ++j)
    if (!ZZ.areEqual(B[j], B2[j])) {
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The data communicated is inconsistent                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      return false;
    }
  return true;
}

template <class Field> bool genData(Givaro::Integer q, BlasMatrix<Field>& A, size_t bits)
{
  
  Field ZZ;//typename Field::Element ZZ;
  typedef typename Field::RandIter RandIter;
  RandIter RI(ZZ);
  LinBox::RandomDenseMatrix<RandIter, Field> RDM(ZZ, RI);
  RDM.randomFullRank(A);
  
}

template <class Field> bool genData(Givaro::Integer q, SparseMatrix<Field>& A, size_t bits)
{
  Field ZZ;//typename Field::Element ZZ;
  typedef typename Field::RandIter RandIter;
  RandIter RI(ZZ);
  LinBox::RandomDenseMatrix<RandIter, Field> RDM(ZZ, RI);
  RDM.random(A);
  
}

template <class Field> bool genData(Givaro::Integer q, BlasVector<Field>& B, size_t bits)
{
  
  Field ZZ;//typename Field::Element ZZ;
  typedef typename Field::RandIter RandIter;
  RandIter RI(ZZ);
  B.random(RI);
  
}

template <> bool genData(Givaro::Integer q, BlasMatrix<Givaro::ZRing<Integer>>& A, size_t bits)
{
  Givaro::ZRing<Integer> ZZ;
  typedef typename Givaro::ZRing<Integer>::RandIter RandIter;
  RandIter RI(ZZ, bits);
  LinBox::RandomDenseMatrix<RandIter, Givaro::ZRing<Integer>> RDM(ZZ, RI);
  RDM.randomFullRank(A);
}
template <> bool genData(Givaro::Integer q, SparseMatrix<Givaro::ZRing<Integer>>& A, size_t bits)
{
  Givaro::ZRing<Integer> ZZ;
  typedef typename Givaro::ZRing<Integer>::RandIter RandIter;
  RandIter RI(ZZ, bits);
  LinBox::RandomDenseMatrix<RandIter, Givaro::ZRing<Integer>> RDM(ZZ, RI);
  RDM.random(A);
  
}
template <> bool genData(Givaro::Integer q, BlasVector<Givaro::ZRing<Integer>>& B, size_t bits)
{
  Givaro::ZRing<Integer> ZZ;
  typedef typename Givaro::ZRing<Integer>::RandIter RandIter;
  RandIter RI(ZZ, bits);
  B.random(RI);
}


template <class T> 
bool genData(Givaro::Integer q, BlasMatrix<Givaro::Modular<T>>& A, size_t bits)
{
  
  Givaro::Modular<T> ZZ(q);
  typedef typename Givaro::Modular<T>::RandIter RandIter;
  RandIter RI(ZZ, bits);
  LinBox::RandomDenseMatrix<RandIter, Givaro::Modular<T>> RDM(ZZ, RI);
  RDM.random(A);
  
}
template <class T> 
bool genData(Givaro::Integer q, SparseMatrix<Givaro::Modular<T>>& A, size_t bits)
{
  
  Givaro::Modular<T> ZZ(q);
  typedef typename Givaro::Modular<T>::RandIter RandIter;
  RandIter RI(ZZ, bits);
  LinBox::RandomDenseMatrix<RandIter, Givaro::Modular<T>> RDM(ZZ, RI);
  RDM.random(A);
  
}
template <class T> 
bool genData(Givaro::Integer q, BlasVector<Givaro::Modular<T>>& B, size_t bits)
{
  Givaro::Modular<T> ZZ(q);//typename Field::Element ZZ;
  typedef typename Givaro::Modular<T>::RandIter RandIter;
  RandIter RI(ZZ);
  B.random(RI);
  
}


template <class T> 
bool genData(Givaro::Integer q, BlasMatrix<Givaro::ModularBalanced<T>>& A, size_t bits)
{
  
  Givaro::ModularBalanced<T> ZZ(q);
  typedef typename Givaro::ModularBalanced<T>::RandIter RandIter;
  RandIter RI(ZZ, bits);
  LinBox::RandomDenseMatrix<RandIter, Givaro::ModularBalanced<T>> RDM(ZZ, RI);
  RDM.random(A);
  
}
template <class T> 
bool genData(Givaro::Integer q, SparseMatrix<Givaro::ModularBalanced<T>>& A, size_t bits)
{
  
  Givaro::ModularBalanced<T> ZZ(q);
  typedef typename Givaro::ModularBalanced<T>::RandIter RandIter;
  RandIter RI(ZZ, bits);
  LinBox::RandomDenseMatrix<RandIter, Givaro::ModularBalanced<T>> RDM(ZZ, RI);
  RDM.random(A);
  
}
template <class T> 
bool genData(Givaro::Integer q, BlasVector<Givaro::ModularBalanced<T>>& B, size_t bits)
{
  Givaro::ModularBalanced<T> ZZ(q);//typename Field::Element ZZ;
  typedef typename Givaro::ModularBalanced<T>::RandIter RandIter;
  RandIter RI(ZZ);
  B.random(RI);
  
}


template <class Field> 
void test_with_field(Givaro::Integer q, size_t bits, size_t ni, size_t nj, Communicator* Cptr)
{
  
  Field ZZ(q);
  
  typedef BlasVector<Field> DenseVector;
  
  DenseMatrix<Field> A(ZZ, ni, nj), A2(ZZ, ni, nj);
  DenseVector B2(ZZ, A.coldim()), B(ZZ, A.coldim());
  SparseMatrix<Field> A3(ZZ, ni, nj), A4(ZZ, ni, nj);
  
  //DenseMatrix<Field> A5(ZZ, ni, nj), A6(ZZ, ni, nj);
  //SparseMatrix<Field> A7(ZZ, ni, nj), A8(ZZ, ni, nj);
  
  if (0 == Cptr->rank()) {
    
    genData(q, A, bits);
    genData(q, A3, bits);
    //genData(q, A5, bits);
    //genData(q, A7, bits);
    genData(q, B, bits);
    
  } // End of BLock for process(0)
  
#if 1
  if (0 == Cptr->rank()) {
    
    // double starttime, endtime;
    // MPI_Barrier(MPI_COMM_WORLD);
    // starttime = MPI_Wtime();
    // MPI data distribution for Integer type value
    Cptr->ssend(B, 1);
    // MPI_Barrier(MPI_COMM_WORLD);
    // endtime   = MPI_Wtime();
    // std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
  }
  else {
    
    if (1 == Cptr->rank()) Cptr->recv(B2, 0);
    
  }
  Cptr->bcast(B, 0);
  if (1 == Cptr->rank()) checkResult(ZZ, B, B2);
  
  if (0 == Cptr->rank()) {

    // double starttime, endtime;
    // MPI_Barrier(MPI_COMM_WORLD);
    // starttime = MPI_Wtime();
    // MPI data distribution for Integer type value
    Cptr->ssend(A, 1);
    // MPI_Barrier(MPI_COMM_WORLD);
    // endtime   = MPI_Wtime();
    // std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
  }
  else {
    
    if (1 == Cptr->rank()) Cptr->recv(A2, 0);
    
  }
  Cptr->bcast(A, 0);
  
  if (1 == Cptr->rank()) checkResult(ZZ, A, A2);
  
  if (0 == Cptr->rank()) {
    
    // double starttime, endtime;
    // MPI_Barrier(MPI_COMM_WORLD);
    // starttime = MPI_Wtime();
    // MPI data distribution for Integer type value
    Cptr->ssend(A3, 1);
    // MPI_Barrier(MPI_COMM_WORLD);
    // endtime   = MPI_Wtime();
    // std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
  }
  else {
    
    if (1 == Cptr->rank()) Cptr->recv(A4, 0);
    
  }
  Cptr->bcast(A3, 0);
  
  if (1 == Cptr->rank()) checkResult(ZZ, A3, A4);
#endif
  
#if 0
  if (0 == Cptr->rank()) {
    // double starttime, endtime;
    // MPI_Barrier(MPI_COMM_WORLD);
    // starttime = MPI_Wtime();
    // MPI data distribution for Integer type value
    Cptr->isend(A5, 1);
    // MPI_Barrier(MPI_COMM_WORLD);
    // endtime   = MPI_Wtime();
    // std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
  }
  else {
    
    if (1 == Cptr->rank()) Cptr->recv(A6, 0);
    
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Cptr->bcast(A5, 0);
  
  if (1 == Cptr->rank()) checkResult(ZZ, A5, A6);
  
  if (0 == Cptr->rank()) {
    // double starttime, endtime;
    // MPI_Barrier(MPI_COMM_WORLD);
    // starttime = MPI_Wtime();
    // MPI data distribution for Integer type value
    Cptr->isend(A7, 1);
    // MPI_Barrier(MPI_COMM_WORLD);
    // endtime   = MPI_Wtime();
    // std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
  }
  else {
    if (1 == Cptr->rank()) Cptr->recv(A8, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Cptr->bcast(A7, 0);
  
  if (1 == Cptr->rank()) checkResult(ZZ, A7, A8);

#endif
}

int main(int argc, char** argv)
{
  
  Communicator* Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
  size_t bits, ni, niter, nj, n;
  bool loop=false;
  
  bits = 10, niter = 1, ni = 3, nj = 3;
  Givaro::Integer q=101;
  static Argument args[] = {
    {'b', "-b B", "Set the mxaimum number of digits of integers to generate.", TYPE_INT, &bits},
    {'i', "-i I", "Set the number of iteration over unit test sets.", TYPE_INT, &niter},
    { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
    {'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT, &n},
    { 'q', "-q Q", "Set the field cardinality.",  TYPE_INTEGER , &q },
    END_OF_ARGUMENTS};
  parseArguments(argc, argv, args);
  
  
  bool peak=false;
  srand(time(NULL));
  do{
    for (int j = 0; j < niter; j++) {
      if (0 == Cptr->rank()) {
	ni=rand() % n + 1; 
	if (ni<n/2 && ni%2==0 && !peak) ni=1;
	
	peak=!peak;
	std::cout<<" Test with dimension: " << ni <<" x "<< ni<<std::endl;
      }
      MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD);
      nj=ni;
      MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      test_with_field<Givaro::ZRing<Integer>>            (q, bits, ni, nj, Cptr);
      
      test_with_field<Givaro::Modular<float> >           (q, bits, ni, nj, Cptr);
      test_with_field<Givaro::Modular<double> >          (q, bits, ni, nj, Cptr);
      test_with_field<Givaro::Modular<int32_t> >   	   (q, bits, ni, nj, Cptr);
      
      
      test_with_field<Givaro::ZRing<float>>		   (q, bits, ni, nj, Cptr);
      test_with_field<Givaro::ZRing<double>>             (q, bits, ni, nj, Cptr);
      test_with_field<Givaro::ZRing<int32_t>>            (q, bits, ni, nj, Cptr);
      
      test_with_field<Givaro::ModularBalanced<float> >    (q, bits, ni, nj, Cptr);
      
      test_with_field<Givaro::ModularBalanced<double> >   (q, bits, ni, nj, Cptr);
      test_with_field<Givaro::ModularBalanced<int32_t> >  (q, bits, ni, nj, Cptr);
      
      
      //test_with_field<Givaro::Modular<int64_t> >   	   (q, bits, ni, nj, Cptr);
      //test_with_field<Givaro::ModularBalanced<int64_t> >  (q, bits, ni, nj, Cptr);
      
      
      test_with_field<Givaro::Modular<Givaro::Integer> >  (q, bits, ni, nj, Cptr);
      
      
    }
  }while(loop);
  MPI_Finalize();
  return 0;
}
