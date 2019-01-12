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
#endif
#include <linbox/solutions/methods.h>
#include <linbox/solutions/solve.h>

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

std::string gmp_rand ( long maxNdigits)
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
#ifdef __LINBOX_HAVE_MPI

template <typename T > class chooseMPItype;
template <> struct chooseMPItype<unsigned int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED;};
template <> struct chooseMPItype<unsigned long long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG_LONG;};
template <> struct chooseMPItype<unsigned long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG;};


void printMPItype(MPI_Datatype type){ 
  if(MPI_UNSIGNED==type)  std::cout<<"MPI_UNSIGNED is to be used"<< std::endl;
  if(MPI_UNSIGNED_LONG_LONG==type)   std::cout<<"MPI_UNSIGNED_LONG_LONG is to be used"<< std::endl;
  if(MPI_UNSIGNED_LONG==type) std::cout<<"MPI_UNSIGNED_LONG is to be used"<<std::endl;
};




template <class LinMat>
inline void MPIgmpBcast(LinMat& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
{ 
  int A_mp_alloc[ni*nj], A_a_size[ni*nj];
  unsigned lenA;
  std::vector<mp_limb_t> A_mp_data;
  //MPI_Barrier(MPI_COMM_WORLD);
  if(src==Cptr->rank()){
    
    //Split Matrix A into arrays 
    //std::cerr << "A=:= " << std::endl; 
    __mpz_struct * ptr;
    for (size_t i = 0; i < ni; ++i){
      for (size_t j = 0; j < nj; ++j){
	
	//std::cerr << A.getEntry(i,j)<< "\t" ; std::cerr<< std::endl;
	ptr = const_cast<__mpz_struct*>(A.getEntry(i,j).get_mpz());
	A_mp_alloc[j+i*nj] = ptr->_mp_alloc;
	A_a_size[j+i*nj] = ptr->_mp_size;
	mp_limb_t * a_array = ptr->_mp_d;
	for(long k=0; k< ptr->_mp_alloc; ++k)
	  A_mp_data.push_back(a_array[k]);
	
	
      }
    }
    lenA = A_mp_data.size();
  }
  //Distribut Givaro::Integer through its elementary parts
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&A_mp_alloc[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A_a_size[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&lenA, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if(src!=Cptr->rank()) A_mp_data.resize(lenA);
  MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  if(src!=Cptr->rank()){
    //Reconstruction of matrix A
    __mpz_struct * ptr2;
    size_t count=0;
    Givaro::Integer temp; 
    //std::cerr << "received A:= " << std::endl;
    for (size_t i = 0; i < ni; ++i){
      for (size_t j = 0; j < nj; ++j){

	ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
	ptr2->_mp_alloc = A_mp_alloc[j+i*nj];
	ptr2->_mp_size = A_a_size[j+i*nj];
	_mpz_realloc(ptr2,ptr2->_mp_alloc);
	for(long k=0; k< ptr2->_mp_alloc; ++k){
	  ptr2->_mp_d[k] = (A_mp_data[k+count]);	
	}count+=ptr2->_mp_alloc;
	A.setEntry(i,j,temp);
	//std::cerr << A.getEntry(i,j) << "\t" ; std::cerr<< std::endl;
      }
    }
    
  }//MPI_Barrier(MPI_COMM_WORLD);
  
}


//template<typename T>
//void MPIvecBcast(LinBox::SparseMatrix< T > & A, int ni, int nj, Communicator *Cptr)
template <class LinMat>
void MPIgmpRMA(LinMat& A, int ni, int nj, Communicator *Cptr)
{ 
  int A_mp_alloc[ni*nj], A_a_size[ni*nj];
  int lenA;
  mp_limb_t* _mp_data; 

  if(0==Cptr->rank()){

    std::vector<mp_limb_t> A_mp_data;
    //Split Matrix A into arrays 
    //std::cerr << "A=:= " << std::endl; 
    __mpz_struct * ptr;
    for (long i = 0; i < ni; ++i){
      for (long j = 0; j < nj; ++j){	
	//std::cerr << A.getEntry(i,j)<< "\t" ; std::cerr<< std::endl;
	ptr = const_cast<__mpz_struct*>(A.getEntry(i,j).get_mpz());
	A_mp_alloc[j+i*nj] = ptr->_mp_alloc;
	A_a_size[j+i*nj] = ptr->_mp_size;
//std::cerr<<"A_mp_alloc["<<j+i*nj<< "]:="<<A_mp_alloc[j+i*nj]<<"   "<<"A_a_size["<<j+i*nj<<"]:="<<A_a_size[j+i*nj]<<std::endl;
	mp_limb_t * a_array = ptr->_mp_d;
	for(size_t k=0; k< ptr->_mp_alloc; ++k){
//std::cerr<<"A_mp_data["<<k<<"]="<<a_array[k]<< std::endl;
	  A_mp_data.push_back(a_array[k]);
        }
      }
    }
lenA = A_mp_data.size(); std::cerr<<"lenA:="<<lenA<< std::endl;
_mp_data = (mp_limb_t*) malloc(lenA * sizeof(mp_limb_t)); 
    
for(size_t k=0; k< lenA; ++k)  _mp_data[k]=A_mp_data[k];
  }MPI_Barrier(MPI_COMM_WORLD); //All processes must wait until the initialization finishes

MPI_Win win;
if(Cptr->rank() == 0) MPI_Win_create(A_mp_alloc, ni*nj, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
else MPI_Win_create(NULL,0,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&win); 
MPI_Win_fence(MPI_MODE_NOPRECEDE,win);
if(0!=Cptr->rank()){
MPI_Get(A_mp_alloc, ni*nj, MPI_INT, 0, 0, ni*nj, MPI_INT, win);
}
MPI_Win_fence(0,win);

if(Cptr->rank() == 0) MPI_Win_create(A_a_size, ni*nj, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win); 
else MPI_Win_create(NULL,0,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&win); 
MPI_Win_fence(MPI_MODE_NOPRECEDE,win);
if(0!=Cptr->rank()){
MPI_Get(A_a_size, ni*nj, MPI_INT, 0, 0, ni*nj, MPI_INT, win);
}
MPI_Win_fence(0,win);

if(Cptr->rank() == 0) MPI_Win_create(&lenA, 1, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
else MPI_Win_create(NULL,0,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&win); 
MPI_Win_fence(MPI_MODE_NOPRECEDE,win);
if(0!=Cptr->rank()){
MPI_Get(&lenA, 1, MPI_INT, 0, 0, 1, MPI_INT, win);
}MPI_Win_fence(0,win);

if(Cptr->rank() == 0) MPI_Win_create(_mp_data, lenA, sizeof(mp_limb_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
else MPI_Win_create(NULL,0,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&win); 
MPI_Win_fence(MPI_MODE_NOPRECEDE,win);
if(0!=Cptr->rank()){
_mp_data=(mp_limb_t*) malloc(lenA * sizeof(mp_limb_t)); //std::cerr<<"********** lenB:="<<lenB<< std::endl;
MPI_Get(_mp_data, lenA, chooseMPItype<mp_limb_t>::val, 0, 0, lenA, chooseMPItype<mp_limb_t>::val, win);
}
MPI_Win_fence(0,win);

MPI_Barrier(MPI_COMM_WORLD);
if(0!=Cptr->rank()){
    //Reconstruction of matrix A
    __mpz_struct * ptr2;
    size_t count=0;
    Givaro::Integer temp; 
    //std::cerr << "received A:= " << std::endl;
    for (long i = 0; i < ni; ++i){
      for (long j = 0; j < nj; ++j){
	ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
	ptr2->_mp_alloc = A_mp_alloc[j+i*nj];
	ptr2->_mp_size = A_a_size[j+i*nj];
	_mpz_realloc(ptr2,ptr2->_mp_alloc);
	for(size_t k=0; k< ptr2->_mp_alloc; ++k){
	  ptr2->_mp_d[k] = (_mp_data[k+count]);
	}count+=ptr2->_mp_alloc;
	A.setEntry(i,j,temp);
	//std::cerr << A.getEntry(i,j) << "\t" ; std::cerr<< std::endl;
      }
    }
}

MPI_Win_free(&win);
  
}
#endif



//#include <map>
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char ** argv)
{
  
  
#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
  long dg,ni,nj;  
  long tmp;

   dg=10, ni=11,nj=11;

 	static Argument args[] = {
		{ 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
		{ 'd', "-d M", "Set the mxaimum number of digits of integers to generate.", TYPE_INT,     &dg },
		END_OF_ARGUMENTS
	};	
parseArguments (argc, argv, args); nj = ni;
MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); nj=ni;


//std::cerr << "proc("<< Cptr->rank() <<") init>> ni = nj = "<< nj << std::endl;

  Givaro::ZRing<Integer> ZZ;
  SparseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);

  typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;
  DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
  Givaro::ZRing<Integer>::Element d;



  if(0==Cptr->rank()){

const char * c;
std::string result;
long unsigned int r=0;
/*
  ofstream myfile;  ofstream myfile2;
  myfile.open ("matrix.txt"); myfile2.open ("vec.txt");
  myfile <<ni<<" "<<nj<<" M\n"; myfile2 <<nj<< " 1"<<" M\n";
*/
while(r!=ni){
    for (long i = 0; i < ni; ++i)
      for (long j = 0; j < nj; ++j){
result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
A.setEntry(i,j,res);
result.clear();
//myfile <<i+1<<" "<<j+1<<" "<<A.getEntry(i,j)<<"\n";
      }
/*  
    for (long j = 0; j < nj; ++j){
result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
B.setEntry(j,res);
result.clear();
myfile2 <<j+1<< " 1 "<<B.getEntry(j)<<"\n";
    }
*/
LinBox::rank (r, A);

}
/*
myfile <<"0  0  0";
myfile2 <<"0  0  0";
 myfile.close();
 myfile2.close();
*/
LinBox::rank (r, A); //std::cout<<"The rank of generated matrix A is:"<<r<<std::endl;  
 std::cout << "The rank of generated matrix A is: "<<r<<std::endl; 


/*
      std::cerr << "Compute with B: " << std::endl;
      for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
      
      std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
      if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
*/

  }//End of BLock for process(0)


/*
  //MPI data distribution for Integer type value
  MPIgmpBcast(A, ni, nj, Cptr);
  MPIgmpBcast(B, nj, Cptr);
*/


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/*
  int A_mp_alloc[ni*nj], A_a_size[ni*nj];
  int lenA;
  mp_limb_t* _mp_data; 

int B_mp_alloc[ni*nj], B_a_size[ni*nj];
int lenB;
mp_limb_t* B_mp_data;


MPI_Barrier(MPI_COMM_WORLD);
if(0!=Cptr->rank()){
std::cerr<<"proc("<<Cptr->rank() << ") lenB:="<<lenB<< std::endl;
    for (long i = 0; i < ni; ++i)
      for (long j = 0; j < nj; ++j){	
//std::cerr<<"B_mp_alloc["<<j+i*nj<< "]:="<<B_mp_alloc[j+i*nj]<<std::endl;
//std::cerr<<"B_mp_alloc["<<j+i*nj<< "]:="<<B_mp_alloc[j+i*nj]<<"   "<<"B_a_size["<<j+i*nj<< "]:="<<B_a_size[j+i*nj]<<std::endl;
      }std::cerr<<std::endl;
//for (long i = 0; i < lenB; ++i) std::cerr<<"B_mp_data["<<i<< "]:="<<B_mp_data[i]<<std::endl;
}*/

  double starttime, endtime; 
  MPI_Barrier(MPI_COMM_WORLD);
  starttime = MPI_Wtime(); 
  MPIgmpRMA(A, ni, nj, Cptr);
  MPI_Barrier(MPI_COMM_WORLD);	
  endtime   = MPI_Wtime(); 
  std::cout<<"MPIgmpRMA used CPU time: "<< endtime-starttime <<std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  starttime = MPI_Wtime(); 
  MPIgmpBcast(A, ni, nj, 0, Cptr);
  MPI_Barrier(MPI_COMM_WORLD);	
  endtime   = MPI_Wtime(); 
  std::cout<<"MPIgmpBcast used CPU time: "<< endtime-starttime <<std::endl;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


std::cout<<"Computation is done over Q"<<std::endl;
  /***********************
   Results verification 
  ***********************/
/*    
	std::cerr <<"proc("<<Cptr->rank()<<") Computes with B: " << std::endl;
	for(long j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 

	
 MPI_Barrier(MPI_COMM_WORLD);
if(0!=Cptr->rank()){
	std::cerr <<"proc("<<Cptr->rank()<<") Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
	if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A::=",Tag::FileFormat::Maple) << ';' << std::endl;
}
*/

 MPI_Barrier(MPI_COMM_WORLD);
  if(0==Cptr->rank()){    
   
  MPI_Finalize();
  return 0;    
  }
  MPI_Finalize(); 
  
    
#else ///////////////////////////////////SEQUENNTIAL//////////////////////////////////////////////////////
   std::cout<<"Computation is done over Q"<<std::endl;
  srand48(time(NULL));
  long dg=20, ni=3,nj=3,max=100;
 
  Givaro::ZRing<Integer> ZZ;
  DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);
  
  typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;
  DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
  Givaro::ZRing<Integer>::Element d;

const char * c;
std::string result;
long unsigned int r=0;
  ofstream myfile;  ofstream myfile2;
  myfile.open ("matrix.txt"); myfile2.open ("vec.txt");
  myfile <<ni<<" "<<nj<<" M\n"; myfile2 <<nj<< " 1"<<" M\n";
while(r!=ni){
    for (long i = 0; i < ni; ++i)
      for (long j = 0; j < nj; ++j){
result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
A.setEntry(i,j,res);
result.clear();
myfile <<i+1<<" "<<j+1<<" "<<A.getEntry(i,j)<<"\n";
      }
  
    for (long j = 0; j < nj; ++j){
result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
B.setEntry(j,res);
result.clear();
myfile2 <<j+1<< " 1 "<<B.getEntry(j)<<"\n";
    }

LinBox::rank (r, A);
}
myfile <<"0  0  0";
myfile2 <<"0  0  0";
 myfile.close();
 myfile2.close();

LinBox::rank (r, A); //std::cout<<"The rank of generated matrix A is:"<<r<<std::endl;  
 std::cout << "\033[1;33mThe rank of generated matrix A is:\033[0m"<<r<<std::endl;   
 /***********************
  Results verification 
 ***********************/
/*    
	std::cerr << "Compute with B: " << std::endl;
	for(long j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
	
	std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
	if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A::=",Tag::FileFormat::Maple) << ';' << std::endl;
*/


  return 0;
  
#endif
  
}


