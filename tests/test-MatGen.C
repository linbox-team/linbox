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

template <class Vector>
//void MPIgmpBcast(Vector& B, int ni, int nj, Communicator *Cptr)
void MPIgmpBcast(Vector& B, int nj, Communicator *Cptr)
{ 
int B_mp_alloc[nj], B_a_size[nj]; 
unsigned lenB;  Givaro::Integer temp; 
std::vector<mp_limb_t> B_mp_data;
MPI_Barrier(MPI_COMM_WORLD);
  if(0==Cptr->rank()){
//std::cerr << "B=:= " << std::endl;
	//Split vector B into arrays
	__mpz_struct * ptr;
	for(int j=0;j<nj;j++){
	//std::cerr << B.getEntry(j)<< "\t" ; std::cerr<< std::endl;
	ptr = const_cast<__mpz_struct*>(B.getEntry(j).get_mpz());
	B_mp_alloc[j] = ptr->_mp_alloc;
	B_a_size[j] = ptr->_mp_size;
	    mp_limb_t * a_array = ptr->_mp_d;
	    for(size_t i=0; i< ptr->_mp_alloc; ++i){
		B_mp_data.push_back(a_array[i]);
	    }
	}

	lenB = B_mp_data.size();

	}
//Distribut Givaro::Integer through its elementary parts
MPI_Barrier(MPI_COMM_WORLD);
MPI_Bcast(&B_mp_alloc[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&B_a_size[0], nj, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
if(0!=Cptr->rank()) B_mp_data.resize(lenB);
MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);

 if(0!=Cptr->rank()){
	//Reconstruction of vector B
	//std::cerr << "received B::= " << std::endl;
	__mpz_struct * ptr2;
	size_t count=0;

	for(int j=0;j<nj;j++){ 
	ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
	ptr2->_mp_alloc = B_mp_alloc[j];
	ptr2->_mp_size = B_a_size[j];
	    _mpz_realloc(ptr2,ptr2->_mp_alloc);
	    for(size_t i=0; i< ptr2->_mp_alloc; ++i){
		ptr2->_mp_d[i] = (B_mp_data[i+count]);
	    }count+=ptr2->_mp_alloc;
	B.setEntry(j,temp);
	//std::cerr << B.getEntry(j) << "\t" ; std::cerr<< std::endl; 

	}
	
	}MPI_Barrier(MPI_COMM_WORLD);

}

//template<typename T>
//void MPIvecBcast(LinBox::SparseMatrix< T > & A, int ni, int nj, Communicator *Cptr)
template <class LinMat>
void MPIgmpBcast(LinMat& A, int ni, int nj, Communicator *Cptr)
{ 
int A_mp_alloc[ni*nj], A_a_size[ni*nj];
unsigned lenA;
std::vector<mp_limb_t> A_mp_data;
MPI_Barrier(MPI_COMM_WORLD);
  if(0==Cptr->rank()){

//Split Matrix A into arrays 
//std::cerr << "A=:= " << std::endl; 
__mpz_struct * ptr;
   for (long i = 0; i < ni; ++i){
     for (long j = 0; j < nj; ++j){

//std::cerr << A.getEntry(i,j)<< "\t" ; std::cerr<< std::endl;
ptr = const_cast<__mpz_struct*>(A.getEntry(i,j).get_mpz());
A_mp_alloc[j+i*nj] = ptr->_mp_alloc;
A_a_size[j+i*nj] = ptr->_mp_size;
    mp_limb_t * a_array = ptr->_mp_d;
    for(size_t k=0; k< ptr->_mp_alloc; ++k)
        A_mp_data.push_back(a_array[k]);


}
}
lenA = A_mp_data.size();
}
//Distribut Givaro::Integer through its elementary parts
MPI_Barrier(MPI_COMM_WORLD);
MPI_Bcast(&A_mp_alloc[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&A_a_size[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&lenA, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
if(0!=Cptr->rank()) A_mp_data.resize(lenA);
MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
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
        ptr2->_mp_d[k] = (A_mp_data[k+count]);	
    }count+=ptr2->_mp_alloc;
A.setEntry(i,j,temp);
//std::cerr << A.getEntry(i,j) << "\t" ; std::cerr<< std::endl;
}
}

}MPI_Barrier(MPI_COMM_WORLD);

}
#endif



#include <gmp++/gmp++.h>
#include <string>


int main(int argc, char ** argv)
{
  
std::cout<<"Computation is done over Q"<<std::endl;
#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);MPI_Barrier(MPI_COMM_WORLD);
  int ni=21,nj=21,max=10000;  
  Givaro::ZRing<Integer> ZZ;
  SparseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);

typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;
DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
Givaro::ZRing<Integer>::Element d;




  if(0==Cptr->rank()){
   srand48(time(NULL));
   long tmp;
Givaro::Integer valeur("123456789123456789123456789123456789123456789123456789123456789123456789"); 
   printf("%ld %ld M\n",1, ni, nj);
   for (long i = 0; i < ni; ++i)
     for (long j = 0; j < nj; ++j){
       A.setEntry(i,j,myrand(tmp, max));
     }
//A.setEntry(nj-1,nj-1,valeur);
//   std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
//   if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cerr << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
Givaro::Integer valeur2("448832189123456789123456789123456789123456789123456789123456789123456789"); 
   for (long j = 0; j < nj-1; ++j)
	B.setEntry(j,myrand(tmp, max));

  }//End of BLock for process(0)


MPIgmpBcast(A, ni, nj, Cptr);
MPIgmpBcast(B, nj, Cptr);


/***********************
 Results verification 
***********************/
RingCategories::IntegerTag tg;
Timer chrono;
   double starttime, endtime; 
   starttime = MPI_Wtime(); 
solveCRA (X2, d, A, B, tg, Method::Hybrid(*Cptr),Cptr);	
   endtime   = MPI_Wtime(); 
	

if(0==Cptr->rank()){


std::cerr << "Compute with B: " << std::endl;
for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 

std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
    if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A::=",Tag::FileFormat::Maple) << ';' << std::endl;



		// BlasElimination
		//DenseVector X(ZZ, A.coldim());
                std::cout << "BlasElimination" << std::endl;
                chrono.start();
                solve (X, d, A, B, Method::BlasElimination());
                chrono.stop();

 		std::cout << "(BlasElimination) Solution is [";
                for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
 			ZZ.write(cout, *it) << " ";
                std::cout << "] / ";
                ZZ.write(std::cout, d)<< std::endl;
                std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;


		// solveCRA
		std::cout << "solveCRA" << std::endl;

			std::cout << "MPI CRA Solution is [";
			for(DenseVector::const_iterator it=X2.begin();it != X2.end(); ++it)
				ZZ.write(cout, *it) << " ";
			std::cout << "] / ";
			ZZ.write(std::cout, d) << std::endl;
			std::cout << "CPU time (seconds): " << endtime-starttime << std::endl;


	for (size_t j = 0 ; j < nj ; ++j){
		if(!Givaro::ZRing<Integer>().areEqual(X[j],X2[j])){
			std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			std::cerr << "                  ***Failed***                   " << std::endl;
			std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
			break;
		}
	}


}


MPI_Barrier(MPI_COMM_WORLD);
MPI_Finalize();
return 0;
  
  
#else

  
  srand48(time(NULL));
  long ni=3,nj=1,max=100;
  //	int offset = 0;
    
  Givaro::ZRing<Integer> ZZ;
  SparseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);


  long tmp;
  printf("%ld %ld M\n", ni, nj);
  for (long i = 0; i < ni; ++i)
    for (long j = 0; j < nj; ++j){
      A.setEntry(i,j,myrand(tmp, max));
      std::cout<<i+1<<" "<<j+1<<" "<<A.getEntry(i,j)<<std::endl;
    }
  
  std::cout<<"*****************************************************************"<<std::endl;
  
  std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
  if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cerr << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
  

  return 0;
  
#endif
  
}


