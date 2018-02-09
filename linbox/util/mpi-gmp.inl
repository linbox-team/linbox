#include "mpicpp.h"
#include <gmp++/gmp++.h>

template <typename T > class chooseMPItype;
template <> struct chooseMPItype<unsigned int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED;};
template <> struct chooseMPItype<unsigned long long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG_LONG;};
template <> struct chooseMPItype<unsigned long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG;};

void printMPItype(MPI_Datatype type){ 
  if(MPI_UNSIGNED==type)  std::cout<<"MPI_UNSIGNED is to be used"<< std::endl;
  if(MPI_UNSIGNED_LONG_LONG==type)   std::cout<<"MPI_UNSIGNED_LONG_LONG is to be used"<< std::endl;
  if(MPI_UNSIGNED_LONG==type) std::cout<<"MPI_UNSIGNED_LONG is to be used"<<std::endl;
};


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
template<class Vector>
void MPIintBcast(Vector& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr){
long long*A_data=(long long*)malloc(ni*nj*sizeof(long long));

  if(src==Cptr->rank()){

    for (size_t i = 0; i < ni; ++i)
      for (size_t j = 0; j < nj; ++j){
	A_data[j+i*nj] = A.getEntry(i,j);
      }
  }
  MPI_Bcast(&A_data[0], ni*nj, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  if(src!=Cptr->rank()){
    

    for (size_t i = 0; i < ni; ++i)
      for (size_t j = 0; j < nj; ++j){
	A.setEntry(i,j,A_data[j+i*nj]);
      }
  }

}


template<class Vector>
void MPIintBcast(Vector& A, size_t ni, int src, LinBox::Communicator *Cptr){
  long long*A_data=(long long*)malloc(ni*sizeof(long long));

  if(src==Cptr->rank()) for (size_t i = 0; i < ni; ++i) A_data[i] = A.getEntry(i);
  MPI_Bcast(&A_data[0], ni, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  if(src!=Cptr->rank())  for (size_t i = 0; i < ni; ++i) A.setEntry(i,A_data[i]);

}


template<class Vector, typename T>
struct gmpBcast2
{
 void operator()(Vector& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
 { }
};
template<class Vector>
struct gmpBcast2<Vector, long>
{
 void operator()(Vector& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
 {

MPIintBcast(A, ni, nj, src, Cptr);

 }
};
template<class Vector>
struct gmpBcast2<Vector, int>
{
 void operator()(Vector& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
 {

MPIintBcast(A, ni, nj, src, Cptr);

 }
};
template<class Vector>
struct gmpBcast2<Vector, LinBox::Integer>
{
 void operator()(Vector& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
 {
  //int A_mp_alloc[ni*nj], A_a_size[ni*nj];
int *A_mp_alloc=(int*)malloc(ni*nj*sizeof(int));
int *A_a_size=(int*)malloc(ni*nj*sizeof(int));
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
};
//template<class Vector, typename T>struct MPIgmpBcast;

template<class Vector, typename T>
struct gmpBcast
{
 void operator()(Vector& B, size_t ni, int src, LinBox::Communicator *Cptr)
 { }
};
template<class Vector>
struct gmpBcast<Vector, long>
{
 void operator()(Vector& B, size_t ni, int src, LinBox::Communicator *Cptr)
 {
MPIintBcast(B, ni, src, Cptr);
 }
};
template<class Vector>
struct gmpBcast<Vector, int>
{
 void operator()(Vector& B, size_t ni, int src, LinBox::Communicator *Cptr)
 { 
MPIintBcast(B, ni, src, Cptr);
 }
};
template<class Vector>
struct gmpBcast<Vector, LinBox::Integer>
{
 void operator()(Vector& B, size_t nj, int src, LinBox::Communicator *Cptr)
 {
//  int B_mp_alloc[nj], B_a_size[nj]; 
int *B_mp_alloc=(int*)malloc(nj*sizeof(int));
int *B_a_size=(int*)malloc(nj*sizeof(int));
  unsigned lenB;  Givaro::Integer temp; 
  std::vector<mp_limb_t> B_mp_data;
  //MPI_Barrier(MPI_COMM_WORLD);
  if(src==Cptr->rank()){
    //std::cerr << "B=:= " << std::endl;
    //Split vector B into arrays
    __mpz_struct * ptr;
    for(size_t j=0;j<nj;j++){
      //std::cerr << B.getEntry(j)<< "\t" ; std::cerr<< std::endl;
      ptr = const_cast<__mpz_struct*>(B.getEntry(j).get_mpz());
      B_mp_alloc[j] = ptr->_mp_alloc;
      B_a_size[j] = ptr->_mp_size;
      mp_limb_t * a_array = ptr->_mp_d;
      for(long i=0; i< ptr->_mp_alloc; ++i){
	B_mp_data.push_back(a_array[i]);
      }
    }
    
    lenB = B_mp_data.size();
    
  }
  //Distribut Givaro::Integer through its elementary parts
 // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&B_mp_alloc[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&B_a_size[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if(src!=Cptr->rank()) B_mp_data.resize(lenB);
  MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  if(src!=Cptr->rank()){
    //Reconstruction of vector B
    //std::cerr << "received B::= " << std::endl;
    __mpz_struct * ptr2;
    size_t count=0;
    
    for(size_t j=0;j<nj;j++){ 
      ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
      ptr2->_mp_alloc = B_mp_alloc[j];
      ptr2->_mp_size = B_a_size[j];
      _mpz_realloc(ptr2,ptr2->_mp_alloc);
      for(long i=0; i< ptr2->_mp_alloc; ++i){
	ptr2->_mp_d[i] = (B_mp_data[i+count]);
      }count+=ptr2->_mp_alloc;
      B.setEntry(j,temp);
      //std::cerr << B.getEntry(j) << "\t" ; std::cerr<< std::endl; 
      
    }
    
  }//MPI_Barrier(MPI_COMM_WORLD);
 }
};

template<class Vector, typename T>
void MPIgmpBcast(Vector& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
{
gmpBcast2<Vector, T> gB;
gB(A, ni, nj, src, Cptr);
}
template<class Vector, typename T>
void MPIgmpBcast(Vector& A, size_t ni, int src, LinBox::Communicator *Cptr)
{
gmpBcast<Vector, T>gB;
gB(A, ni, src, Cptr);  
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

/*
template <class Vector>
inline void MPIgmpBcast(Vector& B, size_t nj, int src, LinBox::Communicator *Cptr)
{
//  int B_mp_alloc[nj], B_a_size[nj]; 
int *B_mp_alloc=(int*)malloc(nj*sizeof(int));
int *B_a_size=(int*)malloc(nj*sizeof(int));
  unsigned lenB;  Givaro::Integer temp; 
  std::vector<mp_limb_t> B_mp_data;
  //MPI_Barrier(MPI_COMM_WORLD);
  if(src==Cptr->rank()){
    //std::cerr << "B=:= " << std::endl;
    //Split vector B into arrays
    __mpz_struct * ptr;
    for(size_t j=0;j<nj;j++){
      //std::cerr << B.getEntry(j)<< "\t" ; std::cerr<< std::endl;
      ptr = const_cast<__mpz_struct*>(B.getEntry(j).get_mpz());
      B_mp_alloc[j] = ptr->_mp_alloc;
      B_a_size[j] = ptr->_mp_size;
      mp_limb_t * a_array = ptr->_mp_d;
      for(long i=0; i< ptr->_mp_alloc; ++i){
	B_mp_data.push_back(a_array[i]);
      }
    }
    
    lenB = B_mp_data.size();
    
  }
  //Distribut Givaro::Integer through its elementary parts
 // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&B_mp_alloc[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&B_a_size[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if(src!=Cptr->rank()) B_mp_data.resize(lenB);
  MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  
  if(src!=Cptr->rank()){
    //Reconstruction of vector B
    //std::cerr << "received B::= " << std::endl;
    __mpz_struct * ptr2;
    size_t count=0;
    
    for(size_t j=0;j<nj;j++){ 
      ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
      ptr2->_mp_alloc = B_mp_alloc[j];
      ptr2->_mp_size = B_a_size[j];
      _mpz_realloc(ptr2,ptr2->_mp_alloc);
      for(long i=0; i< ptr2->_mp_alloc; ++i){
	ptr2->_mp_d[i] = (B_mp_data[i+count]);
      }count+=ptr2->_mp_alloc;
      B.setEntry(j,temp);
      //std::cerr << B.getEntry(j) << "\t" ; std::cerr<< std::endl; 
      
    }
    
  }//MPI_Barrier(MPI_COMM_WORLD);
  
}

*/


/*
//template<typename T>
//void MPIvecBcast(LinBox::SparseMatrix< T > & A, int ni, int nj, Communicator *Cptr)
template <class LinMat>
inline void MPIgmpBcast(LinMat& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
{ 
  //int A_mp_alloc[ni*nj], A_a_size[ni*nj];
int *A_mp_alloc=(int*)malloc(ni*nj*sizeof(int));
int *A_a_size=(int*)malloc(ni*nj*sizeof(int));
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

*/












template <class Vector>
inline void MPIgmpBcast2(Vector& B, size_t nj, int src, LinBox::Communicator *Cptr)
{ 
  int B_mp_alloc[nj], B_a_size[nj]; 
  unsigned lenB;  Givaro::Integer temp; 
  std::vector<mp_limb_t> B_mp_data;
  //MPI_Barrier(MPI_COMM_WORLD);
  if(src==Cptr->rank()){
    //std::cerr << "Split B := " << std::endl;
    //Split vector B into arrays
    __mpz_struct * ptr;
    for(size_t j=0;j<nj;j++){
      //std::cerr << B.getEntry(j)<< "\t" ; std::cerr<< std::endl;
      ptr = const_cast<__mpz_struct*>(B.getEntry(j).get_mpz());
      B_mp_alloc[j] = ptr->_mp_alloc;
      B_a_size[j] = ptr->_mp_size;
      mp_limb_t * a_array = ptr->_mp_d;
      for(long i=0; i< ptr->_mp_alloc; ++i){
	B_mp_data.push_back(a_array[i]);
      }
    }
    
    lenB = B_mp_data.size();
 mp_limb_t *_mp_data;
_mp_data = (mp_limb_t*) malloc(sizeof(mp_limb_t)*lenB);
for(long i=0;i<lenB;i++) _mp_data[i]= B_mp_data[i];
for(int i=1; i<Cptr->size(); i++){
MPI_Send(&B_mp_alloc,nj,MPI_INT,i,0,MPI_COMM_WORLD);
MPI_Send(&B_a_size,nj,MPI_INT,i,0,MPI_COMM_WORLD);
MPI_Send(&lenB,1,MPI_UNSIGNED,i,0,MPI_COMM_WORLD);
MPI_Send(_mp_data,lenB,chooseMPItype<mp_limb_t>::val,i,0,MPI_COMM_WORLD);
}
  }else{
MPI_Recv(&B_mp_alloc,nj,MPI_INT,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Recv(&B_a_size,nj,MPI_INT,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Recv(&lenB,1,MPI_UNSIGNED,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
mp_limb_t *_mp_data;
_mp_data = (mp_limb_t*) malloc(sizeof(mp_limb_t)*lenB);

MPI_Recv(_mp_data,lenB,chooseMPItype<mp_limb_t>::val,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
B_mp_data.resize(lenB);
for(long i=0;i<lenB;i++) B_mp_data[i]= _mp_data[i];


    //Reconstruction of vector B
    //std::cerr << "Reconstructed B:= " << std::endl;
    __mpz_struct * ptr2;
    size_t count=0;
    
    for(size_t j=0;j<nj;j++){ 
      ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
      ptr2->_mp_alloc = B_mp_alloc[j];
      ptr2->_mp_size = B_a_size[j];
      _mpz_realloc(ptr2,ptr2->_mp_alloc);
      for(long i=0; i< ptr2->_mp_alloc; ++i){
	ptr2->_mp_d[i] = (B_mp_data[i+count]);
      }count+=ptr2->_mp_alloc;
      B.setEntry(j,temp);
      //std::cerr << B.getEntry(j) << "\t" ; std::cerr<< std::endl; 
      
    }
    
  }//MPI_Barrier(MPI_COMM_WORLD);
  
}





template <class LinMat>
inline void MPIgmpBcast2(LinMat& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
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
MPI_Barrier(MPI_COMM_WORLD);

if(src==Cptr->rank()){

for(int i=1; i<Cptr->size(); i++){
MPI_Send(&A_mp_alloc,ni*nj,MPI_INT,i,0,MPI_COMM_WORLD);
MPI_Send(&A_a_size,ni*nj,MPI_INT,i,0,MPI_COMM_WORLD);
MPI_Send(&lenA,1,MPI_UNSIGNED,i,0,MPI_COMM_WORLD);
std::cerr << "lenA:="<<lenA << std::endl;
}
mp_limb_t *_mp_data;
_mp_data = (mp_limb_t*) malloc(lenA*sizeof(mp_limb_t));
for(long i=0;i<lenA;i++) _mp_data[i]= A_mp_data[i];
for(int i=1; i<Cptr->size(); i++) MPI_Send(_mp_data,lenA,chooseMPItype<mp_limb_t>::val,i,0,MPI_COMM_WORLD);

  }else{

MPI_Recv(&A_mp_alloc,ni*nj,MPI_INT,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Recv(&A_a_size,ni*nj,MPI_INT,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Recv(&lenA,1,MPI_UNSIGNED,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
mp_limb_t *_mp_data;
_mp_data = (mp_limb_t*) malloc(lenA*sizeof(mp_limb_t));

MPI_Recv(_mp_data,lenA,chooseMPItype<mp_limb_t>::val,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
A_mp_data.resize(lenA);
for(long i=0;i<lenA;i++) A_mp_data[i]= _mp_data[i];

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










template <class Vector>
inline void MPIgmpBcast3(Vector& B, size_t nj, int src, LinBox::Communicator *Cptr)
{ 
  //int B_mp_alloc[nj], B_a_size[nj]; 
int *B_mp_alloc=(int*)malloc(nj*sizeof(int));
int *B_a_size=(int*)malloc(nj*sizeof(int));
  unsigned lenB;  Givaro::Integer temp; 
  std::vector<mp_limb_t> B_mp_data;
  //MPI_Barrier(MPI_COMM_WORLD);
  if(src==Cptr->rank()){
    //std::cerr << ">>>> B:= " << std::endl;
    //Split vector B into arrays
    __mpz_struct * ptr;
    for(size_t j=0;j<nj;j++){
      //std::cerr << B.getEntry(j)<< "\t" ; std::cerr<< std::endl;
      ptr = const_cast<__mpz_struct*>(B.getEntry(j).get_mpz());
      B_mp_alloc[j] = ptr->_mp_alloc; //std::cerr <<"B_mp_alloc["<<j<<"] >>> "<<B_mp_alloc[j]<< std::endl;
      B_a_size[j] = ptr->_mp_size; //std::cerr <<"B_a_size>>>["<<j<<"] >>> "<<B_a_size[j]<< std::endl;
      mp_limb_t * a_array = ptr->_mp_d;
      for(long i=0; i< ptr->_mp_alloc; ++i){
	B_mp_data.push_back(a_array[i]);
      }
    }
    
    lenB = B_mp_data.size();
    
  }



int dl=500;
if(nj>dl){
  //Distribut Givaro::Integer through its elementary parts
for(int i=0;i<nj/dl;i++){
  MPI_Bcast(&B_mp_alloc[i*dl], dl, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&B_a_size[i*dl], dl, MPI_INT, 0, MPI_COMM_WORLD);
}
if(nj%dl>0){
  MPI_Bcast(&B_mp_alloc[nj/dl*dl], nj-nj/dl*dl, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&B_a_size[nj/dl*dl], nj-nj/dl*dl, MPI_INT, 0, MPI_COMM_WORLD);
}

}else{
  MPI_Bcast(&B_mp_alloc[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&B_a_size[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
}

MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
if(src!=Cptr->rank()) B_mp_data.resize(lenB);

if(lenB>dl){
for(int i=0;i<lenB/dl;i++){
  MPI_Bcast(&B_mp_data[i*dl], dl, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
}
if(lenB%dl>0){
  MPI_Bcast(&B_mp_data[lenB/dl*dl], lenB-lenB/dl*dl, chooseMPItype<mp_limb_t>::val, 0, MPI_COMM_WORLD);
}
}else{
  MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
}



  if(src!=Cptr->rank()){
    //Reconstruction of vector B
    //std::cerr << "<<<< B:= " << std::endl;
    __mpz_struct * ptr2;
    size_t count=0;
    
    for(size_t j=0;j<nj;j++){ 
      ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
      ptr2->_mp_alloc = B_mp_alloc[j]; //std::cerr <<"B_mp_alloc["<<j<<"] >>> "<<B_mp_alloc[j]<< std::endl;
      ptr2->_mp_size = B_a_size[j]; //std::cerr <<"B_a_size["<<j<<"] >>> "<<B_a_size[j]<< std::endl;

      _mpz_realloc(ptr2,ptr2->_mp_alloc);

      for(long i=0; i< ptr2->_mp_alloc; ++i){
	ptr2->_mp_d[i] = (B_mp_data[i+count]);
      }count+=ptr2->_mp_alloc;

      B.setEntry(j,temp);
      //std::cerr << B.getEntry(j) << "\t" ; std::cerr<< std::endl;       

    }
    
  }//MPI_Barrier(MPI_COMM_WORLD);
 free(B_mp_alloc);
 free(B_a_size);
}







//template<typename T>
//void MPIvecBcast(LinBox::SparseMatrix< T > & A, int ni, int nj, Communicator *Cptr)
template <class LinMat>
inline void MPIgmpBcast3(LinMat& A, size_t ni, size_t nj, int src, LinBox::Communicator *Cptr)
{ 

//  int A_mp_alloc[ni*nj], A_a_size[ni*nj];
int *A_mp_alloc=(int*)malloc(ni*nj*sizeof(int));
int *A_a_size=(int*)malloc(ni*nj*sizeof(int));

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

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int dl=500;
if(ni*nj>dl){
  //Distribut Givaro::Integer through its elementary parts
for(int i=0;i<ni*nj/dl;i++){
  MPI_Bcast(&A_mp_alloc[i*dl], dl, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A_a_size[i*dl], dl, MPI_INT, 0, MPI_COMM_WORLD);
}
if(ni*nj%dl>0){
  MPI_Bcast(&A_mp_alloc[ni*nj/dl*dl], ni*nj-ni*nj/dl*dl, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A_a_size[ni*nj/dl*dl], ni*nj-ni*nj/dl*dl, MPI_INT, 0, MPI_COMM_WORLD);
}

}else{
  MPI_Bcast(&A_mp_alloc[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A_a_size[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
}

MPI_Bcast(&lenA, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
if(src!=Cptr->rank()) A_mp_data.resize(lenA);

if(lenA>dl){
for(int i=0;i<lenA/dl;i++){
  MPI_Bcast(&A_mp_data[i*dl], dl, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
}
if(lenA%dl>0){
  MPI_Bcast(&A_mp_data[lenA/dl*dl], lenA-lenA/dl*dl, chooseMPItype<mp_limb_t>::val, 0, MPI_COMM_WORLD);
}
}else{
  MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 
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
 free(A_mp_alloc);
 free(A_a_size);
}











