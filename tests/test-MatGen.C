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
/*
template <typename T > class chooseMPItype;
template <> struct chooseMPItype<unsigned int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED;};
template <> struct chooseMPItype<unsigned long long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG_LONG;};
template <> struct chooseMPItype<unsigned long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG;};


void printMPItype(MPI_Datatype type){ 
  if(MPI_UNSIGNED==type)  std::cout<<"MPI_UNSIGNED is to be used"<< std::endl;
  if(MPI_UNSIGNED_LONG_LONG==type)   std::cout<<"MPI_UNSIGNED_LONG_LONG is to be used"<< std::endl;
  if(MPI_UNSIGNED_LONG==type) std::cout<<"MPI_UNSIGNED_LONG is to be used"<<std::endl;
};
*/



template <class Vector>
//void MPIgmpBcast(Vector& B, int ni, int nj, Communicator *Cptr)
void MPIgmpBcast(Vector& B, int nj, Communicator *Cptr)
{ 
  int B_mp_alloc[nj], B_a_size[nj]; 
  unsigned lenB;  Givaro::Integer temp; 
  std::vector<mp_limb_t> B_mp_data;
  //MPI_Barrier(MPI_COMM_WORLD);
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
 // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&B_mp_alloc[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&B_a_size[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if(0!=Cptr->rank()) B_mp_data.resize(lenB);
  MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  
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
    
  }//MPI_Barrier(MPI_COMM_WORLD);
  
}

//template<typename T>
//void MPIvecBcast(LinBox::SparseMatrix< T > & A, int ni, int nj, Communicator *Cptr)
template <class LinMat>
void MPIgmpBcast(LinMat& A, int ni, int nj, Communicator *Cptr)
{ 
  int A_mp_alloc[ni*nj], A_a_size[ni*nj];
  unsigned lenA;
  std::vector<mp_limb_t> A_mp_data;
  //MPI_Barrier(MPI_COMM_WORLD);
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
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&A_mp_alloc[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A_a_size[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&lenA, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if(0!=Cptr->rank()) A_mp_data.resize(lenA);
  MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  
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
    
  }//MPI_Barrier(MPI_COMM_WORLD);
  
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
  long dg,ni,nj,max;  
  long tmp;

   dg=20, ni=3,nj=3,max=100; 

/*
//std::cout<<"argv size:"<<argc<<std::endl;
std::map<char,int> input = {{'n',1},{'m',2}};
std::map<char,int>::iterator it;
if(argc>1){

	if(argc-1==2){
	//std::cout<<"argv[1][1]:"<<argv[1][1]<<std::endl;
	it = input.find(argv[1][1]);
	//std::cout<<"idx:"<<it->second<<std::endl;
switch(it->second) {
    case 1 : ni=atoi(argv[2]); break;
    case 2 : max=atoi(argv[2]); break;
}
	}

	if(argc-1==4){
	//std::cout<<"argv[1][1]:"<<argv[1][1]<<std::endl;
	it = input.find(argv[1][1]);
	//std::cout<<"idx:"<<it->second<<std::endl;
switch(it->second) {
    case 1 : ni=atoi(argv[2]); break;
    case 2 : max=atoi(argv[2]); break;
}
	//std::cout<<"argv[3][1]:"<<argv[3][1]<<std::endl;
	it = input.find(argv[3][1]);
	//std::cout<<"idx:"<<it->second<<std::endl;
switch(it->second) {
    case 1 : ni=atoi(argv[4]); break;
    case 2 : max=atoi(argv[4]); break;
}
	}
}
*/

//std::cerr << "proc("<< Cptr->rank() <<") defaut: ni = nj = "<< nj << std::endl;
 	static Argument args[] = {
		{ 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
		{ 'm', "-m M", "Set the mxaimum for the range of integers to generate.", TYPE_INT,     &max },
		{ 'd', "-d M", "Set the mxaimum number of digits of integers to generate.", TYPE_INT,     &dg },
		END_OF_ARGUMENTS
	};	
parseArguments (argc, argv, args); nj = ni;
//std::cerr << "proc("<< Cptr->rank() <<") defaut: ni = nj = "<< nj << std::endl;


 MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); nj=ni;


//std::cerr << "proc("<< Cptr->rank() <<") init>> ni = nj = "<< nj << std::endl;

  Givaro::ZRing<Integer> ZZ;
  SparseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);

  typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;
  DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
  Givaro::ZRing<Integer>::Element d;



  if(0==Cptr->rank()){
/*
Givaro::Integer valeur("1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111");
 printf("%ld %ld M\n", ni, nj);
  for (long i = 0; i < ni; ++i)
    for (long j = 0; j < nj; ++j){
      if(i==j)A.setEntry(i,j,valeur); else A.setEntry(i,j,0); 
    }
  //A.setEntry(nj-1,nj-1,valeur);
long tmp;
  for (long j = 0; j < nj; ++j) B.setEntry(j,valeur*(myrand(tmp, 10)));




long unsigned int r=0;
std::cout<<"The rank of generated matrix A is:"<<r<<std::endl;
while(r!=ni){
    for (long i = 0; i < ni; ++i)
      for (long j = 0; j < nj; ++j){
	//A.setEntry(i,j,myrand(tmp, max));

result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
A.setEntry(i,j,res);
result.clear();

      }
    //A.setEntry(nj-1,nj-1,valeur);

    //Givaro::Integer valeur2("448832189123456789123456789123456789123456789123456789123456789123456789"); 
    for (long j = 0; j < nj; ++j){
      //B.setEntry(j,myrand(tmp, max));
    //B.setEntry(nj-1,valeur);

result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
B.setEntry(j,res);
result.clear();

    }

LinBox::rank (r, A);
}
*/





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
 std::cout << "The rank of generated matrix A is: "<<r<<std::endl; 


/*
      std::cerr << "Compute with B: " << std::endl;
      for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
      
      std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
      if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A::=",Tag::FileFormat::Maple) << ';' << std::endl;
*/

  }//End of BLock for process(0)

  //MPI data distribution for Integer type value
  MPIgmpBcast(A, ni, nj, Cptr);
  MPIgmpBcast(B, nj, Cptr);
MPI_Barrier(MPI_COMM_WORLD);  
std::cout<<"Computation is done over Q"<<std::endl;
  /***********************
   Results verification 
  ***********************/
  RingCategories::IntegerTag tg;
  Timer chrono;
  double starttime, endtime; 
  starttime = MPI_Wtime(); 
  solveCRA (X2, d, A, B, tg, Method::Hybrid(*Cptr),Cptr);	
  endtime   = MPI_Wtime(); 

 MPI_Barrier(MPI_COMM_WORLD);
DenseVector B2(ZZ, A.coldim());
  if(0==Cptr->rank()){
    
/*
      std::cerr << "Compute with B: " << std::endl;
      for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
      
      std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
      if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A::=",Tag::FileFormat::Maple) << ';' << std::endl;
    
*/    

    
/* 
    // BlasElimination
    //DenseVector X(ZZ, A.coldim());
    std::cout << "BlasElimination" << std::endl;
    chrono.start();
    solve (X, d, A, B, Method::BlasElimination());
    chrono.stop();
*/    
/*
    std::cout << "BlasElimination Solution is [";
    for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
      ZZ.write(cout, *it) << " ";
    std::cout << "] / ";
    ZZ.write(std::cout, d)<< std::endl;
    std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;
*/


    // solveCRA
    std::cout << "MPI solveCRA" << std::endl;
/*  
    std::cout << "MPI CRA Solution is [";
    for(DenseVector::const_iterator it=X2.begin();it != X2.end(); ++it)
      ZZ.write(cout, *it) << " ";
    std::cout << "] / ";
    ZZ.write(std::cout, d) << std::endl;
*/
    std::cout << "CPU time (seconds): " << endtime-starttime << std::endl;
    
/*    
    for (long j = 0 ; j < nj ; ++j){
      if(!Givaro::ZRing<Integer>().areEqual(X[j],X2[j])){
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "                  ***Failed***                   " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	break;
      }
    }
*/


DenseVector B2(ZZ, A.coldim());
DenseVector B3(ZZ, A.coldim());
/*A.apply(B2,X);
for (long j = 0 ; j < nj ; ++j) B3.setEntry(j,d*B.getEntry(j));

  for (long j = 0 ; j < nj ; ++j){
    if(!ZZ.areEqual(B2[j],B3[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "        The solution of solve using BlasElimination is incorrect         " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      break;
    }
  }*/

A.apply(B2,X2);
for (long j = 0 ; j < nj ; ++j) B3.setEntry(j,d*B.getEntry(j));

  for (long j = 0 ; j < nj ; ++j){
    if(!ZZ.areEqual(B3[j],B2[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "       The solution of solveCRA using BlasElimination is incorrect       " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      break;
    }
  }


   
  MPI_Finalize();
  return 0;    
  }
  MPI_Finalize();
  
  
  
  
#else
   std::cout<<"Computation is done over Q"<<std::endl;
  srand48(time(NULL));
  long dg=20, ni=3,nj=3,max=100;
  //	int offset = 0;
 	static Argument args[] = {
		{ 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
		{ 'm', "-m M", "Set the mxaimum for the range of integers to generate.", TYPE_INT,     &max },
		{ 'd', "-d M", "Set the mxaimum number of digits of integers to generate.", TYPE_INT,     &dg },
		END_OF_ARGUMENTS
	};	
	parseArguments (argc, argv, args); nj = ni;
 
  Givaro::ZRing<Integer> ZZ;
  DenseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);
  
  typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;
  DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
  Givaro::ZRing<Integer>::Element d;
/*
  long tmp;
Givaro::Integer valeur("11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111");
printf("%ld %ld M\n", ni, nj);
  for (long i = 0; i < ni; ++i)
    for (long j = 0; j < nj; ++j){
      if(i==j)A.setEntry(i,j,valeur); else A.setEntry(i,j,0); 
    }
  //A.setEntry(nj-1,nj-1,valeur);

  for (long j = 0; j < nj; ++j) B.setEntry(j,valeur*(j+1));
 */  
/*
const char * c;
std::string result;
    for (long i = 0; i < ni; ++i)
      for (long j = 0; j < nj; ++j){
result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
A.setEntry(i,j,res);
result.clear();
      }

    for (long j = 0; j < nj; ++j){

result = gmp_rand(dg);
c = result.c_str();
Givaro::Integer res(c);
B.setEntry(j,res);
result.clear();
    }
*/


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
  Timer chrono;
/*
  // Sequential BlasElimination
  //DenseVector X(ZZ, A.coldim());
  std::cout << "BlasElimination" << std::endl;
  chrono.start();
  solve (X, d, A, B, Method::BlasElimination());
  chrono.stop();
  
  std::cout << "BlasElimination Solution is  [";
  for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
    ZZ.write(cout, *it) << " ";
  std::cout << "] / ";
  ZZ.write(std::cout, d)<< std::endl;
  std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;
  */
  // Sequential solveCRA
  std::cout << "Sequential solveCRA" << std::endl;
  RingCategories::IntegerTag tg;
  chrono.start();
  solveCRA (X2, d, A, B, tg, Method::BlasElimination());	
  chrono.stop();
/*
  std::cout << "Sequential CRA Solution is  [";
  for(DenseVector::const_iterator it=X2.begin();it != X2.end(); ++it)
    ZZ.write(cout, *it) << " ";
  std::cout << "] / ";
  ZZ.write(std::cout, d) << std::endl;
*/
  std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DenseVector B2(ZZ, A.coldim());
DenseVector B3(ZZ, A.coldim());
A.apply(B2,X);
/*
for (long j = 0 ; j < nj ; ++j) B3.setEntry(j,d*B.getEntry(j));

  for (long j = 0 ; j < nj ; ++j){
    if(!ZZ.areEqual(B2[j],B3[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "        The solution of solve using BlasElimination is incorrect         " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      break;
    }
  }
*/
A.apply(B2,X2);
for (long j = 0 ; j < nj ; ++j) B3.setEntry(j,d*B.getEntry(j));

  for (long j = 0 ; j < nj ; ++j){
    if(!ZZ.areEqual(B3[j],B2[j])){   
      std::cerr << "\033[1;31m>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\033[0m" << std::endl;
      std::cerr << "\033[1;31m       The solution of solveCRA using BlasElimination is incorrect       \033[0m" << std::endl;
      std::cerr << "\033[1;31m<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\033[0m" << std::endl;
      break;
    }
  }

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  return 0;
  
#endif
  
}


