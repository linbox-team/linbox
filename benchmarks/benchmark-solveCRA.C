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
#include "linbox/util/mpi-gmp.inl"
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



#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char ** argv)
{
  
  
#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
  size_t dg,ni,nj,max;  
  
  dg=10, ni=3,nj=3,max=10000; 
  
  static Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
    { 'm', "-m M", "Set the mxaimum for the range of integers to generate.", TYPE_INT,     &max },
    { 'd', "-d M", "Set the mxaimum number of digits of integers to generate.", TYPE_INT,     &dg },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); nj = ni;
  
  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); nj=ni;
  
  Givaro::ZRing<Integer> ZZ;
  SparseMatrix<Givaro::ZRing<Integer>> A (ZZ,ni,nj);
  
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;
  DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
  Givaro::ZRing<Integer>::Element d;
  
  
  
  if(0==Cptr->rank()){
    
    
    //ofstream myfile;  ofstream myfile2;
    //myfile.open ("matrix.txt"); myfile2.open ("vec.txt");
    //myfile <<ni<<" "<<nj<<" M\n"; myfile2 <<nj<< " 1"<<" M\n";
    long unsigned int r=0;
    
    if(dg>20){
      const char * c;
      std::string result;
      
      while(r!=ni){

	for (size_t i = 0; i < ni; ++i)
	  for (size_t j = 0; j < nj; ++j){
	    result = gmp_rand(dg);
	    c = result.c_str();
	    Givaro::Integer res(c);
	    A.setEntry(i,j,res);
	    result.clear();
	    //myfile <<i+1<<" "<<j+1<<" "<<A.getEntry(i,j)<<"\n";
	  }
	
	
	for (size_t j = 0; j < nj; ++j){
	  result = gmp_rand(dg);
	  c = result.c_str();
	  Givaro::Integer res(c);
	  B.setEntry(j,res);
	  result.clear();
	  //myfile2 <<j+1<< " 1 "<<B.getEntry(j)<<"\n";
	}
	
	LinBox::rank (r, A);
      }
      
    }else{
      
      long tmp; 
      while(r!=ni){
	
	for (size_t i = 0; i < ni; ++i)
	  for (size_t j = 0; j < nj; ++j){
	    A.setEntry(i,j,myrand(tmp,max));
	    //myfile <<i+1<<" "<<j+1<<" "<<A.getEntry(i,j)<<"\n";
	  }
	
	
	for (size_t j = 0; j < nj; ++j){
	  B.setEntry(j,myrand(tmp,max));
	  //myfile2 <<j+1<< " 1 "<<B.getEntry(j)<<"\n";
	}
	
	LinBox::rank (r, A); // Check if the generated matrix is invertible
      }
    }
    //myfile <<"0  0  0";
    //myfile2 <<"0  0  0";
    // myfile.close();
    // myfile2.close();
    LinBox::rank (r, A); std::cout<<"The rank of generated matrix A is:"<<r<<std::endl;  
    
    
    /*
      std::cerr << ">>>>Compute with B: " << std::endl;      
      for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
      
      std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
      if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
    */
    
  }//End of BLock for process(0)

  //distribute big integer compatible data
  {
    double starttime, endtime; 
    MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();
    //MPI data distribution for Integer type value
    MPIgmpBcast(A, ni, nj, 0, Cptr);
    MPIgmpBcast(B, ni, 0, Cptr);
    MPI_Barrier(MPI_COMM_WORLD);
    endtime   = MPI_Wtime(); 
    //std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
  }
  
  
  
  //Check if data are correctly distributed to all processes
  /*
  if(0!=Cptr->rank()){std::cerr <<"process("<<Cptr->rank()<< ")<<<<Compute with B: " << std::endl;
  for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; }
  
  std::cerr <<"process("<<Cptr->rank()<< ")<<<<Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
  if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
  */
  
  
  
  
  

  MPI_Barrier(MPI_COMM_WORLD);  
  std::cout<<"Computation is done over Q"<<std::endl;
  
  /***********************
   Results verification 
  ***********************/
  RingCategories::IntegerTag tg;
  Timer chrono;
  double starttime, endtime; 
  starttime = MPI_Wtime(); 
  solveCRA (X2, d, A, B, tg, /*Method::BlasElimination() */Method::Hybrid(*Cptr),Cptr);	
  endtime   = MPI_Wtime(); 
  
  MPI_Barrier(MPI_COMM_WORLD);
  DenseVector B2(ZZ, A.coldim());
  if(0 == Cptr->rank()){  
    
    /*
      std::cerr << "Compute with B: " << std::endl;
      for(int j=0;j<nj;j++) std::cerr << B.getEntry(j) << std::endl; 
      
      std::cout << "Compute with A: " << A.rowdim() << " by " << A.coldim() << std::endl;
      if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cout << "A::=",Tag::FileFormat::Maple) << ';' << std::endl;
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
    
    
    
    
    DenseVector B2(ZZ, A.coldim());
    DenseVector B3(ZZ, A.coldim());
    
    A.apply(B2,X2);
    for (size_t j = 0 ; j < nj ; ++j) B3.setEntry(j,d*B.getEntry(j));
    
    for (size_t j = 0 ; j < nj ; ++j){
      if(!ZZ.areEqual(B3[j],B2[j])){
	std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
	std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
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
  
  const char * c;
  std::string result;
  long unsigned int r=0;
  /*
    ofstream myfile;  ofstream myfile2;
    myfile.open ("matrix.txt"); myfile2.open ("vec.txt");
    myfile <<ni<<" "<<nj<<" M\n"; myfile2 <<nj<< " 1"<<" M\n";
  */
  
  if(dg>20){
    const char * c;
    std::string result;
    
    while(r!=ni){
      
      for (size_t i = 0; i < ni; ++i)
	for (size_t j = 0; j < nj; ++j){
	  result = gmp_rand(dg);
	  c = result.c_str();
	  Givaro::Integer res(c);
	  A.setEntry(i,j,res);
	  result.clear();
	  //myfile <<i+1<<" "<<j+1<<" "<<A.getEntry(i,j)<<"\n";
	}

      
      for (size_t j = 0; j < nj; ++j){
	result = gmp_rand(dg);
	c = result.c_str();
	Givaro::Integer res(c);
	B.setEntry(j,res);
	result.clear();
	//myfile2 <<j+1<< " 1 "<<B.getEntry(j)<<"\n";
    }
      
      LinBox::rank (r, A);
    }
    
  }else{
    
    long tmp; 
    while(r!=ni){
      
      for (size_t i = 0; i < ni; ++i)
	for (size_t j = 0; j < nj; ++j){
	  A.setEntry(i,j,myrand(tmp,max));
	  //myfile <<i+1<<" "<<j+1<<" "<<A.getEntry(i,j)<<"\n";
	}
      
      
      for (size_t j = 0; j < nj; ++j){
	B.setEntry(j,myrand(tmp,max));
	//myfile2 <<j+1<< " 1 "<<B.getEntry(j)<<"\n";
      }
      
      LinBox::rank (r, A); // Check if the generated matrix is invertible
    }
  }
  //myfile <<"0  0  0";
  //myfile2 <<"0  0  0";
  // myfile.close();
  // myfile2.close();
  
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
  
  
  DenseVector B2(ZZ, A.coldim());
  DenseVector B3(ZZ, A.coldim());
  A.apply(B2,X);
  
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
  

  
  return 0;
  
#endif
  
}


