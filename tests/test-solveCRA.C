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

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template <class Field, class Matrix>
static bool checkInput (const Field  &F,
			 Matrix &A, Matrix &A2,
			 BlasVector<Field> &B,			 
			 BlasVector<Field> &B2
			 )
{
for (size_t i = 0 ; i < A.rowdim() ; ++i)
  for (size_t j = 0 ; j < A.coldim() ; ++j){
    if(!F.areEqual(A2.getEntry(i,j),A.getEntry(i,j))){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The input matrix is inconsistent                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      return false;
    }
  }
  for (size_t j = 0 ; j < B.size() ; ++j){
    if(!F.areEqual(B2.getEntry(j),B.getEntry(j))){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The input vector is inconsistent                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      return false;
    }
  }
  return true;
}


template <class Field>
static bool checkOutput (const Field  &F,
			 BlasVector<Field> &X,
			 BlasVector<Field> &X2
			 )
{


  for (size_t j = 0 ; j < X.size() ; ++j){
    if(!F.areEqual(X2.getEntry(j),X.getEntry(j))){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The output vector is inconsistent                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      return false;
    }
  }

  return true;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





template <class Field, class Matrix>
static bool checkResult (const Field  &F,
			 Matrix &A,
			 BlasVector<Field> &B,
			 BlasVector<Field> &X,
			 typename Field::Element &d)
{

  BlasVector<Field> B2(F, A.coldim());
  BlasVector<Field> B3(F, A.coldim());
  A.apply(B2,X);

Integer tmp;
  for (size_t j = 0 ; j < B.size() ; ++j){
  //tmp=B.getEntry(j);F.mulin(tmp,d);B3.setEntry(j,tmp);
    B3.setEntry(j,d*B.getEntry(j));
  }
  for (size_t j = 0 ; j < A.coldim() ; ++j){
    if(!F.areEqual(B2[j],B3[j])){
      std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
      std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      /*
      std::cerr << " B2["<<j<<"] := "<< B2[j] << std::endl;
      std::cerr << " B3["<<j<<"] := "<< B3[j] << std::endl;
      std::cerr << " d*B["<<j<<"] := "<< d*B.getEntry(j) << std::endl;
      */
      return false;
    }
    }
  return true;
}



template <class Field, class Matrix>
void genData (const Field  &F, Matrix  &Mat, size_t bits, int seed){
  typedef typename Field::RandIter RandIter;  
  
  RandIter RI(F, bits, seed) ;
  LinBox::RandomDenseMatrix<RandIter,Field>  RDM(F,RI);
  RDM.randomFullRank(Mat);


}


template <class Field>
void genData (const Field  &F, BlasVector<Field>  &Vec, size_t bits, int seed){
  typedef typename Field::RandIter RandIter; 
      
  RandIter RI(F, bits, seed);
  Vec.random(RI);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN ROI
/////////////////////////////////////////////////////////////////////////////////////////////////////
template <class Field>
bool test_set_with_field(BlasVector<Field> &X2,
	      BlasMatrix<Field> &A,
	      BlasVector<Field> &B
	      ,int bits
#ifdef __LINBOX_HAVE_MPI
	      , Communicator *Cptr   //, BlasVector<Field> &X2_cp
#endif

	      ){
  bool tag = false;
  Givaro::ZRing<Integer> F;
  Givaro::ZRing<Integer>::Element d;
  std::cerr<<"Computation is done over Q"<<std::endl;
  
#ifdef __LINBOX_HAVE_MPI
  std::cerr << "MPI solveCRA" << std::endl;
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
if(0 == Cptr->rank())
#endif
//  checkOutput (F, X2, X2_cp) ;//<<---------------------------------
  
  
#ifdef __LINBOX_HAVE_MPI
  endtime   = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
#else
  end = omp_get_wtime();
#endif
  //  DenseVector B2(F, A.coldim());

  
#ifdef __LINBOX_HAVE_MPI
  if(0 == Cptr->rank()){  

    std::cout << "Total CPU time (seconds): " << endtime-starttime << std::endl;
#else
    std::cout << "Total CPU time (seconds): " << end-start << std::endl;
#endif

    tag=checkResult (F, A, B, X2, d);/*
    if(!tag){
        B.write(std::cout << " >>>> Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
	    A.write(std::cout << " >>>> Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;
	    std::cout << " >>>> Found d:="<<d<< std::endl;
	    std::cout << " >>>> Found X:=";
	    for(long i=0; i<X2.size();i++)std::cout << X2[i] << std::endl;
    }*/
#ifdef __LINBOX_HAVE_MPI
  }
#endif 
  
#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&tag, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif

  return tag;	      
	      }
/////////////////////////////////////////////////////////////////////////////////////////////////////
// END ROI
/////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
bool test_set(BlasVector<Givaro::ZRing<Integer> > &X2,
	      BlasMatrix<Givaro::ZRing<Integer> > &A,
	      BlasVector<Givaro::ZRing<Integer> > &B
,	      int bits
#ifdef __LINBOX_HAVE_MPI
	      , Communicator *Cptr
#endif
	      ){
  bool tag = false;
  Givaro::ZRing<Integer> F(bits);
  Givaro::ZRing<Integer>::Element d;
  std::cerr<<"Computation is done over Q"<<std::endl;
  
#ifdef __LINBOX_HAVE_MPI
  std::cerr << "MPI solveCRA" << std::endl;
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
  //  DenseVector B2(F, A.coldim());

  
#ifdef __LINBOX_HAVE_MPI
  if(0 == Cptr->rank()){  

    std::cout << "Total CPU time (seconds): " << endtime-starttime << std::endl;
#else
    std::cout << "Total CPU time (seconds): " << end-start << std::endl;
#endif

    tag=checkResult (F, A, B, X2, d);/*
    if(!tag){
        B.write(std::cout << " >>>> Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
	    A.write(std::cout << " >>>> Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;
	    std::cout << " >>>> Found d:="<<d<< std::endl;
	    std::cout << " >>>> Found X:=";
	    for(long i=0; i<X2.size();i++)std::cout << X2[i] << std::endl;
    }*/
#ifdef __LINBOX_HAVE_MPI
  }
#endif 
  
#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&tag, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif
  return tag;
}

#endif

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//subroutine for defaut seed if no user provided seed as parameter
	uint64_t getSeed(){
		struct timeval tp;
		gettimeofday(&tp, 0) ;
        return static_cast<uint64_t> (tp.tv_usec + tp.tv_sec*1000000);
	}


	void rand_param_vary(size_t &n, size_t &ni, size_t &bits, size_t &bitsize, bool &peak){

n = ni;

    n = rand() % ni + 1;
    if (n < ni / 2 && n % 2 == 0 && !peak) n = 1;

    bits = rand() % bitsize + 1;
    if (bits < bitsize / 2 && bitsize % 2 == 0 && !peak) bits = 1;
    
    peak = !peak;
    
    std::cout << " Test with dimension: " << n << " x " << n << std::endl;
    std::cout << " Test with bitsize: " << bits << std::endl;



	}

template <class Field>
void prepare_data_with_field(size_t bits, int seed,
          BlasVector<Field> &X2,
	      BlasMatrix<Field> &A,
	      BlasVector<Field> &B
	      
#ifdef __LINBOX_HAVE_MPI
	      , Communicator *Cptr
#endif
	      ){
Field F;

std::cout << " Proc("<<Cptr->rank()<<") <<<<<<<<<<<<<<<<< seed:= "<<seed<<std::endl;
std::cout << " Proc("<<Cptr->rank()<<") <<<<<<<<<<<<<<<<< n:= "<<A.coldim()<<std::endl;
std::cout << " Proc("<<Cptr->rank()<<") <<<<<<<<<<<<<<<<< bits:= "<<bits<<std::endl;

#ifdef __LINBOX_HAVE_MPI  
    if(0==Cptr->rank()){
#endif
    //std::cerr << " seed := "<<seed<<std::endl;

      genData (F, A, bits, seed);
      genData (F, B, bits, seed);

#ifdef __LINBOX_HAVE_MPI 	
    }//End of BLock for process(0)
#endif   



//	B.write(std::cout << " Proc("<<Cptr->rank()<<") >>>> Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
//	A.write(std::cout << " Proc("<<Cptr->rank()<<") >>>> Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;



//std::cout << " Proc("<<Cptr->rank()<<") >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
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
//std::cout << " Proc("<<Cptr->rank()<<") <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << std::endl;

    
		     //Check if data are correctly distributed to all processes
//	B.write(std::cout << " Proc("<<Cptr->rank()<<") <<<< Compute with B:=",Tag::FileFormat::Maple) << ';' << std::endl;
//	A.write(std::cout << " Proc("<<Cptr->rank()<<") <<<< Compute with A:=",Tag::FileFormat::Maple) << ';' << std::endl;


	      }


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



int main(int argc, char ** argv)
{
  
#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
#endif
  size_t bits=10;
  size_t bitsize=10;
  size_t niter=1;
  size_t ni=1;

  size_t n=1;
  int q=-1;
  bool peak = false;
  bool loop=false;
  
//@ Use the following subroutine for defaut seed if no user provided seed as parameter 
    int seed = 0; 
//@ print out used seed before each test


 Argument args[] = {
    { 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
    { 'b', "-b B", "Set the bitsize for input value.", TYPE_INT,     &bitsize },
    { 'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT,     &niter },
    { 'q', "-q Q", "Set the randomness of test (<0 for random and >0 for derterministic).",    TYPE_INT , &q },
    { 'l', "-loop Y/N", "Set if run the test in an infinite loop.", TYPE_BOOL , &loop },
    { 's', "-s S", "Set the seed to fill the input matrices.", TYPE_INT,     &seed },
    END_OF_ARGUMENTS
  };	
  parseArguments (argc, argv, args); 


  

///////////////////////////////////////////////////////////////////////////////////////////////////////
#if 0

#ifdef __LINBOX_HAVE_MPI

  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

#else
  srand (time(NULL));
#endif

    
  Givaro::ZRing<Integer> F(bitsize);
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;

  DenseMatrix<Givaro::ZRing<Integer> > A (F,ni,ni);
  DenseVector X(F, A.rowdim()), X2(F, A.coldim()),  B(F, A.coldim());
  
  for(long j=0;loop ||j<(long)niter;j++){  
    
#ifdef __LINBOX_HAVE_MPI  
    if(0==Cptr->rank()){
#endif


    std::cout << " Test with dimension: " << ni << " x " << ni << std::endl;
    std::cout << " Test with bitsize: " << bitsize << std::endl;


      genData (F, A, bitsize, seed);//genData (A, bits);
      genData (F, B, bitsize, seed);//genData (B, bits);






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



    if(!test_set(X2, A, B, bitsize
#ifdef __LINBOX_HAVE_MPI
		 , Cptr
#endif
		 )){

		 break;
		 }

  }
#else /////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&loop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
//  MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&q, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif


    
  Givaro::ZRing<Integer> F;
  typedef BlasVector<Givaro::ZRing<Integer> > DenseVector;

//Givaro::ZRing<Integer> F(bits);

    
    

#ifdef __LINBOX_HAVE_MPI
    if(0==Cptr->rank()){
#endif

        if(q<0){   
            seed = getSeed();//seed = 1247341449;
            rand_param_vary(n, ni, bits, bitsize, peak);
        }else{
            bits = bitsize;  
            n = ni;
        }

std::cout << " >>>>>>>>>>>>>>>>>> seed:= "<<seed<<std::endl;
std::cout << " >>>>>>>>>>>>>>>>>> n:= "<<n<<std::endl;
std::cout << " >>>>>>>>>>>>>>>>>> bits:= "<<bits<<std::endl;
    
#ifdef __LINBOX_HAVE_MPI 	
    }
#endif    


#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bits, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif



    DenseMatrix<Givaro::ZRing<Integer> > A (F,n,n);
    DenseVector X(F, A.coldim()), X2(F, A.coldim()),  B(F, A.coldim());


        
 prepare_data_with_field(bits, seed, X2, A, B	      //<-------------------instant modification
#ifdef __LINBOX_HAVE_MPI
	      , Cptr
#endif
	      );


//////////////////////////////Make a copy to compare at each iteration/////////////////////////////////////////
    DenseMatrix<Givaro::ZRing<Integer> > A_cp(F,n,n);
    DenseVector X2_cp(F, A.coldim()),  B_cp(F, A.coldim());	



    
 prepare_data_with_field(bits, seed, X2_cp, A_cp, B_cp	      //<-------------------instant modification
#ifdef __LINBOX_HAVE_MPI
	      , Cptr
#endif
	      );
	      
checkInput (F, A, A_cp, B, B_cp) ;


test_set_with_field<Givaro::ZRing<Integer>>(X2_cp, A, B, bits
#ifdef __LINBOX_HAVE_MPI
		 , Cptr          //, X2_cp
#endif

		 );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////      
	      
long computtIMES= 0;

  for(size_t j=0;loop || j<niter;j++){  

//-------------------------Coould be wrapped into a tempalted test set subroutine------------------------->

//------------------Coould be wrapped into a tempalted subroutine--------------------->


#ifdef __LINBOX_HAVE_MPI
    if(0==Cptr->rank()){
#endif

        if(q<0){   
            seed = getSeed();
            rand_param_vary(n, ni, bits, bitsize, peak);
        }else{
            bits = bitsize;  
            n = ni;
        }

std::cout << " >>>>>>>>>>>>>>>>>> seed:= "<<seed<<std::endl;
std::cout << " >>>>>>>>>>>>>>>>>> n:= "<<n<<std::endl;
std::cout << " >>>>>>>>>>>>>>>>>> bits:= "<<bits<<std::endl;
    
#ifdef __LINBOX_HAVE_MPI 	
    }//End of BLock for process(0)
#endif    




#ifdef __LINBOX_HAVE_MPI
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bits, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    DenseMatrix<Givaro::ZRing<Integer> > A (F,n,n);
    DenseVector X(F, A.coldim()), X2(F, A.coldim()),  B(F, A.coldim());


        
 prepare_data_with_field(bits, seed, X2, A, B	      //<-------------------instant modification
#ifdef __LINBOX_HAVE_MPI
	      , Cptr
#endif
	      );
	      size_t r;
std::cout << " A rank:= "<<LinBox::rank(r,A)<<std::endl;
//<--------------------Coould be wrapped into a tempalted subroutine----------------------	
	
//@FutherImprovement: Use the templated test function with field as a template at compile time: 	    ok = ok && run_with_field<Modular<int32_t> >(q,b,n,iters,seed)     

std::cout << " ################################ Until now computed "<<computtIMES<<" times! "<<std::endl;

if(!checkInput (F, A, A_cp, B, B_cp)) std::cerr<<" Proc("<<Cptr->rank()<<") >>>>>>>>>>>>>>>>>>>>> Input is inconsistent !!!!!!!!!!!!!!!!! "<<std::endl;;
      
    if(!test_set_with_field<Givaro::ZRing<Integer>>(X2, A, B, bits
#ifdef __LINBOX_HAVE_MPI
		 , Cptr        // , X2_cp
#endif

		 )){
		 break;
		 }
computtIMES++;
MPI_Barrier(MPI_COMM_WORLD);
//<-------------------------Coould be wrapped into a tempalted test set subroutine-------------------------

  }
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __LINBOX_HAVE_MPI
  delete Cptr;  
#endif
    return 0;
    
}

