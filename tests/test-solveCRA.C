
/* tests/test-solve.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by
 *
 * --------------------------------------------------------
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
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

#endif



#include <gmp++/gmp++.h>
#include <string>


int main (int argc, char **argv)
{
	bool pass = true;

	size_t ni = 10;
	size_t nj = ni;

	size_t r = 20;
	size_t N = 16;

#ifdef __LINBOX_HAVE_MPI
  Communicator *Cptr = NULL;
  Cptr = new Communicator(&argc, &argv);
  int max=10000;  
  Givaro::ZRing<Integer> ZZ;
  SparseMatrix<Givaro::ZRing<Integer> > A (ZZ,ni,nj);

typedef DenseVector<Givaro::ZRing<Integer> > DenseVector;
DenseVector X(ZZ, A.coldim()), X2(ZZ, A.coldim()),  B(ZZ, A.coldim());
Givaro::ZRing<Integer>::Element d;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
int B_mp_alloc[nj], B_a_size[nj]; 
int A_mp_alloc[ni*nj], A_a_size[ni*nj];
unsigned lenA,lenB;
std::vector<mp_limb_t> B_mp_data;
std::vector<mp_limb_t> A_mp_data;

  if(0==Cptr->rank()){
	static Argument args[] = {
		{ 'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT,     &ni },
		{ 'm', "-m M", "Set the mxaimum for the range of integers to generate.", TYPE_INT,     &max },
		END_OF_ARGUMENTS
	};	nj = ni;

std::cout<<"Computation is done over Q"<<std::endl;

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

//Split vector B into arrays
Givaro::Integer valeur2("448832189123456789123456789123456789123456789123456789123456789123456789"); 
   for (long j = 0; j < nj-1; ++j)
	B.setEntry(j,myrand(tmp, max));
//B.setEntry(nj-1,valeur2);
//std::cerr << "B:= " << std::endl; 
//__mpz_struct * ptr;
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


  }//End of BLock for process(0)


//Distribut Givaro::Integer through its elementary parts
MPI_Barrier(MPI_COMM_WORLD);
MPI_Bcast(&A_mp_alloc[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&A_a_size[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Bcast(&lenA, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
if(0!=Cptr->rank()) A_mp_data.resize(lenA);
MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD);
MPI_Bcast(&B_mp_alloc[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&B_a_size[0], nj, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
if(0!=Cptr->rank()) B_mp_data.resize(lenB);
MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);

 if(0!=Cptr->rank()){


//Reconstruction of matrix A
__mpz_struct * ptr2;
mp_limb_t * a_array;
size_t count=0;
 Givaro::Integer temp; 
//std::cerr << "received A:= " << std::endl;
   for (long i = 0; i < ni; ++i){
     for (long j = 0; j < nj; ++j){

ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());//ptr2 = const_cast<__mpz_struct*>(A.getEntry(i,j).get_mpz());
ptr2->_mp_alloc = A_mp_alloc[j+i*nj];
ptr2->_mp_size = A_a_size[j+i*nj];
   _mpz_realloc(ptr2,ptr2->_mp_alloc);
   for(size_t k=0; k< ptr2->_mp_alloc; ++k){
        ptr2->_mp_d[k] = (A_mp_data[k+count]);	
    }count+=ptr2->_mp_alloc;
A.setEntry(i,j,temp);//<----------------------------<<<<<
//std::cerr << A.getEntry(i,j) << "\t" ; std::cerr<< std::endl;
}
}


//Reconstruction of vector B
//std::cerr << "received B::= " << std::endl;
//__mpz_struct * ptr2;
count=0;
for(int j=0;j<nj;j++){ 
ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());//ptr2 = const_cast<__mpz_struct*>(B.getEntry(j).get_mpz());
ptr2->_mp_alloc = B_mp_alloc[j];
ptr2->_mp_size = B_a_size[j];
    _mpz_realloc(ptr2,ptr2->_mp_alloc);
    for(size_t i=0; i< ptr2->_mp_alloc; ++i){
        ptr2->_mp_d[i] = (B_mp_data[i+count]);
    }count+=ptr2->_mp_alloc;
B.setEntry(j,temp);//<----------------------------<<<<<
//std::cerr << B.getEntry(j) << "\t" ; std::cerr<< std::endl; 

}

  }


/***********************
* Results verification 
***********************/
RingCategories::IntegerTag tg;
Timer chrono;
   double starttime, endtime; 
   starttime = MPI_Wtime(); 
solveCRA (X2, d, A, B, tg, Method::Hybrid(*Cptr),Cptr);	
   endtime   = MPI_Wtime(); 
	

if(0==Cptr->rank()){


		// BlasElimination
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

			/*std::cout << "MPI CRA Solution is [";
			for(DenseVector::const_iterator it=X2.begin();it != X2.end(); ++it)
				ZZ.write(cout, *it) << " ";
			std::cout << "] / ";
			ZZ.write(std::cout, d) << std::endl;*/
			std::cout << "CPU time (seconds): " << endtime-starttime << std::endl;


	for (size_t j = 0 ; j < nj ; ++j){
		if(!Givaro::ZRing<Integer>().areEqual(X[j],X2[j])){
			std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			std::cerr << "                  ***Failed***                   " << std::endl;
			std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
			pass = false;
			break;
		}
	}


}

  MPI_Finalize();
return pass ? 0 : -1;
	
#else
Timer chrono;
                std::cout << "BlasElimination" << std::endl;
                chrono.start();
                solve (X, d, A, B, Method::BlasElimination());
                chrono.stop();
                std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;

		
		std::cout << "CRA" << std::endl;
                chrono.start();
		solveCRA (X2, d, A, B, tg, Method::Hybrid(*Cptr),Cptr);	
                chrono.stop();
		std::cout << "CPU time (seconds): " <<  chrono.usertime() << std::endl;
	for (size_t j = 0 ; j < nj ; ++j){
		if(!Givaro::ZRing<Integer>().areEqual(X[j],X2[j])){
			std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			std::cerr << "                  ***Failed***                   " << std::endl;
			std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
			pass = false;
			break;
		}
	}

	return pass ? 0 : -1;
#endif
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


	
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
