/* Copyright (C) 2018 The LinBox group
 * Updated by Hongguang Zhu <zhuhongguang2014@gmail.com>
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


template <typename T > class chooseMPItype;
template <> struct chooseMPItype<int>{ static constexpr MPI_Datatype val = MPI_INT;};
template <> struct chooseMPItype<unsigned int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED;};
template <> struct chooseMPItype<unsigned long long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG_LONG;};
template <> struct chooseMPItype<unsigned long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG;};
template <> struct chooseMPItype<float>{ static constexpr MPI_Datatype val = MPI_FLOAT;};
template <> struct chooseMPItype<double>{ static constexpr MPI_Datatype val = MPI_DOUBLE;};



template <class Matrix>
void gmp_unpackMat(Matrix& M, int *A_mp_alloc, int *A_a_size, std::vector<mp_limb_t>& A_mp_data){
size_t ni=M.rowdim(), nj=M.rowdim();

        __mpz_struct * ptr;
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){

                ptr = const_cast<__mpz_struct*>(M.getEntry(i,j).get_mpz());
                A_mp_alloc[j+i*nj] = ptr->_mp_alloc;
                A_a_size[j+i*nj] = ptr->_mp_size;
                mp_limb_t * a_array = ptr->_mp_d;
                for(long k=0; k< ptr->_mp_alloc; ++k)
                    A_mp_data.push_back(a_array[k]);
                
            }
        }
}


template <class Matrix>
void gmp_packMat(Matrix& M, int *A_mp_alloc, int *A_a_size, std::vector<mp_limb_t>& A_mp_data){
    size_t ni=M.rowdim(), nj=M.rowdim();
    __mpz_struct * ptr2;
    size_t count=0;
    Givaro::Integer temp; 
    
    for (size_t i = 0; i < ni; ++i){
        for (size_t j = 0; j < nj; ++j){
            
            ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
            ptr2->_mp_alloc = A_mp_alloc[j+i*nj];
            ptr2->_mp_size = A_a_size[j+i*nj];
            _mpz_realloc(ptr2,ptr2->_mp_alloc);
            for(long k=0; k< ptr2->_mp_alloc; ++k){
                ptr2->_mp_d[k] = (A_mp_data[k+count]);	
            }count+=ptr2->_mp_alloc;
            M.setEntry(i,j,temp);
            
        }
    }
}


template <class Matrix>
void gmp_unpackSparseMat(Matrix& M, std::vector<int>& A_mp_alloc, std::vector<int>&  A_a_size, std::vector<mp_limb_t>& A_mp_data, std::vector<long>&  A_index){
    size_t ni=M.rowdim(), nj=M.rowdim();
    Givaro::ZRing<Givaro::Integer> ZZ;
    
    __mpz_struct * ptr;
    for (size_t i = 0; i < ni; ++i){
        for (size_t j = 0; j < nj; ++j){
            if(!ZZ.areEqual(M.getEntry(i,j),ZZ.zero)){
                
                ptr = const_cast<__mpz_struct*>(M.getEntry(i,j).get_mpz());
                A_mp_alloc.push_back( ptr->_mp_alloc );
                A_a_size.push_back( ptr->_mp_size );
                mp_limb_t * a_array = ptr->_mp_d;
                A_index.push_back(i);A_index.push_back(j);
                for(long k=0; k< ptr->_mp_alloc; ++k)
                    A_mp_data.push_back(a_array[k]);
            }
            
        }
    }
}


template <class Matrix>
void gmp_packSparseMat(Matrix& M, std::vector<int>& A_mp_alloc, std::vector<int>&  A_a_size, std::vector<mp_limb_t>& A_mp_data, std::vector<long>&  A_index){
    size_t ni=M.rowdim(), nj=M.rowdim();
    Givaro::ZRing<Givaro::Integer> ZZ;
    __mpz_struct * ptr2;
    size_t count=0;
    Givaro::Integer temp; 
    
    for (size_t i = 0; i < ni; ++i){
        for (size_t j = 0; j < nj; ++j){
            
            M.setEntry(i,j,ZZ.zero);
            
        }
    }
    
    for (size_t i = 0; i < A_a_size.size(); ++i){
        
        //ptr = const_cast<__mpz_struct*>(M.getEntry().get_mpz());
        ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
        ptr2->_mp_alloc = A_mp_alloc[i];
        ptr2->_mp_size = A_a_size[i];
        _mpz_realloc(ptr2,ptr2->_mp_alloc);
        
        for(long k=0; k< ptr2->_mp_alloc; ++k){
            ptr2->_mp_d[k] = (A_mp_data[k+count]);	
        }count+=ptr2->_mp_alloc;
        
        M.setEntry(A_index[i*2],A_index[i*2+1],temp);
        
    }
    
}


template <class Vector>
void gmp_unpackVec(Vector& V, int *B_mp_alloc, int *B_a_size, std::vector<mp_limb_t>& B_mp_data){
    size_t nj = V.size();
    Givaro::Integer temp; 
    
    //Split vector B into arrays
    __mpz_struct * ptr;
    for(size_t j=0;j<nj;j++){

        ptr = const_cast<__mpz_struct*>(V.getEntry(j).get_mpz());
        B_mp_alloc[j] = ptr->_mp_alloc;
        B_a_size[j] = ptr->_mp_size;
        mp_limb_t * a_array = ptr->_mp_d;
                for(long i=0; i< ptr->_mp_alloc; ++i){
                    B_mp_data.push_back(a_array[i]);
                }
    }
}

template <class Vector>
void gmp_packVec(Vector& V, int *B_mp_alloc, int *B_a_size, std::vector<mp_limb_t>& B_mp_data){
    size_t nj = V.size();
    //Reconstruction of vector B
    __mpz_struct * ptr2;
    size_t count=0;
    Givaro::Integer temp;
    for(size_t j=0;j<nj;j++){ 
        ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
        ptr2->_mp_alloc = B_mp_alloc[j];
        ptr2->_mp_size = B_a_size[j];
        _mpz_realloc(ptr2,ptr2->_mp_alloc);
        for(long i=0; i< ptr2->_mp_alloc; ++i){
            ptr2->_mp_d[i] = (B_mp_data[i+count]);
        }count+=ptr2->_mp_alloc;
        V.setEntry(j,temp);
        
    }
}


#ifndef __LINBOX_mpicpp_INL
#define __LINBOX_mpicpp_INL
//  BRYAN - changed the iterator specific lines to simply
//  sizeof (int *)'s...  also updated the .h file to not
//  define these functions as having default parameters,
//  as mpicxx compiler was having trouble with this

namespace LinBox
{
    
	Communicator::Communicator(MPI_Comm comm = MPI_COMM_NULL) :
		_mpi_comm(comm), _mpi_boss(false)
	{}
    
	// MPI_initializing constructor
	// When this communicator is destroyed MPI is shut down (finalized).
	Communicator::Communicator(int* ac, char*** av) :
		_mpi_comm(MPI_COMM_WORLD), _mpi_boss(true)
	{	MPI_Init(ac, av); }
    
	// copy constructor
	Communicator::Communicator(const Communicator& D) :
		_mpi_comm(D._mpi_comm), _mpi_boss(false), stat(D.stat)
	{}
    
	Communicator::~Communicator()
	{	if (_mpi_boss) MPI_Finalize(); }
    
	// accessors
	int Communicator::size()
	{	int s; MPI_Comm_size(_mpi_comm, &s); return s; }
    
	int Communicator::rank()
	{	int r; MPI_Comm_rank(_mpi_comm, &r); return r; }
    
	MPI_Status Communicator::status()
	{	return stat; }
    
	MPI_Comm Communicator::mpi_communicator()
	{	return _mpi_comm; }
    
	// peer to peer communication
    
	template < class Ptr >
	void Communicator::send( Ptr b, Ptr e, int dest, int tag)
	{
		MPI_Send( &*b,
                  (e - b)*sizeof(typename Ptr::value_type),
                  MPI_BYTE,
                  dest,
                  tag,
                  _mpi_comm);
	}
    
	template < class Ptr >
	void Communicator::ssend( Ptr b, Ptr e, int dest, int tag)
	{	MPI_Ssend( &b[0],
                   //(e - b)*sizeof(iterator_traits<Ptr>::value_type),
                   (e-b)*sizeof(int *),
                   MPI_BYTE,
                   dest,
                   tag,
                   _mpi_comm);
	}
    
    
	template < class Ptr >
	void Communicator::recv( Ptr b, Ptr e, int dest, int tag)
	{
		MPI_Recv( &b[0],
                  (e - b)*sizeof(typename Ptr::value_type),
                  MPI_BYTE,
                  dest,
                  tag,
                  _mpi_comm,
                  &stat);
	}
    /*
      template < class X >
      void Communicator::send( X *b, X *e, int dest, int tag)
      {
      MPI_Send( b,
      (e - b)*sizeof(X),
      MPI_BYTE,
      dest,
      tag,
      _mpi_comm);
      }
    */
	template < class X >
	void Communicator::recv( X *b, X *e, int dest, int tag)
	{
		MPI_Recv( b,
                  (e - b)*sizeof(X),
                  MPI_BYTE,
                  dest,
                  tag,
                  _mpi_comm,
                  &stat);
	}
    
    
	// whole object send and recv
    //	void Communicator::send( Ptr b, Ptr e, int dest, int tag)
    
    
	template < class X >
	void Communicator::send( X& b, int dest/*, int tag = 0*/)
	{	MPI_Send(&b,
                 sizeof(X),
                 MPI_BYTE,
                 dest,
                 0,
                 _mpi_comm);
	}
    
    template <class Field>
    void Communicator::send (BlasMatrix<Field>& M, int dest){
        
        size_t ni=M.rowdim(), nj=M.rowdim(); 
        std::vector<typename Field::Element>A_mp_data; A_mp_data.resize(ni*nj); 
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){
                
                A_mp_data[j+i*nj] = M.getEntry(i,j);
                
            }
        }
        
        MPI_Send(&A_mp_data[0], ni*nj, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD);
        
    }
    template <class Field>
    void Communicator::send (SparseMatrix<Field>& M, int dest){
        Field ZZ;
        size_t ni=M.rowdim(), nj=M.rowdim(); 
        
        std::vector<typename Field::Element> datum;
        std::vector<long> index;
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){
                if(!ZZ.areEqual(M.getEntry(i,j),ZZ.zero)){
                    datum.push_back( M.getEntry(i,j)); index.push_back(i);index.push_back(j); 
                };
                
            }
        }
        
        int len = datum.size();
        MPI_Send(&len, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        MPI_Send(&datum[0], len, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD);
        
        MPI_Send(&index[0], 2*len, MPI_LONG , dest, 0, MPI_COMM_WORLD);
        
    }
    
    
    template <class Field>
    void Communicator::send (BlasVector<Field>& V, int dest){
        
        size_t nj=V.size();
        
        std::vector<typename Field::Element>B_mp_data; B_mp_data.resize(nj);

        for(size_t j=0;j<nj;j++){
            
            B_mp_data[j] = V.getEntry(j);
            
        }
        
        MPI_Send(&B_mp_data[0], nj, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD);              
        
    }
    
    template <class Matrix>
    void Communicator::send_integerMat (Matrix& A, int dest){
        size_t ni=A.rowdim(), nj=A.rowdim();
        
        std::vector<int>A_mp_alloc; A_mp_alloc.resize(ni*nj);
        std::vector<int>A_a_size; A_a_size.resize(ni*nj);
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        
        gmp_unpackMat(A, &A_mp_alloc[0], &A_a_size[0], A_mp_data);
        lenA = A_mp_data.size();
        
        MPI_Send(&A_mp_alloc[0], ni*nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        MPI_Send(&A_a_size[0], ni*nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Send(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD);
        
    }
    
    template <class Matrix>
    void Communicator::send_integerSparseMat (Matrix& A, int dest){
        
        std::vector<int>A_mp_alloc;
        std::vector<int>A_a_size; 
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        std::vector<long> index;
        gmp_unpackSparseMat(A, A_mp_alloc, A_a_size, A_mp_data, index);
        
        lenA = A_mp_alloc.size();        
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Send(&A_mp_alloc[0], lenA, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        lenA = A_a_size.size();
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Send(&A_a_size[0], lenA, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        lenA = index.size();
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Send(&index[0], lenA, MPI_LONG,  dest, 0, MPI_COMM_WORLD);
        lenA = A_mp_data.size();
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Send(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD);
        
    }
    
    template <class Vector>
    void Communicator::send_integerVec (Vector& B, int dest){
        size_t nj=B.size();
        std::vector<int>B_mp_alloc; B_mp_alloc.resize(nj);
        std::vector<int>B_a_size; B_a_size.resize(nj);
        
        unsigned lenB;  Givaro::Integer temp; 
        std::vector<mp_limb_t> B_mp_data;
        
        gmp_unpackVec(B, &B_mp_alloc[0], &B_a_size[0], B_mp_data); 
        lenB = B_mp_data.size();
        
        MPI_Send(&B_mp_alloc[0], nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        
        MPI_Send(&B_a_size[0], nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        MPI_Send(&lenB, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Send(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD);              
        
    }
    
    template <>
    void Communicator::send (BlasMatrix<Givaro::ZRing<Integer> >& M, int dest){
        send_integerMat(M, dest);
    }
    template <>
    void Communicator::send (SparseMatrix<Givaro::ZRing<Integer> >& M, int dest){
        send_integerSparseMat(M, dest);
    }
    template <>
    void Communicator::send (DenseVector<Givaro::ZRing<Integer> >& V, int dest){
        
        send_integerVec(V, dest);
        
    }
    
    

#if 1
	template < class X >
	void Communicator::isend( X& b, int dest/*, int tag = 0*/)
	{	MPI_Isend(&b,
                  sizeof(X),
                  MPI_BYTE,
                  dest,
                  0,
                  _mpi_comm,&req);
	}
    
    template <class Field>
    void Communicator::isend (BlasMatrix<Field>& M, int dest){
        
        size_t ni=M.rowdim(), nj=M.rowdim(); 
        std::vector<typename Field::Element>A_mp_data; A_mp_data.resize(ni*nj); 
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){

                A_mp_data[j+i*nj] = M.getEntry(i,j);
                
            }
        }
        
        MPI_Send(&A_mp_data[0], ni*nj, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD);

    }
    template <class Field>
    void Communicator::isend (SparseMatrix<Field>& M, int dest){
        
        size_t ni=M.rowdim(), nj=M.rowdim(); 
        
        std::vector<typename Field::Element> datum;
        std::vector<long> index;
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){
                if(0!=M.getEntry(i,j)){datum.push_back( M.getEntry(i,j)); index.push_back(i);index.push_back(j);};
                
            }
        }

        int len = datum.size();
        MPI_Send(&len, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        MPI_Isend(&datum[0], len, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD,&req);

        MPI_Isend(&index[0], 2*len, MPI_LONG , dest, 0, MPI_COMM_WORLD,&req);

    }
    
    
    template <class Field>
    void Communicator::isend (BlasVector<Field>& V, int dest){
        
        size_t nj=V.size();
        
        std::vector<typename Field::Element>B_mp_data; B_mp_data.resize(nj); 

        for(size_t j=0;j<nj;j++){

            B_mp_data[j] = V.getEntry(j);
            
        }
        
        MPI_Isend(&B_mp_data[0], nj, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD,&req);              
  
        
    }

#if 1    
    template <class Matrix>
    void Communicator::isend_integerMat (Matrix& A, int dest){

        size_t ni=A.rowdim(), nj=A.rowdim();
        
        std::vector<int>A_mp_alloc; A_mp_alloc.resize(ni*nj); //int *A_mp_alloc=(int*)malloc(ni*nj*sizeof(int));
        std::vector<int>A_a_size; A_a_size.resize(ni*nj); //int *A_a_size=(int*)malloc(ni*nj*sizeof(int));  
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;

        gmp_unpackMat(A, &A_mp_alloc[0], &A_a_size[0], A_mp_data);
        lenA = A_mp_data.size();        


        
        MPI_Isend(&A_mp_alloc[0], ni*nj, MPI_INT,  dest, 0, MPI_COMM_WORLD,&req);
        MPI_Isend(&A_a_size[0], ni*nj, MPI_INT,  dest, 0, MPI_COMM_WORLD,&req);

        MPI_Isend(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD,&req);
        MPI_Isend(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD,&req);
        MPI_Wait(&req,  MPI_STATUS_IGNORE);
  
        
    }
    
    template <class Matrix>
    void Communicator::isend_integerSparseMat (Matrix& A, int dest){
        //        size_t ni=A.rowdim(), nj=A.rowdim();
        
        std::vector<int>A_mp_alloc;
        std::vector<int>A_a_size;
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        std::vector<long> index;
        gmp_unpackSparseMat(A, A_mp_alloc, A_a_size, A_mp_data, index);
        
        lenA = A_mp_alloc.size();        
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Isend(&A_mp_alloc[0], lenA, MPI_INT,  dest, 0, MPI_COMM_WORLD,&req);
        lenA = A_a_size.size(); 
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Isend(&A_a_size[0], lenA, MPI_INT,  dest, 0, MPI_COMM_WORLD,&req);
        lenA = index.size();
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Isend(&index[0], lenA, MPI_LONG,  dest, 0, MPI_COMM_WORLD,&req);
        lenA = A_mp_data.size();
        MPI_Send(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);

//        MPI_Send(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD);
        MPI_Isend(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD,&req);
        MPI_Wait(&req,  MPI_STATUS_IGNORE);
        
    }
#endif

    template <class Vector>
    void Communicator::isend_integerVec (Vector& B, int dest){
        size_t nj=B.size();
        std::vector<int>B_mp_alloc; B_mp_alloc.resize(nj); 
        std::vector<int>B_a_size; B_a_size.resize(nj);
        
        unsigned lenB;  Givaro::Integer temp; 
        std::vector<mp_limb_t> B_mp_data;
        
        gmp_unpackVec(B, &B_mp_alloc[0], &B_a_size[0], B_mp_data); 
        lenB = B_mp_data.size();
        
        MPI_Isend(&B_mp_alloc[0], nj, MPI_INT,  dest, 0, MPI_COMM_WORLD,&req);
        
        MPI_Isend(&B_a_size[0], nj, MPI_INT,  dest, 0, MPI_COMM_WORLD,&req);
        MPI_Send(&lenB, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Isend(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD,&req); 
        MPI_Wait(&req,  MPI_STATUS_IGNORE);            

        
    }
    
    template <>
    void Communicator::isend (BlasMatrix<Givaro::ZRing<Integer> >& M, int dest){
        isend_integerMat(M, dest);
    }
    template <>
    void Communicator::isend (SparseMatrix<Givaro::ZRing<Integer> >& M, int dest){
        isend_integerSparseMat(M, dest);
    }
    template <>
    void Communicator::isend (DenseVector<Givaro::ZRing<Integer> >& V, int dest){        
        isend_integerVec(V, dest);
    }

#endif


	template < class X >
	void Communicator::ssend( X& b, int dest/*, int tag = 0*/)
	{	MPI_Ssend(&b,
                 sizeof(X),
                 MPI_BYTE,
                 dest,
                 0,
                 _mpi_comm);
	}
    
    template <class Field>
    void Communicator::ssend (BlasMatrix<Field>& M, int dest){
        
        size_t ni=M.rowdim(), nj=M.rowdim(); 
        std::vector<typename Field::Element>A_mp_data; A_mp_data.resize(ni*nj); 
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){
                
                A_mp_data[j+i*nj] = M.getEntry(i,j);
                
            }
        }
        
        MPI_Ssend(&A_mp_data[0], ni*nj, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD);
        
    }
    template <class Field>
    void Communicator::ssend (SparseMatrix<Field>& M, int dest){
        Field ZZ;
        size_t ni=M.rowdim(), nj=M.rowdim(); 
        
        std::vector<typename Field::Element> datum;
        std::vector<long> index;
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){
                if(!ZZ.areEqual(M.getEntry(i,j),ZZ.zero)){
                    datum.push_back( M.getEntry(i,j)); index.push_back(i);index.push_back(j); 
                };
                
            }
        }
        
        int len = datum.size();
        MPI_Ssend(&len, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&datum[0], len, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD);
        
        MPI_Ssend(&index[0], 2*len, MPI_LONG , dest, 0, MPI_COMM_WORLD);
        
    }
    
    
    template <class Field>
    void Communicator::ssend (BlasVector<Field>& V, int dest){
        
        size_t nj=V.size();
        
        std::vector<typename Field::Element>B_mp_data; B_mp_data.resize(nj);

        for(size_t j=0;j<nj;j++){
            
            B_mp_data[j] = V.getEntry(j);
            
        }
        
        MPI_Ssend(&B_mp_data[0], nj, chooseMPItype<typename Field::Element>::val, dest, 0, MPI_COMM_WORLD);              
        
    }
    
    template <class Matrix>
    void Communicator::ssend_integerMat (Matrix& A, int dest){
        size_t ni=A.rowdim(), nj=A.rowdim();
        
        std::vector<int>A_mp_alloc; A_mp_alloc.resize(ni*nj);
        std::vector<int>A_a_size; A_a_size.resize(ni*nj);
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        
        gmp_unpackMat(A, &A_mp_alloc[0], &A_a_size[0], A_mp_data);
        lenA = A_mp_data.size();
        
        MPI_Ssend(&A_mp_alloc[0], ni*nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&A_a_size[0], ni*nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD);
        
    }
    
    template <class Matrix>
    void Communicator::ssend_integerSparseMat (Matrix& A, int dest){
        
        std::vector<int>A_mp_alloc;
        std::vector<int>A_a_size; 
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        std::vector<long> index;
        gmp_unpackSparseMat(A, A_mp_alloc, A_a_size, A_mp_data, index);
        
        lenA = A_mp_alloc.size();        
        MPI_Ssend(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&A_mp_alloc[0], lenA, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        lenA = A_a_size.size();
        MPI_Ssend(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&A_a_size[0], lenA, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        lenA = index.size();
        MPI_Ssend(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&index[0], lenA, MPI_LONG,  dest, 0, MPI_COMM_WORLD);
        lenA = A_mp_data.size();
        MPI_Ssend(&lenA, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD);
        
    }
    
    template <class Vector>
    void Communicator::ssend_integerVec (Vector& B, int dest){
        size_t nj=B.size();
        std::vector<int>B_mp_alloc; B_mp_alloc.resize(nj);
        std::vector<int>B_a_size; B_a_size.resize(nj);
        
        unsigned lenB;  Givaro::Integer temp; 
        std::vector<mp_limb_t> B_mp_data;
        
        gmp_unpackVec(B, &B_mp_alloc[0], &B_a_size[0], B_mp_data); 
        lenB = B_mp_data.size();
        
        MPI_Ssend(&B_mp_alloc[0], nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        
        MPI_Ssend(&B_a_size[0], nj, MPI_INT,  dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&lenB, 1, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD);
        MPI_Ssend(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val, dest, 0, MPI_COMM_WORLD);              
        
    }
    
    template <>
    void Communicator::ssend (BlasMatrix<Givaro::ZRing<Integer> >& M, int dest){
        ssend_integerMat(M, dest);
    }
    template <>
    void Communicator::ssend (SparseMatrix<Givaro::ZRing<Integer> >& M, int dest){
        ssend_integerSparseMat(M, dest);
    }
    template <>
    void Communicator::ssend (DenseVector<Givaro::ZRing<Integer> >& V, int dest){
        
        ssend_integerVec(V, dest);
        
    }



	template < class X >
	void Communicator::bsend(X& b, int dest /*, int tag = 0*/)
	{	MPI_Bsend( &b,
                   sizeof(X),
                   MPI_BYTE,
                   dest,
                   0,
                   _mpi_comm);
	}
    
    
	template < class X >
	void Communicator::recv( X& b, int src /*, int tag = 0*/)
	{	MPI_Recv( &b,
                  sizeof(X),
                  MPI_BYTE,
                  src,
                  0,
                  _mpi_comm,
                  &stat);
	}
    
    template <class Field>
    void Communicator::recv (BlasMatrix<Field>& M, int src){
        
        size_t ni=M.rowdim(), nj=M.rowdim();   
        std::vector<typename Field::Element>A_mp_data; A_mp_data.resize(ni*nj); 

        MPI_Recv(&A_mp_data[0], ni*nj, chooseMPItype<typename Field::Element>::val,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){
                
                M.setEntry(i,j,A_mp_data[j+i*nj]);

            }
        }
        
    }
    template <class Field>
    void Communicator::recv (SparseMatrix<Field>& M, int src){
        
        size_t ni=M.rowdim(), nj=M.rowdim();   
        
        std::vector<typename Field::Element> datum;
        std::vector<long> index;
        int len;
        MPI_Recv(&len, 1, MPI_INT,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        datum.resize(len);
        MPI_Recv(&datum[0], len, chooseMPItype<typename Field::Element>::val,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        index.resize(2*len);
        MPI_Recv(&index[0], 2*len, MPI_LONG,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (size_t i = 0; i < ni; ++i){
            for (size_t j = 0; j < nj; ++j){
                
                M.setEntry(i,j,0);
                
            }
        }
        
        size_t i =0;
        for (size_t j = 0; j < len; ++j) {
            
            M.setEntry(index[i],index[i+1],datum[j]);
            i+=2;
        }
        
    }
    
    template <class Field>
    void Communicator::recv (BlasVector<Field>& V, int src){
        
        size_t nj=V.size();   
        std::vector<typename Field::Element>B_mp_data; B_mp_data.resize(nj); 
        MPI_Recv(&B_mp_data[0], nj, chooseMPItype<typename Field::Element>::val,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (size_t j = 0; j < nj; ++j){
            V.setEntry(j,B_mp_data[j]);
        }
    }
    
    template <class Matrix>
    void Communicator::recv_integerMat (Matrix& A, int src){
        size_t ni=A.rowdim(), nj=A.rowdim();
        
        std::vector<int>A_mp_alloc; A_mp_alloc.resize(ni*nj); 
        std::vector<int>A_a_size; A_a_size.resize(ni*nj); 
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        std::vector<int> index;
        
        MPI_Recv(&A_mp_alloc[0], ni*nj, MPI_INT,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&A_a_size[0], ni*nj, MPI_INT,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&lenA, 1, MPI_UNSIGNED, 0,  src, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        A_mp_data.resize(lenA);
        MPI_Recv(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        gmp_packMat(A, &A_mp_alloc[0], &A_a_size[0], A_mp_data);
  
    }
    
    template <class Matrix>
    void Communicator::recv_integerSparseMat (Matrix& A, int src){
        
        std::vector<int>A_mp_alloc; 
        std::vector<int>A_a_size;  
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        std::vector<long> index;
        
        MPI_Recv(&lenA, 1, MPI_UNSIGNED, 0,  src, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
        A_mp_alloc.resize(lenA);
        MPI_Recv(&A_mp_alloc[0], lenA, MPI_INT,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 
        MPI_Recv(&lenA, 1, MPI_UNSIGNED, 0,  src, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       
        A_a_size.resize(lenA);
        MPI_Recv(&A_a_size[0], lenA, MPI_INT,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        MPI_Recv(&lenA, 1, MPI_UNSIGNED, 0,  src, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
        index.resize(lenA);
        MPI_Recv(&index[0], lenA, MPI_LONG,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        MPI_Recv(&lenA, 1, MPI_UNSIGNED, 0,  src, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
        A_mp_data.resize(lenA);
        MPI_Recv(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
        gmp_packSparseMat(A, A_mp_alloc, A_a_size, A_mp_data,index);
        
    }
    
    template <class Vector>
    void Communicator::recv_integerVec (Vector& B, int src){
        size_t nj=B.size();
        
        std::vector<int>B_mp_alloc; B_mp_alloc.resize(nj); 
        std::vector<int>B_a_size; B_a_size.resize(nj); 
        unsigned lenB;  Givaro::Integer temp; 
        std::vector<mp_limb_t> B_mp_data;
        
        MPI_Recv(&B_mp_alloc[0], nj, MPI_INT,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&B_a_size[0], nj, MPI_INT,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&lenB, 1, MPI_UNSIGNED, 0,  src, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        B_mp_data.resize(lenB);
        MPI_Recv(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        gmp_packVec(B, &B_mp_alloc[0], &B_a_size[0], B_mp_data);         
        
    }

    template<>
    void Communicator::recv (BlasMatrix<Givaro::ZRing<Integer> > & M, int src){
        recv_integerMat(M, src);
    }
    template<>
    void Communicator::recv (SparseMatrix<Givaro::ZRing<Integer> > & M, int src){
        recv_integerSparseMat(M, src);
    }
    template<>
    void Communicator::recv (DenseVector<Givaro::ZRing<Integer> > & V, int src){
        recv_integerVec(V, src);
    }
    
    // Broadcasts
    template <class X>
    void Communicator::bcast (X& x, int src){
        MPI_Bcast( &x, sizeof(X), MPI_BYTE, src, _mpi_comm);
    }

    template <class Field>
    void Communicator::bcast (BlasMatrix<Field>& M, int src){
        size_t ni = M.rowdim(), nj = M.rowdim();
        
        std::vector<typename Field::Element>A_mp_data; A_mp_data.resize(ni*nj); 
        
        unsigned lenA = ni*nj;
        
        if(src==rank()){
            
            for (size_t i = 0; i < ni; ++i){
                for (size_t j = 0; j < nj; ++j){
                    
                    A_mp_data[j+i*nj] = M.getEntry(i,j);                
                    
                }
            }
            
        }
        //Distribut Givaro::Integer through its elementary parts
        MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<typename Field::Element>::val,  0, MPI_COMM_WORLD);
        
        if(src!=rank()){
            //Reconstruction of matrix A part by part
            for (size_t i = 0; i < ni; ++i){
                for (size_t j = 0; j < nj; ++j){
                    
                    M.setEntry(i,j,A_mp_data[j+i*nj]);
                    
                }
            }
            
        }
        
    }
    
    template <class Field>
    void Communicator::bcast (SparseMatrix<Field>& M, int src){
        size_t ni = M.rowdim(), nj = M.rowdim();
        
        std::vector<typename Field::Element>A_mp_data;
        int lenA =0;
        
        if(src==rank()){

            for (size_t i = 0; i < ni; ++i){
                for (size_t j = 0; j < nj; ++j){
                    
                    if(M.getEntry(i,j)!=0){ 
                        A_mp_data.push_back( M.getEntry(i,j) );
                        A_mp_data.push_back( i );
                        A_mp_data.push_back( j );
                    }
                }
            }
            lenA = A_mp_data.size();
        }
        //Distribut Givaro::Integer through its elementary parts
        MPI_Bcast(&lenA, 1, MPI_INT,  0, MPI_COMM_WORLD);
        if(src!=rank()) A_mp_data.resize(lenA);
        
        MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<typename Field::Element>::val,  0, MPI_COMM_WORLD);
        
        if(src!=rank()){
            //Reconstruction of matrix A
            
            for (size_t i = 0; i < ni; ++i){
                for (size_t j = 0; j < nj; ++j){
                    
                    M.setEntry(i,j,0);
                    
                }
            }
            for (size_t i = 0; i < A_mp_data.size(); i+=3) M.setEntry(A_mp_data[i+1],A_mp_data[i+2],A_mp_data[i]);
            
        }
        
    }
    
    template <class Field>
    void Communicator::bcast (BlasVector<Field>& V, int src){
        
        size_t nj = V.size();      
        std::vector<typename Field::Element>B_mp_data; B_mp_data.resize(nj); 
        unsigned lenA = nj;
        
        if(src==rank()){
            
            //Split Matrix A into arrays 
            for (size_t j = 0; j < nj; ++j){
                B_mp_data[j] = V.getEntry(j);
            }
        }
        //Distribut Givaro::Integer through its elementary parts
        MPI_Bcast(&B_mp_data[0], lenA, chooseMPItype<typename Field::Element>::val,  0, MPI_COMM_WORLD);
        
        if(src!=rank()){
            //Reconstruction of matrix A
               for (size_t j = 0; j < nj; ++j){
                   
                   V.setEntry(j,B_mp_data[j]);
               }
         
        }
    }
    
    template <class X>
    void Communicator::bcast (X* b, X* e, int src){
        MPI_Bcast (b, (e-b)*sizeof(X), MPI_BYTE, src, _mpi_comm);
    }
    
    template <class Matrix>
    void Communicator::bcast_integerMat (Matrix& A, int src){
        size_t ni=A.rowdim(), nj=A.rowdim();
        std::vector<int>A_mp_alloc; A_mp_alloc.resize(ni*nj);
        std::vector<int>A_a_size; A_a_size.resize(ni*nj);
        
        unsigned lenA;
        std::vector<mp_limb_t> A_mp_data;
        if(src==rank()){

            gmp_unpackMat(A, &A_mp_alloc[0], &A_a_size[0], A_mp_data);
            lenA = A_mp_data.size();
        }
        //Distribut Givaro::Integer through its elementary parts
        MPI_Bcast(&A_mp_alloc[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&A_a_size[0], ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&lenA, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if(src!=rank()) A_mp_data.resize(lenA);
        MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
        
        if(src!=rank()){
            
            gmp_packMat(A, &A_mp_alloc[0], &A_a_size[0], A_mp_data);
            
        }
        
    }
    
    
    template <class Matrix>
    void Communicator::bcast_integerSparseMat(Matrix& A, int src){
        
        std::vector<int>A_mp_alloc;
        std::vector<int>A_a_size;       
        std::vector<long>index;
        unsigned lenA,lenB,lenC,lenD;
        std::vector<mp_limb_t> A_mp_data;
        
        if(src==rank()){
            
            gmp_unpackSparseMat(A, A_mp_alloc, A_a_size, A_mp_data,index);
            lenA = A_mp_alloc.size();
            lenB = A_a_size.size();
            lenC = index.size();
            lenD = A_mp_data.size();
        }
        
        //Distribut Givaro::Integer through its elementary parts
        MPI_Bcast(&lenA, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if(src!=rank()) A_mp_alloc.resize(lenA);
        MPI_Bcast(&A_mp_alloc[0], lenA, MPI_INT, 0, MPI_COMM_WORLD);
        
        MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if(src!=rank()) A_a_size.resize(lenB);
        MPI_Bcast(&A_a_size[0], lenB, MPI_INT, 0, MPI_COMM_WORLD);
        
        MPI_Bcast(&lenC, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if(src!=rank()) index.resize(lenC);
        MPI_Bcast(&index[0], lenC, MPI_LONG, 0, MPI_COMM_WORLD);
        
        MPI_Bcast(&lenD, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if(src!=rank()) A_mp_data.resize(lenD);
        MPI_Bcast(&A_mp_data[0], lenD, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);

        if(src!=rank()){

            gmp_packSparseMat(A, A_mp_alloc, A_a_size, A_mp_data,index);
 
        }
        
    }    
    
    template <class Vector>
    void Communicator::bcast_integerVec (Vector& B, int src){
        size_t nj=B.size();
        
        std::vector<int>B_mp_alloc; B_mp_alloc.resize(nj); 
        std::vector<int>B_a_size; B_a_size.resize(nj); 
        unsigned lenB;  Givaro::Integer temp; 
        std::vector<mp_limb_t> B_mp_data;
        
        if(src==rank()){
            
            gmp_unpackVec(B, &B_mp_alloc[0], &B_a_size[0], B_mp_data); 
            lenB = B_mp_data.size();            
        }
        
        //Distribut Givaro::Integer through its elementary parts
        MPI_Bcast(&B_mp_alloc[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&B_a_size[0], nj, MPI_INT, 0, MPI_COMM_WORLD);
        
        MPI_Bcast(&lenB, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if(src!=rank()) B_mp_data.resize(lenB);
        MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
        
        if(src!=rank()){
            
            gmp_packVec(B, &B_mp_alloc[0], &B_a_size[0], B_mp_data);            
        }
        
    }
    
    template <>
    void Communicator::bcast (BlasMatrix<Givaro::ZRing<Integer> > & M, int src){
        bcast_integerMat(M, src);
    }
    template <>
    void Communicator::bcast (SparseMatrix<Givaro::ZRing<Integer> > & M, int src){
        bcast_integerSparseMat(M, src);
    }
    template<>
    void Communicator::bcast (DenseVector<Givaro::ZRing<Integer> > & V, int src){
        bcast_integerVec(V, src);
    }
    
   
    
	template < class X >
	void Communicator::buffer_attach(X b)
	{
		MPI_Buffer_attach( malloc(sizeof(X) *  60) ,
                           sizeof(X) * 60 );
	}
    
	template < class X >
	int Communicator::buffer_detach(X &b, int *size)
	{
        return MPI_Buffer_detach( &b,
                                  size);
	}
    
	// collective communication
	template < class Ptr, class Function_object >
	void Communicator::reduce( Ptr bloc, Ptr eloc, Ptr bres, Function_object binop, int root)
	{}
    
	// member access
	MPI_Status Communicator::get_stat()
	{
		return stat;
	}
    
} // namespace LinBox
#endif // __LINBOX_mpicpp_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
