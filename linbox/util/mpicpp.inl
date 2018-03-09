/* Copyright (C) 2010 LinBox
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
template <> struct chooseMPItype<unsigned int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED;};
template <> struct chooseMPItype<unsigned long long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG_LONG;};
template <> struct chooseMPItype<unsigned long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG;};

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
	template < class X >
	void Communicator::send( X& b, int dest /*, int tag = 0 */)
	{	MPI_Send(&b,
			 sizeof(X),
			 MPI_BYTE,
			 dest,
			 0,
			 _mpi_comm);
	}

	template < class X >
	void Communicator::ssend( X& b, int dest /*, int tag = 0 */)
	{	MPI_Ssend(&b,
			  sizeof(X),
			  MPI_BYTE,
			  dest,
			  0,
			  _mpi_comm);
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
	void Communicator::recv( X& b, int dest /*, int tag = 0*/)
	{	MPI_Recv( &b,
			  sizeof(X),
			  MPI_BYTE,
			  dest,
			  0,
			  _mpi_comm,
			  &stat);
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Broadcasts
        template <class X>
        void Communicator::bcast (X& x, int src){
                MPI_Bcast( &x, sizeof(X), MPI_BYTE, src, _mpi_comm);
        }
        template <class Field>
        void Communicator::bcast (DenseMatrix<Field>& M, int src){
                MPI_Bcast(M._ptr, M._ptr+M.rowdim()*M.coldim(), src, _mpi_comm);
        }
        template <class Field>
        void Communicator::bcast (SparseMatrix<Field>& M, int src){
                MPI_Bcast(M._ptr, M._ptr+M.rowdim()*M.coldim(), src, _mpi_comm);
        }
        template <class Field>
        void Communicator::bcast (DenseVector<Field>& V, int src){
                // MPI_Bcast(V._ptr, V._ptr+V.size(), src, _mpi_comm);
        }
        
        template <class X>
        void Communicator::bcast (X* b, X* e, int src){
                MPI_Bcast (b, (e-b)*sizeof(X), MPI_BYTE, src, _mpi_comm);
        }

        template <class Matrix>
        void Communicator::bcast_integer2 (Matrix& A, int src){
                size_t ni=A.rowdim(), nj=A.rowdim();

                        int *A_mp_alloc=(int*)malloc(ni*nj*sizeof(int));
                        int *A_a_size=(int*)malloc(ni*nj*sizeof(int));
                        unsigned lenA;
                        std::vector<mp_limb_t> A_mp_data;
                        //MPI_Barrier(MPI_COMM_WORLD);
                        if(src==rank()){
                                
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
                        if(src!=rank()) A_mp_data.resize(lenA);
                        MPI_Bcast(&A_mp_data[0], lenA, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
                        //MPI_Barrier(MPI_COMM_WORLD);
                        
                        if(src!=rank()){
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
        void Communicator::bcast_integer (Vector& B, int src){
                size_t ni=B.rowdim(), nj=B.rowdim();
                nj=nj>ni?ni:nj;
                int *B_mp_alloc=(int*)malloc(nj*sizeof(int));
                int *B_a_size=(int*)malloc(nj*sizeof(int));
                unsigned lenB;  Givaro::Integer temp; 
                std::vector<mp_limb_t> B_mp_data;
                //MPI_Barrier(MPI_COMM_WORLD);
                if(src==rank()){
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
                if(src!=rank()) B_mp_data.resize(lenB);
                MPI_Bcast(&B_mp_data[0], lenB, chooseMPItype<mp_limb_t>::val,  0, MPI_COMM_WORLD);
                //MPI_Barrier(MPI_COMM_WORLD);
                
                if(src!=rank()){
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
        
        
        template<>
        void Communicator::bcast (DenseMatrix<Givaro::ZRing<Integer> > & M, int src){
                bcast_integer2(M, src);
        }
        template<>
        void Communicator::bcast (DenseVector<Givaro::ZRing<Integer> > & V, int src){
                bcast_integer(V, src);
        }
        
        template <>
        void Communicator::bcast (SparseMatrix<Givaro::ZRing<Integer> > & M, int src){
                bcast_integer2(M, src);
        }
        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////        

        
	template < class X >
	void Communicator::buffer_attach(X b)
	{
		MPI_Buffer_attach( malloc(sizeof(X) *  60) ,
				   sizeof(X) * 60 );
	}

	template < class X >
	int Communicator::buffer_detach(X &b, int *size)
	{  return MPI_Buffer_detach( &b,
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


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

