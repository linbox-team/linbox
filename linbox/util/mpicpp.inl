#ifndef __MPICPP_INL_
#define __MPICPP_INL_
//  BRYAN - changed the iterator specific lines to simply
//  sizeof (int *)'s...  also updated the .h file to not
//  define these functions as having default parameters,
//  as mpicxx compiler was having trouble with this

namespace LinBox {

Communicator::Communicator(MPI_Comm comm = MPI_COMM_NULL)
	: _mpi_comm(comm), _mpi_boss(false)
{}

// MPI_initializing constructor
// When this communicator is destroyed MPI is shut down (finalized).
Communicator::Communicator(int* ac, char*** av)
	: _mpi_comm(MPI_COMM_WORLD), _mpi_boss(true)
{	MPI_Init(ac, av); }

// copy constructor
Communicator::Communicator(const Communicator& D) 
	: _mpi_comm(D._mpi_comm), _mpi_boss(false), stat(D.stat)
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
void Communicator::send( Ptr b, Ptr e, int dest, int tag = 0)
{
	MPI_Send( &*b, 
	(e - b)*sizeof(typename Ptr::value_type), 
	MPI_BYTE,
	dest, 
	tag,
	_mpi_comm);
}

template < class Ptr >
void Communicator::ssend( Ptr b, Ptr e, int dest, int tag = 0)
{	MPI_Ssend( &b[0], 
	//(e - b)*sizeof(iterator_traits<Ptr>::value_type), 
	(e-b)*sizeof(int *),
	MPI_BYTE,
	dest, 
	tag,
	_mpi_comm);
}

template < class Ptr >
void Communicator::recv( Ptr b, Ptr e, int dest, int tag = 0)
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
void Communicator::send( X *b, X *e, int dest, int tag = 0)
{
	MPI_Send( b, 
	(e - b)*sizeof(X), 
	MPI_BYTE,
	dest, 
	tag,
	_mpi_comm);
}


template < class X >
void Communicator::recv( X *b, X *e, int dest, int tag = 0)
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

   /// member access

MPI_Status Communicator::get_stat(){
   return stat;
}

}// namespace LinBox
#endif // __MPICPP_INL_
