/* Copyright (C) 2010 LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_mpicpp_H
#define __LINBOX_mpicpp_H

#ifndef __LINBOX_HAVE_MPI
 typedef int Communicator; 
#else
#include <iterator>

// problem of mpi(ch2) in C++
#undef SEEK_SET 
#undef SEEK_CUR 
#undef SEEK_END 
#include <mpi.h>


namespace LinBox {
/* Idea:  Only use ssend for send.
*/
class Communicator
{   public:
    // constructors and destructor

	// constructor from existing communicator 
	//`Communicator(MPI_Comm comm = MPI_COMM_NULL);
   Communicator(MPI_Comm comm);
	// MPI_initializing constructor
	// When this communicator is destroyed MPI is shut down (finalized).
	Communicator(int* ac, char*** av);

	// copy constructor
	Communicator(const Communicator& D);

	~Communicator();

    // accessors
	int size();

int rank();

	MPI_Status status(); 

	MPI_Comm mpi_communicator();

    // peer to peer communication
	template < class Ptr >
	void send( Ptr b, Ptr e, int dest, int tag);

	template < class Ptr >
	void ssend( Ptr b, Ptr e, int dest, int tag);

	template < class Ptr >
	void recv( Ptr b, Ptr e, int dest, int tag);

	template < class X >
	void send( X *b, X *e, int dest, int tag);

	template < class X >
	void recv( X *b, X *e, int dest, int tag);

	// whole object send and recv
	template < class X >
	void send( X& b, int dest /*, int tag = 0 */);

	template < class X >
	void ssend( X& b, int dest /*, int tag = 0 */);

	template < class X >
	void bsend( X& b, int dest);

	template < class X >
	void recv( X& b, int dest /*, int tag = 0*/);
/*
   template < vector < class X > >
	void send( X& b, int dest );
*/

	template < class X >
	void buffer_attach( X b);

	template < class X >
	int buffer_detach( X &b, int *size);


    // collective communication
	template < class Ptr, class Function_object >
	void reduce( Ptr bloc, Ptr eloc, Ptr bres, Function_object binop, int root);

    // member access
	MPI_Status get_stat();
	 
   protected:
	MPI_Comm _mpi_comm; // MPI's handle for the communicator
	bool _mpi_boss; // true of an MPI initializing communicator
                       // There is at most one initializing communicator.
	MPI_Status stat; // status from most recent receive

};
}// namespace LinBox

#include "mpicpp.inl"
#endif // __LINBOX_HAVE_MPI
#endif // __LINBOX_mpicpp_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
