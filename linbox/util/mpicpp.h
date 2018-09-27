/* Copyright (C) 2018 The LinBox group
 * Updated by Hongguang Zhu <zhuhongguang2014@gmail.com>
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

#ifndef __LINBOX_mpicpp_H
#define __LINBOX_mpicpp_H

#ifndef __LINBOX_HAVE_MPI
typedef int Communicator;
#else

// problem of mpi(ch2) in C++
//#undef SEEK_SET
//#undef SEEK_CUR
//#undef SEEK_END


#include <mpi.h>

namespace LinBox {
    /**
     * MPI-based communicator to send/receive LinBox data (like matrices).
     */
    class Communicator {
    public:
        /**
         * Main (boss) communicator.
         * Calls MPI_Init and MPI_Finalize.
         */
        Communicator(int* argc, char*** argv);

        // Non-boss from already existing communicator.
        Communicator(const Communicator& communicator);
        Communicator(MPI_Comm comm = MPI_COMM_NULL);

        ~Communicator();

        // Accessors
        int size() const;
        int rank() const;
        MPI_Status status() const { return _mpi_status; }
        MPI_Comm comm() const { return _mpi_comm; }

        // peer to peer communication
        template <class Ptr> void send(Ptr begin, Ptr end, int dest, int tag);
        template <class Ptr> void ssend(Ptr begin, Ptr end, int dest, int tag);
        template <class Ptr> void recv(Ptr begin, Ptr end, int dest, int tag);
        template <class X> void recv(X* begin, X* end, int dest, int tag);

        // whole object send and recv
        template <class X> void send(X& b, int dest);
        template <class Field> void send(BlasMatrix<Field>& b, int dest);
        template <class Field> void send(SparseMatrix<Field>& b, int dest);
        template <class Field> void send(BlasVector<Field>& b, int dest);

        template <class X> void ssend(X& b, int dest);
        template <class Field> void ssend(BlasMatrix<Field>& b, int dest);
        template <class Field> void ssend(SparseMatrix<Field>& b, int dest);
        template <class Field> void ssend(BlasVector<Field>& b, int dest);

        template <class X> void recv(X& b, int dest);
        template <class Field> void recv(BlasMatrix<Field>& b, int src);
        template <class Field> void recv(SparseMatrix<Field>& b, int src);
        template <class Field> void recv(BlasVector<Field>& b, int src);

        // collective communication
        template <class X> void bcast(X& b, int src);
        template <class X> void bcast(X* b, X* e, int src);
        template <class Field> void bcast(BlasMatrix<Field>& b, int src);
        template <class Field> void bcast(SparseMatrix<Field>& b, int src);
        template <class Field> void bcast(BlasVector<Field>& b, int src);

    protected:
        MPI_Comm _mpi_comm;       // MPI's handle for the communicator
        MPI_Status _mpi_status;   // status from most recent receive
        MPI_Request _mpi_request; // request from most recent send
        bool _mpi_boss = false;   // Whether it's a MPI initializing communicator

        template <class X> void send_integerVec(X& b, int dest);
        template <class X> void send_integerMat(X& b, int dest);
        template <class X> void send_integerSparseMat(X& b, int dest);

        template <class X> void ssend_integerVec(X& b, int dest);
        template <class X> void ssend_integerMat(X& b, int dest);
        template <class X> void ssend_integerSparseMat(X& b, int dest);

        template <class X> void recv_integerVec(X& b, int src);
        template <class X> void recv_integerMat(X& b, int src);
        template <class X> void recv_integerSparseMat(X& b, int src);

        template <class X> void bcast_integerVec(X& b, int src);
        template <class X> void bcast_integerMat(X& b, int src);
        template <class X> void bcast_integerSparseMat(X& b, int src);
    };
}

#include "mpicpp.inl"

#endif
#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
