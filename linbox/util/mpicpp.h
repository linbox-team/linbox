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
namespace LinBox {
    // Dummy declaration when no MPI exists.
    class Communicator {
    public:
        Communicator(int* argc, char*** argv) {}

        inline int size() const { return 1; }
        inline int rank() const { return 0; }
        inline bool master() const { return true; }

        template <class T> inline void send(const T& value, int dest) {}
        template <class T> inline void ssend(const T& value, int dest) {}
        template <class T> inline void recv(T& value, int src) {}
        template <class T> inline void bcast(T& value, int src) {}
    };
}
#else

#include <mpi.h>

namespace LinBox {
    /**
     * MPI-based communicator to send/receive LinBox data (like matrices).
     */
    class Communicator {
    public:
        enum class ThreadMode : int {
            Single = MPI_THREAD_SINGLE,         // Only one thread.
            Funneled = MPI_THREAD_FUNNELED,     // Only main thread will make the MPI calls.
            Serialized = MPI_THREAD_SERIALIZED, // Any thread can make MPI calls, but never concurrently.
            Multiple = MPI_THREAD_MULTIPLE,     // No restriction.
        };

    public:
        /**
         * Main (boss) communicator.
         * Calls MPI_Init and MPI_Finalize.
         */
        Communicator(int* argc, char*** argv);
        Communicator(int* argc, char*** argv, ThreadMode threadMode);

        // Non-boss from already existing communicator.
        Communicator(const Communicator& communicator);

        ~Communicator();

        // Accessors
        int size() const { return _size; }
        int rank() const { return _rank; }
        bool master() const { return _rank == 0; }
        MPI_Status status() const { return _status; }
        MPI_Comm comm() const { return _comm; }

        // peer to peer communication
        template <class Ptr> void send(Ptr begin, Ptr end, int dest, int tag);
        template <class Ptr> void ssend(Ptr begin, Ptr end, int dest, int tag);
        template <class Ptr> void recv(Ptr begin, Ptr end, int dest, int tag);
        template <class X> void recv(X* begin, X* end, int dest, int tag);

        // whole object communication
        template <class T> void send(const T& value, int dest);
        template <class T> void ssend(const T& value, int dest);
        template <class T> void recv(T& value, int src);
        template <class T> void bcast(T& value, int src);

    protected:
        MPI_Comm _comm;       // MPI's handle for the communicator
        MPI_Status _status;   // status from most recent receive
        int _size = 0;
        int _rank = 0;
        bool _boss = false;   // Whether it's a MPI initializing communicator
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
