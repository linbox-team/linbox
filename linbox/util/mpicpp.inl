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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 * ========LICENCE========
 */

#pragma once

#include "./mpicpp.h"

#include "./serialization.h"

namespace LinBox {

    // ----- Constructors

    Communicator::Communicator(int* argc, char*** argv)
        : _comm(MPI_COMM_WORLD)
        , _boss(true)
    {
        MPI_Init(argc, argv);

        MPI_Comm_rank(_comm, &_rank);
        MPI_Comm_size(_comm, &_size);
    }

    Communicator::Communicator(int* argc, char*** argv, ThreadMode threadMode)
        : _comm(MPI_COMM_WORLD)
        , _boss(true)
    {
        int effectiveThreadMode = -1;
        MPI_Init_thread(argc, argv, static_cast<int>(threadMode), &effectiveThreadMode);
        if (effectiveThreadMode != static_cast<int>(threadMode)) {
            std::cerr << "Warning: MPI thread mode cannot be set as required." << std::endl;
        }

        MPI_Comm_rank(_comm, &_rank);
        MPI_Comm_size(_comm, &_size);
    }

    Communicator::Communicator(const Communicator& communicator)
        : _comm(communicator._comm)
        , _status(communicator._status)
        , _size(communicator._size)
        , _rank(communicator._rank)
        , _boss(false)
    {
    }

    Communicator::~Communicator()
    {
        if (_boss) {
            MPI_Finalize();
        }
    }

    // peer to peer communication

    template <class Ptr> void Communicator::send(Ptr b, Ptr e, int dest, int tag)
    {
        MPI_Send(&*b, (e - b) * sizeof(typename Ptr::value_type), MPI_BYTE, dest, tag, _comm);
    }

    template <class Ptr> void Communicator::ssend(Ptr b, Ptr e, int dest, int tag)
    {
        MPI_Ssend(&b[0], (e - b) * sizeof(int*), MPI_BYTE, dest, tag, _comm);
    }

    template <class Ptr> void Communicator::recv(Ptr b, Ptr e, int dest, int tag)
    {
        MPI_Recv(&b[0], (e - b) * sizeof(typename Ptr::value_type), MPI_BYTE, dest, tag, _comm, &_status);
    }

    template <class X> void Communicator::recv(X* b, X* e, int dest, int tag)
    {
        MPI_Recv(b, (e - b) * sizeof(X), MPI_BYTE, dest, tag, _comm, &_status);
    }

    // whole object communication

    template <class T> void Communicator::send(const T& value, int dest)
    {
        std::vector<uint8_t> bytes;
        uint64_t length = serialize(bytes, value);
        MPI_Send(bytes.data(), length, MPI_UINT8_T, dest, 0, _comm);
    }

    template <class T> void Communicator::ssend(const T& value, int dest)
    {
        std::vector<uint8_t> bytes;
        uint64_t length = serialize(bytes, value);
        MPI_Ssend(bytes.data(), length, MPI_UINT8_T, dest, 0, _comm);
    }

    template <class T> void Communicator::recv(T& value, int src)
    {
        int length = 0;
        MPI_Probe(src, 0, _comm, &_status);
        MPI_Get_count(&_status, MPI_UINT8_T, &length);

        std::vector<uint8_t> bytes(length);
        MPI_Recv(bytes.data(), length, MPI_UINT8_T, src, 0, _comm, &_status);
        unserialize(value, bytes);
    }

    template <class T> void Communicator::bcast(T& value, int src)
    {
        uint64_t length = 0;
        std::vector<uint8_t> bytes;

        if (src == _rank) {
            length = serialize(bytes, value);
        }
        MPI_Bcast(&length, 1, MPI_INT64_T, src, _comm);
        if (src != _rank) {
            bytes.resize(length);
        }

        MPI_Bcast(bytes.data(), length, MPI_UINT8_T, src, _comm);
        if (src != _rank) {
            unserialize(value, bytes);
        }
    }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
