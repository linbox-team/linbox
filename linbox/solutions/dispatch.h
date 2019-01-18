/*
 * Copyright (C) LinBox
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
 */

#pragma once

#include <linbox/util/mpicpp.h> // Communicator @fixme Rename file?

namespace LinBox {
    /// Let implementation decide what to use.
    struct AutoDispatch {
    };

    /// All sub-computations are done sequentially.
    struct SequentialDispatch {
        bool master() const { return true; }
    };

    /// Use symmetric multiprocessing (Paladin) to do sub-computations.
    struct SmpDispatch {
        // @fixme Could store thread info (is there any?)

        bool master() const { return true; }
    };

    /// Use MPI to distribute sub-computations accross nodes.
    struct DistributedDispatch {
        Communicator* communicator = nullptr;

        bool master() const { return communicator != nullptr && communicator->master(); }
    };

    /// Use MPI then Paladin on each node.
    struct CombinedDispatch : public DistributedDispatch, public SmpDispatch {
        using DistributedDispatch::master;
    };

    /// How to dispatch sub-computations for CRA.
    struct Dispatch {
        using Auto = AutoDispatch;
        using Sequential = SequentialDispatch;
        using Smp = SmpDispatch;
        using Distributed = DistributedDispatch;
        using Combined = CombinedDispatch;
    };
}
