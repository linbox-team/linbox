/*
 * Copyright(C) LinBox
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

#if !defined(LINBOX_DEFAULT_FAILURE_PROBABILITY)
#define LINBOX_DEFAULT_FAILURE_PROBABILITY 0.001
#endif

#if !defined(LINBOX_DEFAULT_TRIALS_BEFORE_FAILURE)
#define LINBOX_DEFAULT_TRIALS_BEFORE_FAILURE 20
#endif

#if !defined(LINBOX_DEFAULT_BLOCKING_FACTOR)
#define LINBOX_DEFAULT_BLOCKING_FACTOR 16
#endif

#if !defined(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD)
#define LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD 10
#endif

// Used to decide which method to use when using Method::Auto
//      on a Blackbox or Sparse matrix.
#if !defined(LINBOX_USE_BLACKBOX_THRESHOLD)
#define LINBOX_USE_BLACKBOX_THRESHOLD 10000u
#endif
