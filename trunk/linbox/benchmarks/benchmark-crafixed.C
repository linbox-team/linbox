/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2011 LinBox
 * Written by BB <brice.boyer@imag.fr>
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

/*! @file benchmarks/benchmark-crafixed.C
 * @ingroup benchmarks
 * @brief Benchmarking fixed CRA routines.
 * Here we make benchmarks for CRT (Chinese Remaindering Theorem/Algorithm) in the following case.
 * Let \f$\mathbf{v}\f$ be a vector of size \f$n\f$. Suppose that we only know  \f$\mathbf{v} \mod p_i\f$ for
 * many primes \f$p_i\f$. We try and reconstruct  \f$\mathbf{v}\f$ from these residues.
 *
 * We benchmark for one vector or \f$m\f$ repetitions on different vectors.
 *
 * We use the implementations in LinBox, Givaro, IML and NTL (if the latter two are available).
 */

#include "benchmarks/benchmark.h"
#include "linbox/util/error.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-balanced.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blas-domain.h"

/*   */

int main(int ac, char** av)
{
	return EXIT_SUCCESS ;
}
