/* linbox/matrix/sparse-formats.h
 * Copyright (C) 2013 the LinBox group
 *
 * Written by BB <bbboyer@ncsu.edu>,
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_matrix_sparse_formats_H
#define __LINBOX_matrix_sparse_formats_H

#include "linbox/linbox-tags.h"

namespace LinBox {

	/** Exception class for invalid matrix input
	*/
	namespace Exceptions {
		class InvalidMatrixInput {};
	}




	//! Sparse matrix format (memory storage)
	namespace SparseMatrixFormat {
		class ANY {} ;
		// template<typename Row_t>
		class COO         : public ANY {} ; //!< Cordinate
		// template<typename Row_t>
		class COO1        : public ANY {} ; //!< implicit value COO (with only ones, or mones, or..)
		// template<typename Row_t>
		class CSR         : public ANY {} ; //!< compressed row
		// template<typename Row_t>
		class CSR1        : public ANY {} ; //!< implicit value CSR (with only ones, or mones, or..)
		// template<typename Row_t>
		class ELL         : public ANY {} ; //!< ellpack
		// template<typename Row_t>
		class ELL_R       : public ANY {} ; //!< ellpack fixed row
		// template<typename Row_t>
		class ELL_R1      : public ANY {} ; // ELL_R with only ones (or mones, or..)
		class DIA         : public ANY {} ; //!< Diagonal
		class BCSR        : public ANY {} ; //!< Block CSR
		class HYB         : public ANY {} ; //!< hybrid
		class TPL         : public ANY {} ; //!< vector of triples
		class TPL_omp     : public ANY {} ; //!< triplesbb for openmp
		class LIL         : public ANY {} ; //!< vector of pairs

		// the old sparse matrix reps.
		// class VVP : public ANY {} ; // vector of vector of pairs
		// class VPV : public ANY {} ; // vector of pair of vectors
		// class VMap : public ANY {} ; // vector of index to value maps.
		// template<class Row_t>
		class SparseSequence    /* CoP  */ : public ANY {} ;//!< vector/list of pairs (Container of Pairs). SparseSequence.
		// template<class Row_t>
		class SparseAssociative /* CoM  */ : public ANY {} ;//!< vector/list of pairs (Container of Maps). SparseAssociative.
		// template<class Row_t>
		class SparseParallel    /* PoC  */ : public ANY {} ;//!< pair of vector/list (Pair of Containers). SparseParallel.

	} // SparseMatrixFormat

	//! Sparse matrix format (file storage) @bug use the enum!!!!
	namespace SparseFileFormat {
		class SMS {} ; // JG format
		class CSR {} ; // Raleigh 11 format
		class COO {} ; // more standard coo format
	} // SparseFileFormat
}

#endif // __LINBOX_matrix_sparse_formats_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
