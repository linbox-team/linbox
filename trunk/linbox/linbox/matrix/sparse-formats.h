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



	// add here the sparse matrix types

	namespace SparseMatrix2Format {
		class ANY {} ;
		class COO : public ANY {} ; // Cordinate
		// template<typename Row_t>
		class CSR : public ANY {} ; // compressed row
		class ELL : public ANY {} ; // ellpack
		class ELL_R : public ANY {} ; // ellpack
		class HYB : public ANY {} ; // hybrid
		// class TPL : public ANY {} ; // triples
	} // SparseMatrix2Format

	namespace SparseFileFormat {
		class SMS {} ; // JG format
		class CSR {} ; // Raleigh 11 format
		class COO {} ; // more standard coo format
	} // SparseFileFormat
}

#endif

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
