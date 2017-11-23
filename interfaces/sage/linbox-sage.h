/* linbox-sage.h
 * Copyright (C) 2007 Martin Albrecht
 *               2008 Clement Pernet
 *
 * Written by Martin Albrecht
 *            Clement Pernet <clement.pernet@gmail.com>
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

#ifndef __LINBOX_sage_H
#define __LINBOX_sage_H

#include <stddef.h>

#include <cstdlib>
#include <vector>

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif
/*****************************************************************

  Sparse over Z/nZ

*****************************************************************/


unsigned long linbox_modn_sparse_matrix_rank(unsigned int modulus,
					     size_t numrows,
					     size_t numcols,
					     void *rows,
					     int reorder);

std::vector<unsigned int> linbox_modn_sparse_matrix_solve (unsigned int modulus, size_t numrows,
						     size_t numcols, void *a, void *b,
						     int method);

#endif // __LINBOX_SAGE_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

