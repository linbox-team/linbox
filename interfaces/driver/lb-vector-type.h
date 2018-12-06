/* lb-vector-type.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
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

#ifndef __LINBOX_lb_vector_type_H
#define __LINBOX_lb_vector_type_H

#include "linbox/vector/blas-vector.h"

/**************************************
 * Define the list of all Vector Type *
 **************************************/

// (NEED TO USE ENVELOPE TO DEFINE A CONCRETE TYPE)
//typedef LinBoxTypelist < VectorEnvelope< std::vector > , LinBoxDumbType> VL1;
typedef LinBoxTypelist < VectorEnvelope< LinBox::BlasVector > , LinBoxDumbType> VL1;


// define the vector typelist
typedef VL1 VectorList;




#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
