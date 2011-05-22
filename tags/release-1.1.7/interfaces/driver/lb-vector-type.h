/* lb-vector-type.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
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

#ifndef __LINBOX_lb_vector_type_H
#define __LINBOX_lb_vector_type_H


/**************************************
 * Define the list of all Vector Type *
 **************************************/

// (NEED TO USE ENVELOPE TO DEFINE A CONCRETE TYPE)
typedef LinBoxTypelist < VectorEnvelope< std::vector > , LinBoxDumbType> VL1;


// define the vector typelist
typedef VL1 VectorList;


/*******************************************
 * Update the Factory with all vector type *
 *******************************************/
extern Vector_Factory linbox_vector;

void UpdateVector() {
	linbox_vector.add("linbox_dense", Vector_Factory::CallBackMap::value_type::second_type( constructVector_from_size<std::vector>, 
												constructVector_from_stream<std::vector> ));
}



/***************************
 * Default type for vector *
 ***************************/

// definition of the default type vector
#define default_vector  "linbox_dense"



#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
