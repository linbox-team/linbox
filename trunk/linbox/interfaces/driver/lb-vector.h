/* lb-vector.h
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

#ifndef __LINBOX_lb_vector_H
#define __LINBOX_lb_vector_H

#include <lb-domain-collection.h>
#include <lb-vector-collection.h>


/*************************
 * Initializer of Vector *
 *************************/
void UpdateVector();

/*********************************************************
 * API to contruct a n dimensional vector  over a domain *
 *********************************************************/
const VectorKey& createVector(const DomainKey &k, size_t n, const char *name=NULL);


/******************************************
 * API to contruct a vector from a stream *
 ******************************************/
const VectorKey& createVector(const DomainKey &k, std::istream &in, const char *name=NULL);


/**********************************
 * API to copy an existing vector *
 **********************************/
const VectorKey& copyVector(const VectorKey &k);


/******************************************
 * API to get the dimensions of a vector  *
 ******************************************/
size_t getVectorDimension(const VectorKey &key);


/*****************************************
 * API to set a vector with random value *
 *****************************************/
void setVectorAtRandom(const VectorKey &k);


/********************************************
 * API to rebind a vector over a new domain *
 ********************************************/
void rebindVector(const VectorKey &Vkey, const DomainKey &Dkey);


/****************************************
 * API to write a vector over an stream *
 ****************************************/
void writeVector (const VectorKey &key, std::ostream &os);


/*******************************************
 * API to modify the current vector type *
 *******************************************/
void setVector(const char* t);


/*********************************
 * API to write info on a vector *
 *********************************/
void writeVectorInfo(const VectorKey &k, std::ostream& os);


/******************************
 * API to serialize a  vector *
 ******************************/
void SerializeVector (SerialVector &s, const VectorKey &key);

#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
