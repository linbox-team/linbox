/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* lb-element.h
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


#ifndef __LINBOX_LB_ELEMENT_H
#define __LINBOX_LB_ELEMENT_H


#include <lb-domain-collection.h>
#include <lb-element-collection.h>


/*******************************************
 * API to contruct a element over a domain * 
 *******************************************/
const EltKey& createElement(const DomainKey &key);


/*********************************************
 * API to write an a element over its domain * 
 *********************************************/
void writeElement (const EltKey &key, std::ostream &os);

/*******************************
 * API to serialize an element *
 *******************************/
void SerializeElement (SerialElement &s, const EltKey &key);


#endif
