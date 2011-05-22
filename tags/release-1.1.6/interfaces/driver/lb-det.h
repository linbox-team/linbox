/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* lb-det.h
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

#ifndef __LINBOX_LB_DET_H
#define __LINBOX_LB_DET_H

#include <lb-element-collection.h>
#include <lb-blackbox-collection.h>


/*******************************************************
 * API for determinant computation                     *
 * determinant is returned through a given element key *
 *******************************************************/
void lb_determinant(const EltKey &Ekey, const BlackboxKey &Bkey, const char *method="hybrid");


/**************************************************
 * API for determinant computation                *
 * determinant is returned through an element key *
 **************************************************/
const EltKey& lb_determinant(const BlackboxKey &key, const char *method="hybrid");


/*************************************************************
 * API to print available method for determinant computation *
 *************************************************************/
const char* lb_determinant_methods();


#endif

