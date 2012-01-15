/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* lb-garbage.h
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

#ifndef __LINBOX_lb_garbage_H
#define __LINBOX_lb_garbage_H


#include <lb-domain-collection.h>
#include <lb-blackbox-collection.h>
#include <lb-vector-collection.h>
#include <lb-element-collection.h>

/***************************************
 * API to delete a domain from ist key *
 ***************************************/
void deleteElement(const EltKey &key);

/*****************************
 * API to collect all domain *
 *****************************/
void collectElement();

/***************************************
 * API to delete a domain from ist key *
 ***************************************/
void deleteDomain(const DomainKey &key);

/*****************************
 * API to collect all domain *
 *****************************/
void collectDomain();

/*****************************************
 * API to delete a blackbox from its key *
 *****************************************/
void deleteBlackbox (const BlackboxKey &key);

/*******************************
 * API to collect all blackbox *
 *******************************/
void collectBlackbox();

/***************************************
 * API to delete a vector from its key *
 ***************************************/
void deleteVector (const VectorKey &key);

/*****************************
 * API to collect all vector *
 *****************************/
void collectVector();

/***********************************************
 * API to collect all data allocated by LinBox *
 **********************************************/
void LinBoxCollect();


#endif
