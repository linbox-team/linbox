/* lb-maple.h
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


#ifndef __LINBOX_lb_maple_H
#define __LINBOX_lb_maple_H


#include <lb-driver.h>
#include <maplec.h>

extern "C" {

/***********************************
 * Initializer of LinBox interface *
 ***********************************/
        ALGEB lbStart                (MKernelVector kv, ALGEB *argv);
	ALGEB lbStop                 (MKernelVector kv, ALGEB *argv);


/**************************************
 * Information from the LinBox driver *
 **************************************/
	ALGEB lbDataInfo             (MKernelVector kv, ALGEB *argv);


/**************************************
 * Function to create LinBox's object *
 **************************************/
	ALGEB lbCreateElement        (MKernelVector kv, ALGEB *argv);
	ALGEB lbCreateDomain         (MKernelVector kv, ALGEB *argv);
	ALGEB lbCreateBlackbox       (MKernelVector kv, ALGEB *argv);
	ALGEB lbCreateVector         (MKernelVector kv, ALGEB *argv);


/**********************
 * Domain's Functions *
 **********************/
	ALGEB lbSetPrimeField         (MKernelVector kv, ALGEB *argv);
	ALGEB lbSetRationalField      (MKernelVector kv, ALGEB *argv);
	ALGEB lbSetIntegerRing        (MKernelVector kv, ALGEB *argv);

/************************
 * Blackbox's Functions *
 ************************/
	ALGEB lbCopyBlackbox          (MKernelVector kv, ALGEB *argv);
	ALGEB lbGetBlackboxDimension  (MKernelVector kv, ALGEB *argv);
	ALGEB lbSetBlackboxAtRandom   (MKernelVector kv, ALGEB *argv);	
	ALGEB lbRebindBlackbox        (MKernelVector kv, ALGEB *argv);
	ALGEB lbWriteBlackbox         (MKernelVector kv, ALGEB *argv);
	ALGEB lbSetBlackbox           (MKernelVector kv, ALGEB *argv);

/************************
 * Vector's Functions   *
 ************************/
	ALGEB lbCopyVector            (MKernelVector kv, ALGEB *argv);
	ALGEB lbGetVectorDimension    (MKernelVector kv, ALGEB *argv);
	ALGEB lbSetVectorAtRandom     (MKernelVector kv, ALGEB *argv);
	ALGEB lbRebindVector          (MKernelVector kv, ALGEB *argv);
	ALGEB lbWriteVector           (MKernelVector kv, ALGEB *argv);
	ALGEB lbSetVector             (MKernelVector kv, ALGEB *argv);

/**************************
 * Polynomial's Functions *
 **************************/
	ALGEB lbWritePolynomial        (MKernelVector kv, ALGEB *argv);

/******************************
 * LinBox available solutions *
 ******************************/
	ALGEB lbDeterminant            (MKernelVector kv, ALGEB *argv);
	ALGEB lbRank                   (MKernelVector kv, ALGEB *argv);
	ALGEB lbMinpoly                (MKernelVector kv, ALGEB *argv);
	ALGEB lbCharpoly               (MKernelVector kv, ALGEB *argv);
	ALGEB lbSolve                  (MKernelVector kv, ALGEB *argv);


} // end of extern "C"

#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
