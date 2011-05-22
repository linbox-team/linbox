/* lb-interface.h
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


#include <lb-domain.h>
#include <lb-element.h>
#include <lb-vector.h>
#include <lb-blackbox.h>
#include <lb-polynomial.h>
#include <lb-garbage.h>

#include <lb-utilities.h>

#include <lb-det.h>
#include <lb-rank.h> //problem with givaro-extension (UTT type is not always consistent) disabled in the code
#include <lb-minpoly.h>
#include <lb-charpoly.h> 
#include <lb-solve.h>

// overload PreconditionFailed to be a real exception
std::ostringstream PrecondStream;
void initException(){
	LinBox::PreconditionFailed::setErrorStream(PrecondStream);
}


/***********************************
 * Initialization of LinBox driver *
 ***********************************/
void LinBoxInit(){
	initException();
	UpdateDomain();
	UpdateBlackbox();
	UpdateVector();
}

/***************************
 * Close the LinBox Driver *
 ***************************/
void LinBoxEnd(){
	LinBoxCollect();
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
