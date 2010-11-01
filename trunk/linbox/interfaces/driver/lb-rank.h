/* lb-rank.h
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

#ifndef __LINBOX_lb_rank_H
#define __LINBOX_lb_rank_H


#include <lb-blackbox-collection.h>


/**********************************************************
 * API for rank computation                               *
 * rank is returned through a given unsigned long integer *
 **********************************************************/
void lb_rank(unsigned long &res, const BlackboxKey& key);


/*****************************************************
 * API for rank computation                          *
 * rank is returned through an unsigned long integer *
 *****************************************************/
unsigned long lb_rank(const BlackboxKey& key);


/******************************************************
 * API to print available method for rank computation *
 ******************************************************/
const char* lb_rank_methods();


#endif

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
