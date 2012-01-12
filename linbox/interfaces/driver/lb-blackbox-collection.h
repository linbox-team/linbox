/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* lb-blackbox-collection.h
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
 du -h tes *
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

#ifndef __LINBOX_lb_blackbox_collection_H
#define __LINBOX_lb_blackbox_collection_H

#include <map>
#include <lb-blackbox-abstract.h>


/**************************
 * Collection of Blackbox *
 **************************/

// definition of a key
typedef size_t  BlackboxKey;

// comparison functor on key
struct BlackboxKeyLessThan{ bool operator()(const BlackboxKey& k1, const BlackboxKey &k2) {return (k1 < k2);} };

// definition of a hash table type
typedef std::map<BlackboxKey, BlackboxAbstract*, BlackboxKeyLessThan>   BlackboxTable;






#endif
