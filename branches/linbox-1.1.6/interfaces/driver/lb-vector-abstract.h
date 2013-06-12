/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* lb-vector-abstract.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 du -h tes *
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

#ifndef __LINBOX_LB_VECTOR_ABSTRACT_H
#define __LINBOX_LB_VECTOR_ABSTRACT_H

#include <lb-utilities.h>
#include <lb-domain-collection.h>


/*******************************
 * Abstract object for Vectors *
 *******************************/
class VectorAbstract : public LinBoxBaseVisitable {
public:
	virtual ~VectorAbstract(){}
	virtual VectorAbstract*     clone()                 const =0;
	virtual void*               getPtr()                const =0;
	virtual const DomainKey&    getDomainKey()          const =0;
	virtual const char*         info()                  const =0;
	virtual void                rebind(const DomainKey&)      =0;
	LINBOX_VISITABLE();
};


#endif