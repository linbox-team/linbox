/* lb-utilities.h
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

#ifndef __LINBOX_lb_utilities_H
#define __LINBOX_lb_utilities_H


#include <linbox/util/error.h>
#include <linbox/util/debug.h>

#include <iostream>
#include <sstream>


/**************************************
 * definition of LinBox's exceptions  *
 **************************************/
typedef LinBox::LinboxError lb_runtime_error;



/************************************************
 * definition of a functor to compare two char* *
 ************************************************/
struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {return strcmp(s1, s2) < 0;}
};


/*****************************************
 * definition of LinBox's object visitor *
 *****************************************/

class LinBoxBaseVisitor {
public:
	virtual ~LinBoxBaseVisitor() {}
};

template <class T>
class LinBoxVisitor : public virtual LinBoxBaseVisitor {
public:
	virtual void visit(const T&) = 0;
};


class LinBoxBaseVisitable {
public:
  virtual ~LinBoxBaseVisitable(){}
  virtual void  Accept(LinBoxBaseVisitor&) = 0;

protected:
	template <class T>
	static void AcceptImpl(const T& visited, LinBoxBaseVisitor &guest){
	  if (LinBoxVisitor<T>* ptr = dynamic_cast<LinBoxVisitor<T>*> (&guest))
		  ptr->visit(visited);
	  else
		  throw lb_runtime_error("LinBox ERROR: no visitor found (check that all the type are defined in their corresponding type list)\n ");
	}
};

#define LINBOX_VISITABLE() \
	virtual void Accept(LinBoxBaseVisitor &guest) \
	{ AcceptImpl(*this, guest);}


/*******************************************
 * definition of LinBox's object type list *
 *******************************************/

template <class T, class U>
struct LinBoxTypelist {
	typedef T Head;
	typedef U Tail;
};


// Dumb type for ending a type list
struct LinBoxDumbType{};

/*************************************************
 * function to append  LinBox's object type list *
 *************************************************/
namespace LinBoxTL{

	template <class List1, class List2>
	struct Append;
	
	template <>
	struct Append <LinBoxDumbType, LinBoxDumbType>{
		typedef LinBoxDumbType Result; 
	};
	
	template <class T> 
	struct Append <LinBoxDumbType, T> {
		typedef LinBoxTypelist<T, LinBoxDumbType> Result;
	};
	
	template <class T, class U>
	struct Append <LinBoxDumbType, LinBoxTypelist<T, U> > {
		typedef LinBoxTypelist<T,U> Result;
	};
	
	template <class Head, class Tail, class T>
	struct Append <LinBoxTypelist<Head, Tail>, T > {
		typedef LinBoxTypelist<Head, typename Append<Tail, T>::Result> Result;
	};
}

#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
