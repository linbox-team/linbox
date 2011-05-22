/* lb-vector-functor.h
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

#ifndef __LINBOX_lb_vector_functor_H
#define __LINBOX_lb_vector_functor_H

#include <lb-utilities.h>
#include <lb-vector-data.h>
#include <lb-vector-type.h>



/*************************************************************************
 * Base class for apply functor over a vector (used for code generation) *
 *************************************************************************/
template <class Functor, class Result>
class ApplyVectorFunctorBase  {
protected:
	Result        *res;
	const Functor *fct;

public:
	ApplyVectorFunctorBase(){}
	
	template<class Vector>
	void apply(const Vector &b){
		b.launch(*res, *fct);
	}
};


/***************************************
 * Macro for automatic code generation *
 ***************************************/
#define LB_VECTOR_VISIT(B)				\
	void visit(const B &d){apply(d);}


/**********************************************************************************
 * Generate automatically the visitors for all Vector in the vector type list *
 **********************************************************************************/

template <class BList, class Functor, class Result> class LinBoxVectorVisitor;

template <class Functor, class Result> class LinBoxVectorVisitor<LinBoxDumbType, Functor, Result>{};

template <class T, class Functor, class Result> 
class LinBoxVectorVisitor<LinBoxTypelist<T, LinBoxDumbType>, Functor, Result> 
	: public LinBoxVisitor <T>, 
	  public ApplyVectorFunctorBase<Functor, Result> {
public:
	LB_VECTOR_VISIT(T);
};

template <class Head, class Tail, class Functor, class Result>
class LinBoxVectorVisitor<LinBoxTypelist<Head, Tail>, Functor, Result > 
	: public LinBoxVisitor <Head> ,
	  public LinBoxVectorVisitor<Tail, Functor, Result> {
	public:
	LB_VECTOR_VISIT(Head);
};


/****************************************************
 * functionalities to apply functor over a vector *
 ****************************************************/

template<class Functor, class Result>
class ApplyVectorFunctor : public LinBoxVectorVisitor<VectorList, Functor, Result> {
public:
	ApplyVectorFunctor(Result &r, const Functor &f) { 
		this->res = &r;
		this->fct = &f;
	}
};

#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
