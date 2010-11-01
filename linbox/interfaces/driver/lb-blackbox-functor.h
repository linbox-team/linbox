/* lb-blackbox-functor.h
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

#ifndef __LINBOX_lb_blackbox_functor_H
#define __LINBOX_lb_blackbox_functor_H

#include <lb-utilities.h>
#include <lb-blackbox-data.h>
#include <lb-blackbox-type.h>


/*************************************************************************
 * Base class for apply functor over a blackbox (used for code generation) *
 *************************************************************************/
template <class Functor, class Result>
class ApplyBlackboxFunctorBase  {
protected:
	Result        *res;
	const Functor *fct;

public:
	ApplyBlackboxFunctorBase(){}
	
	template<class Blackbox>
	void apply(const Blackbox &b){
		b.launch(*res, *fct);
	}
};


/***************************************
 * Macro for automatic code generation *
 ***************************************/
#define LB_BLACKBOX_VISIT(B)				\
	void visit(const B &d){apply(d);}


/**********************************************************************************
 * Generate automatically the visitors for all Blackbox in the blackbox type list *
 **********************************************************************************/

template <class BList, class Functor, class Result> class LinBoxBlackboxVisitor;

template <class Functor, class Result> class LinBoxBlackboxVisitor<LinBoxDumbType, Functor, Result>{};

template <class T, class Functor, class Result> 
class LinBoxBlackboxVisitor<LinBoxTypelist<T, LinBoxDumbType>, Functor, Result> 
	: public LinBoxVisitor <T>, 
	  public ApplyBlackboxFunctorBase<Functor, Result> {
public:
	LB_BLACKBOX_VISIT(T);
};

template <class Head, class Tail, class Functor, class Result>
class LinBoxBlackboxVisitor<LinBoxTypelist<Head, Tail>, Functor, Result > 
	: public LinBoxVisitor <Head> ,
	  public LinBoxBlackboxVisitor<Tail, Functor, Result> {
	public:
	LB_BLACKBOX_VISIT(Head);
};


/****************************************************
 * functionalities to apply functor over a blackbox *
 ****************************************************/

template<class Functor, class Result>
class ApplyBlackboxFunctor : public LinBoxBlackboxVisitor<BlackboxList, Functor, Result> {
public:
	ApplyBlackboxFunctor(Result &r, const Functor &f) { 
		this->res = &r;
		this->fct = &f;
	}
};

#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
