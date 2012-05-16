/* lb-domain-functor.h
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

#ifndef __LINBOX_lb_domain_functor_H
#define __LINBOX_lb_domain_functor_H

#include <lb-utilities.h>
#include <lb-domain-data.h>
#include <lb-domain-type.h>

/*************************************************************************
 * Base class for apply functor over a domain (used for code generation) *
 *************************************************************************/
template <class Functor, class Result>
class ApplyDomainFunctorBase  {
protected:
	Result        *res;
	const Functor *fct;

public:
	ApplyDomainFunctorBase(){}

	template<class Domain>
	void apply(const Domain &d){
		(*fct)(*res, d.getDomain());
	}
};


/***************************************
 * Macro for automatic code generation *
 ***************************************/
#define LB_DOMAIN_VISIT(D)				\
void visit(const DomainEnvelope<D> &d){apply(d);}


/*******************************************************************************
 * Generate automatically the visitors for all domains in the domain type list *
 *******************************************************************************/

template <class DList, class Functor, class Result> class LinBoxDomainVisitor;

template <class Functor, class Result> class LinBoxDomainVisitor<LinBoxDumbType, Functor, Result>{};

template <class T, class Functor, class Result>
class LinBoxDomainVisitor<LinBoxTypelist<T, LinBoxDumbType>, Functor, Result>
: public LinBoxVisitor <DomainEnvelope<T> >,
public ApplyDomainFunctorBase<Functor, Result> {
public:
	LB_DOMAIN_VISIT(T);
};

template <class Head, class Tail, class Functor, class Result>
class LinBoxDomainVisitor<LinBoxTypelist<Head, Tail>, Functor, Result >
: public LinBoxVisitor <DomainEnvelope<Head> > ,
public LinBoxDomainVisitor<Tail, Functor, Result> {
public:
	LB_DOMAIN_VISIT(Head);
};


/**************************************************
 * functionalities to apply functor over a domain *
 **************************************************/

template <class Functor, class Result>
class ApplyDomainFunctor : public LinBoxDomainVisitor<DomainList, Functor, Result> {
public:
	ApplyDomainFunctor(Result &r, const Functor &f) {
		this->res=&r;
		this->fct=&f;
	}
};

#endif

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

