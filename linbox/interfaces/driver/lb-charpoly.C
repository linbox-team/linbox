/* lb-charpoly.C
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

#ifndef __LINBOX_lb_charpoly_C
#define __LINBOX_lb_charpoly_C


#include <linbox/solutions/charpoly.h>
#ifdef __LINBOX_HAVE_GIVARO
#include "linbox/ring/givaro-polynomial.h"
#endif


#include <lb-charpoly.h>
#include <lb-blackbox-function.h>
#include <lb-vector.h>


extern BlackboxTable blackbox_hashtable;
extern VectorTable   vector_hashtable;

/*************************************
 * Characteristic Polynomial Functor *
 *************************************/

template<typename Method = LinBox::Method::Hybrid>
class CharpolyFunctor{
protected:
	Method meth;
public:
	CharpolyFunctor(Method m= Method()) : meth(m) {}

	template<class Blackbox, class Result>
	void operator() (Result &res, Blackbox *B) const {
		typedef typename Blackbox::Field Field;
		typedef typename Field::Element Element;
#ifdef __LINBOX_HAVE_GIVARO
		// use givpolynomial du to non genericity of charpoly over integer
		typename LinBox::GivPolynomialRing<Field, Dense>::Element pol;
		LinBox::charpoly(pol, *B, meth);
		
		// convert back the result to std::vector
		std::vector<Element> *phi = static_cast<std::vector<Element>*> (res);
		phi->resize(pol.size());
		for (size_t i=0; i< pol.size(); ++i)
			B->field().assign((*phi)[i], pol[i]);
#else
		throw lb_runtime_error("LinBox ERROR: charpoly computation requires Givaro library, computation impossible)\n");
#endif
		
	}
};


/********************************************************************
 * API for characteristic polynomial computation                    *
 * characteristic polynomial is returned through a given vector key *
 ********************************************************************/

void lb_charpoly(const VectorKey &res, const BlackboxKey& key) {
	VectorTable::iterator it = vector_hashtable.find(res);
	if ( it == vector_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: result polynomial does not exist (charpoly computation impossible)\n");
	CharpolyFunctor<> fct;
	void *ret = it->second->getPtr();
	BlackboxFunction::call(ret, key, fct);
}


/**************************************************************
 * API for characteristic polynomial computation              *
 * characteristic polynomial is returned through a vector key *
 **************************************************************/

const VectorKey& lb_charpoly(const BlackboxKey& key) {
	BlackboxTable::iterator it = blackbox_hashtable.find(key);
	if ( it == blackbox_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: blackbox is not defined (charpoly computation impossible)");
		
	const VectorKey *res = & createVector(it->second->getDomainKey(), 0, "linbox_dense");
	lb_charpoly(*res, key);
	return *res;
}


#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
