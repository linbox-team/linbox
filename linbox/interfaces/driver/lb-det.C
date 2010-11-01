/* lb-det.C
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

#ifndef __LINBOX_lb_det_C
#define __LINBOX_lb_det_C

#include <linbox-config.h>
#include <linbox/solutions/det.h>

#include <lb-det.h>
#include <lb-blackbox-function.h>
#include <lb-element.h>
#include <lb-element-data.h>

extern BlackboxTable blackbox_hashtable;
extern DomainTable   domain_hashtable;

/********************************************************
 * list of available method for determinant computation *
 ********************************************************/
const char* lb_determinant_methods();


/***********************
 * Determinant Functor *
 ***********************/

template<typename Method = LinBox::Method::Hybrid>
class DeterminantFunctor{
private:
	Method meth;
public:
	
	DeterminantFunctor(Method m = Method()) : meth(m) {}

	template<class Blackbox>
	void operator() (EltAbstract *&res, Blackbox *B) const {
		typedef typename Blackbox::Field::Element Element;		
		if (Element *d = (dynamic_cast<EltEnvelope<Element>*>(res))->getElement()) 
			LinBox::det(*d, *B, meth);		
		else
			throw lb_runtime_error("LinBox ERROR: incompatible blackbox and element type (determinant computation impossible)");			       
	}
};

/*******************************************************
 * API for determinant computation                     *
 * determinant is returned through a given element key *
 *******************************************************/

void lb_determinant(const EltKey& Ekey, const BlackboxKey& Bkey, const char* method) {
	EltTable::iterator it = element_hashtable.find(Ekey);
	if ( it == element_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: invalid element (determinant computation impossible)");
	
	DeterminantFunctor<> fct;
	BlackboxFunction::call(it->second, Bkey, fct);
}


/**************************************************
 * API for determinant computation                *
 * determinant is returned through an element key *
 **************************************************/

const EltKey& lb_determinant(const BlackboxKey& key, const char *method) {
	BlackboxTable::iterator it = blackbox_hashtable.find(key);
	if (it == blackbox_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: blackbox does not exist (determinant computation impossible)\n");
	
	const DomainKey *d = &(it->second->getDomainKey());
	const EltKey *e = &createElement(*d);
	lb_determinant(*e, key, method);
	return *e;
}




#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
