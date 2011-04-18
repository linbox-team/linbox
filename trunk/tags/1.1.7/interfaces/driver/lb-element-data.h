/* lb-element.h
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


#ifndef __LINBOX_lb_element_data_H
#define __LINBOX_lb_element_data_H

#include <lb-domain-data.h>

#include <lb-element-collection.h>
#include <lb-element-abstract.h>
#include <lb-domain-function.h>

extern EltTable element_hashtable;

/**********************************************************
 * Element Envelope to be compliant with Element Abstract *
 **********************************************************/
template<class Element>
class EltEnvelope :  public EltAbstract {
	DomainKey key;
	Element  *ptr;
public:
	
	EltEnvelope(const DomainKey &k, Element* e) : key(k, true), ptr(e) {}
	
	~EltEnvelope() {delete ptr;}

	Element *getElement() const  {return ptr;}			
	
	const DomainKey& getDomainKey() const {	return key;}
	
	EltAbstract* clone() const { return new EltEnvelope<Element>(key, new Element(*ptr));}
};


/*****************************************
 * construction of Element from a Domain *
 *****************************************/
class CreateEltFunctor{
protected:
	const DomainKey &key;
public:
	CreateEltFunctor(const DomainKey &k) : key(k) {}

	template<class Domain>
	void operator()(EltAbstract *& e, Domain *D) const{
		e = new EltEnvelope<typename Domain::Element> (key, new typename Domain::Element());
	}
};

EltAbstract* constructElt(const DomainKey &key){
	EltAbstract *e;
	CreateEltFunctor Fct(key);
	DomainFunction::call(e, key, Fct);
	return e;
}




#endif

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
