/* lb-element.C
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


#ifndef __LINBOX_lb_element_C
#define __LINBOX_lb_element_C

#include <lb-element.h>
#include <lb-element-data.h>

/********************************
 * Allocate all global variable *
 ********************************/

// global hash table for allocated elements
EltTable element_hashtable;


/*******************************************
 * API to contruct a element over a domain * 
 *******************************************/
const EltKey& createElement(const DomainKey &key) {
	EltAbstract *e = constructElt(key);
	std::pair<EltTable::const_iterator, bool> status;
	status = element_hashtable.insert(std::pair<EltKey, EltAbstract*> (EltKey(e), e));
	if (status.second)
		return status.first->first;
	else
		throw lb_runtime_error("LinBox ERROR: element creation failed \n");// throw an exception
}


/*********************************************
 * API to write an a element over its domain * 
 *********************************************/
class WriteElementFunctor{
protected:
	std::ostream     &os;
	EltAbstract  *elt;
public:	
	WriteElementFunctor(std::ostream &o, EltAbstract* e) : os(o), elt(e) {}

	template<class Domain>
	void operator() (void*, Domain *D) const {	
		if (EltEnvelope<typename Domain::Element> *ptr = dynamic_cast<EltEnvelope<typename Domain::Element>*>(elt))
			D->write(os, *(ptr->getElement()))<<"\n";
		else
			throw lb_runtime_error("LinBox ERROR: incompatible domain and element type (element writing impossible)");
	}
};


void writeElement (const EltKey &key, std::ostream &os){
	EltTable::iterator it = element_hashtable.find(key);
	if ( it == element_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: invalid element (writing impossible)");

	WriteElementFunctor Fct(os, it->second);
	DomainFunction::call(it->second->getDomainKey(), Fct);	
}


/*******************************
 * API to serialize an element *
 *******************************/
template<class Domain, class Category>
void ToSerialElement (SerialElement& s, typename Domain::Element &e, const Domain &D, Category cat){
	s.type ="integer";
	LinBox::integer elt;
	D.convert(elt, e);
	s.list.push_back(elt);
}

template <class Domain>
void ToSerialElement (SerialElement& s, typename Domain::Element &e, const Domain &D, LinBox::RingCategories::RationalTag cat){
	s.type ="rational";
	LinBox::integer num, den;
	D.get_num(num, e);
	D.get_den(den, e);
	s.list.push_back(num);
	s.list.push_back(den);
}

class SerializeElementFunctor{
protected:
	EltAbstract *elt;
public:
	SerializeElementFunctor(EltAbstract *e) : elt(e) {}

	template<class Domain>
	void operator() (SerialElement &s, Domain *D) const {
		if (EltEnvelope<typename Domain::Element> *ptr = dynamic_cast<EltEnvelope<typename Domain::Element>*>(elt)){
			ToSerialElement(s, *(ptr->getElement()), *D, typename LinBox::FieldTraits<Domain>::categoryTag());
		}
		else
			throw lb_runtime_error("LinBox ERROR: incompatible domain and element type (element writing impossible)");
	}
};

void  SerializeElement (SerialElement &s, const EltKey &key) {       
	EltTable::iterator it = element_hashtable.find(key);
	if ( it == element_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: invalid element (serializing impossible)");

	SerializeElementFunctor Fct(it->second);
	DomainFunction::call(s, it->second->getDomainKey(), Fct);	
}

#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
