/* lb-domain-data.h
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

#ifndef __LINBOX_lb_domain_data_H
#define __LINBOX_lb_domain_data_H


#include <linbox-config.h>
#include <linbox/integer.h>
#include <linbox/field/field-traits.h>

#include <map>
#include <utility>

#include <lb-utilities.h>
#include <lb-domain-collection.h>
#include <lb-domain-abstract.h>

#include <fstream>

extern DomainTable domain_hashtable;

/*****************************************
 * Factory to construct Abstract Domains *
 *****************************************/

class Domain_Factory {
public:
	typedef DomainAbstract* (*createDomainCallBack)(const LinBox::integer &);
private:
	typedef std::map<const char*, createDomainCallBack, ltstr> CallBackMap;
	CallBackMap _callback;
public:

	bool add(const char *name, createDomainCallBack createD){
		return _callback.insert(CallBackMap::value_type(name, createD)).second;
	}
	
	bool remove(const char *name){
		return _callback.erase(name) == 1;
	}
	
	DomainAbstract* create(const char *name, const LinBox::integer &p){		
		CallBackMap::iterator it = _callback.find(name);
		if (it != _callback.end()){
			return it->second(p);
		}
		else{
			std::string mes("LinBox ERROR: you are trying to construct a non defined domain << ");
			mes+= std::string(name);
			mes+= std::string(" >>\n");
			mes+= std::string(LinBox::integer(_callback.size()));
			throw lb_runtime_error(mes.c_str());// throw an exception
		}
	}

	size_t size() { return _callback.size(); }
};



/*********************************************************
 * Specific Domain constructor according to the category *
 *********************************************************/
template<class Domain, class Category>
class constructDomainFunctor {
public:
	Domain* operator()(const LinBox::integer &p, Category){
		throw lb_runtime_error("LinBox ERROR: try to construct a Domain of unknown category\n");
	}
};

template<class Domain>
class constructDomainFunctor<Domain, LinBox::RingCategories::ModularTag> {
public:
	Domain* operator()(const LinBox::integer &p, LinBox::RingCategories::ModularTag t) {
		return  new Domain(p);
	}
};


template<class Domain>
class constructDomainFunctor<Domain, LinBox::RingCategories::IntegerTag> {
public:
	Domain* operator()(const LinBox::integer &p, LinBox::RingCategories::IntegerTag t) {
		return  new Domain();
	}
};

template<class Domain>
class constructDomainFunctor<Domain, LinBox::RingCategories::RationalTag> {
public:
	Domain* operator()(const LinBox::integer &p, LinBox::RingCategories::RationalTag t) {
		return  new Domain();
	}
};

/********************************************************
 * Domain Envelope to be compliant with Domain Abstract *
 ********************************************************/
template<class Domain>
class DomainEnvelope :  public DomainAbstract {
	Domain *ptr;

	DomainEnvelope (Domain *D) : ptr(D) {}
public:
	typedef Domain Self_t;
	typedef typename LinBox::FieldTraits<Domain>::categoryTag categoryTag;
	
	DomainEnvelope() {}

	DomainEnvelope(const LinBox::integer &p) {
		ptr = constructDomainFunctor<Domain, categoryTag> ()(p, categoryTag());
	}
	
	~DomainEnvelope() {delete ptr;}

	DomainAbstract* clone() const {
		return new DomainEnvelope<Domain> (new Domain(*ptr));
	}

	LINBOX_VISITABLE();

	Domain *getDomain() const  {return ptr;}			
};


/****************************************************
 * Domain construction function used in the Factory *
 ****************************************************/
template<class Domain>
DomainAbstract* constructDomain(const LinBox::integer &p) {
	//return constructDomainFunctor<Domain,  typename LinBox::FieldTraits<Domain>::categoryTag>()(p, typename LinBox::FieldTraits<Domain>::categoryTag());
	return static_cast<DomainAbstract*>(new DomainEnvelope<Domain>(p));
}





#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
