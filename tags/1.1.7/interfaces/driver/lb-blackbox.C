/* lb-blackbox.C
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

#ifndef __LINBOX_lb_blackbox_C
#define __LINBOX_lb_blackbox_C

#include <linbox-config.h>

#include <lb-blackbox.h>
#include <lb-blackbox-function.h>

#include <lb-domain-collection.h>

//#include <linbox/matrix/matrix-category.h>

/********************************
 * Allocate the global variable *
 ********************************/

// global hashtable for the allocated blackbox
BlackboxTable blackbox_hashtable;

// global variable for the blackbox factory
Blackbox_Factory linbox_blackbox;

// global variable for current blackbox type
const char* current_blackbox  = default_blackbox;


/*******************************************************
 * API to contruct a m x n zero blackbox over a domain * 
 *******************************************************/
const BlackboxKey& createBlackbox(const DomainKey &k, size_t m, size_t n, const char* name){
	const char* type = name;
	if (type == NULL)
		type = current_blackbox;
	       
	BlackboxAbstract* bb = linbox_blackbox.create(type, k, m , n);
	return addBlackbox(bb);
}

/**********************************************************
 * API to contruct a blackbox over a domain from a stream *
 **********************************************************/
const BlackboxKey& createBlackbox(const DomainKey &k, std::istream &in, const char *name){
	const char* type = name;
	if (type == NULL)
		type = current_blackbox;

	BlackboxAbstract* bb = linbox_blackbox.create(type, k, in);
	return addBlackbox(bb);
}

/************************************
 * API to copy an existing blackbox *
 ************************************/
const BlackboxKey& copyBlackbox(const BlackboxKey &k){
	
	BlackboxTable::iterator it = blackbox_hashtable.find(k);
	if (it == blackbox_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: blackbox does not exist (copying impossible)\n");
	
	BlackboxAbstract *v = it->second->clone();
	return addBlackbox(v);       
}

/********************************************
 * API to get the dimensions of  a blackbox *
 ********************************************/
class BlackboxDimensionFunctor{
public:
	template<class Blackbox>
	void operator()(BlackboxDimension &dim, Blackbox *B) const{
		dim.first  = B->rowdim();
		dim.second = B->coldim();
	}
};

BlackboxDimension getBlackboxDimension(const BlackboxKey &key){
	std::pair<size_t, size_t> dim;
	BlackboxDimensionFunctor Fct;
	BlackboxFunction::call(dim, key, Fct);
	return dim;
}


/*******************************************
 * API to set a blackbox with random value *
 *******************************************/
template<class VElement, class DElement, class Category>
class RandomBlackbox{
public:
	template<class Blackbox, class Domain>
	void operator()(Blackbox *V, Domain *D) {
		throw lb_runtime_error("LinBox ERROR: incompatible blackbox type and domain type  (set random value  impossible)\n");
	}
};

template<class Element>
class RandomBlackbox<Element, Element, LinBox::MatrixContainerCategory::Blackbox>{
public:
	template<class Blackbox, class Domain>
	void operator()(Blackbox *V, Domain *D) {
		throw lb_runtime_error("LinBox ERROR: impossible to set random value (storage is a real blackbox)\n");
	}
};

template<class Element, class Category>
class RandomBlackbox<Element, Element, Category>{
public:
	template<class Blackbox, class Domain>
	void operator()(Blackbox *B, Domain *D) {
		typename Domain::RandIter G(*D);
		for (size_t i=0; i< B->coldim();++i) 
			for (size_t j=0;j< B->coldim();++j)
				G.random(B->refEntry(i,j));
	}
};

template<class Blackbox>
class BlackboxAtRandomFunctorSpec{
protected:
	Blackbox *vect;
public:
	BlackboxAtRandomFunctorSpec(Blackbox *V) : vect(V) {}
	
	template<class Domain> 
	void operator()(void *, Domain *D) const {
		RandomBlackbox<typename Blackbox::Field::Element, typename Domain::Element, typename LinBox::MatrixContainerTrait<Blackbox>::Type>()(vect, D);
	}
};

class BlackboxAtRandomFunctor{
protected:
	const BlackboxKey key;
public:
	BlackboxAtRandomFunctor(const BlackboxKey &k) : key(k) {}
	
	template<class Blackbox>
	void operator()(void*, Blackbox *V) const {
		BlackboxTable::iterator it = blackbox_hashtable.find(key);
		if ( it == blackbox_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: invalid blackbox (set random value impossible)");
		
		BlackboxAtRandomFunctorSpec<Blackbox> Fct(V);
		DomainFunction::call(it->second->getDomainKey(), Fct);
	}
};


void setBlackboxAtRandom(const BlackboxKey &k){
	BlackboxAtRandomFunctor Fct(k);
	BlackboxFunction::call(k, Fct);
}


/**********************************************
 * API to rebind a blackbox over a new domain *
 **********************************************/
void rebindBlackbox(const BlackboxKey &Vkey, const DomainKey &Dkey){
	BlackboxTable::iterator it = blackbox_hashtable.find(Vkey);
	if ( it == blackbox_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: invalid blackbox (rebinding to another domain impossible)");
	it->second->rebind(Dkey);
}


/*************************************************
 * API to write a blackbox over an output stream *
 *************************************************/
class WriteBlackboxFunctor{
	std::ostream &os;
public:	
	WriteBlackboxFunctor(std::ostream &o) : os(o) {}
	
	template<class Blackbox>
	void operator() (void*, Blackbox *B) const {
		B->write(os);	
	}
};

void writeBlackbox (const BlackboxKey &key,  std::ostream &os){
	WriteBlackboxFunctor Fct(os);
	BlackboxFunction::call(key, Fct);	
}



/*******************************************
 * API to modify the current blackbox type *
 *******************************************/
void setBlackbox(const char* t){
	current_blackbox= t;
}



/************************************
 * API to write info on a blackbox  *
 ************************************/
void writeBlackboxInfo(const BlackboxKey &k, std::ostream& os){
	BlackboxTable::const_iterator it= blackbox_hashtable.find(k);
	if ( it == blackbox_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: invalid blackbox (writing blackbox information impossible)");

	os<<it->second->info();
}


#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
