/* lb-vector-data.h
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

#ifndef __LINBOX_lb_vector_data_H
#define __LINBOX_lb_vector_data_H

#include <vector>
#include <linbox/field/hom.h>
#include <map>
#include <utility>

#include <lb-utilities.h>
#include <lb-vector-collection.h>
#include <lb-vector-abstract.h>
#include <lb-domain-function.h>

extern VectorTable vector_hashtable;

/***************************************************
 * Functor to determine domain in Abstract Vectors *
 ***************************************************/

template<template<class , class> class Vector, class Functor, template <class> class Alloc=std::allocator>
class VectorSpecFunctor{
	const Functor &fct;
	void *ptr;
public:
	VectorSpecFunctor(const Functor &f, void *p) : fct(f), ptr(p) {}
	
	template<class Domain, class Result>
	void  operator() (Result& res, Domain *d) const {
		fct(res, static_cast<Vector<typename Domain::Element, Alloc<typename Domain::Element> >*> (ptr));
	} 
};

/*****************************************
 * Factory to construct Abstract Vectors *
 *****************************************/

class Vector_Factory {
public:
	typedef VectorAbstract* (*createVector_1_CallBack)(const DomainKey &, size_t, const char*);
	typedef VectorAbstract* (*createVector_2_CallBack)(const DomainKey &, std::istream&, const char*);

	typedef std::map<const char*, std::pair<createVector_1_CallBack, createVector_2_CallBack > , ltstr> CallBackMap;
	
	bool add(const char *name, std::pair<createVector_1_CallBack, createVector_2_CallBack> createD){
		return _callback.insert(CallBackMap::value_type(name, createD)).second;
	}
	
	bool remove(const char *name){
		return _callback.erase(name) == 1;
	}
	
	VectorAbstract* create(const char *name, const DomainKey &k, size_t n){
		CallBackMap::iterator it = _callback.find(name);
		if (it != _callback.end()){
			return it->second.first(k, n, name);
		}
		else {
			std::string mes("LinBox ERROR: you are trying to construct a non defined vector << ");
			mes+= std::string(name);
			mes+= std::string(" >>\n");
			throw lb_runtime_error(mes.c_str());// throw an exception
		}
	}

	VectorAbstract* create(const char *name, const DomainKey &k, std::istream &is){
		CallBackMap::iterator it = _callback.find(name);
		if (it != _callback.end()){ return it->second.second(k, is, name); }
		else {
			std::string mes("LinBox ERROR: you are trying to construct a non defined vector << ");
			mes+= std::string(name);
			mes+= std::string(" >>\n");
			throw lb_runtime_error(mes.c_str());// throw an exception
		}		
		
	}		

	size_t size() { return _callback.size(); }

private:	
	CallBackMap _callback;
};



/***************************
 * Functor to copy Vectors *
 ***************************/
class CopyVectorFunctor {
public:
	template<class Vector>
	void operator()(void *&res, Vector *V) const {
		res= new Vector(*V);
	}
};


/******************************
 * Functors to rebind Vectors *
 ******************************/

template<template<class, class> class Vector, template <class> class Alloc=std::allocator>
class RebindVectorFunctor{
	void            *&ptr;
public:
	RebindVectorFunctor(void *&p) : ptr(p) {}

	template<class Domain>
	void operator()(const DomainKey &key, Domain *D) const {		
		RebindVectorFunctor fct(ptr);
		DomainFunction::call(*D, key, fct);
	}

	
	template<class DomainSource, class DomainTarget>
	void operator()(DomainSource &res, DomainTarget *D) const {
		Vector<typename DomainSource::Element,Alloc<typename DomainSource::Element> > *v_source= static_cast<Vector<typename DomainSource::Element,Alloc<typename DomainSource::Element> > * >  (ptr);				
	Vector<typename DomainTarget::Element,Alloc<typename DomainTarget::Element> > *v_target= new Vector<typename DomainTarget::Element,Alloc<typename DomainTarget::Element> >(v_source->size());

		LinBox::Hom<DomainSource, DomainTarget> hom(res, *D);
		for (size_t i=0;i<v_source->size();++i)
			hom.image((*v_target)[i], (*v_source)[i]);
		
		delete v_source;
		ptr = v_target;
	}
};

/********************************************************
 * Vector Envelope to be compliant with Vector Abstract *
 ********************************************************/

template<template<class Element, class Alloc=std::allocator<Element> > class Vector> 
class VectorEnvelope : public VectorAbstract {
protected:
	void         *ptr;
	DomainKey     key;
	const char* _info;
public:
	VectorEnvelope(void* p, const DomainKey &k, const char* info) : ptr(p), key(k, true), _info(info) {}	

	~VectorEnvelope() {}
	
	LINBOX_VISITABLE();

	template<class Functor, class Result>
	void  launch (Result &res, const Functor &fct) const {
		VectorSpecFunctor<Vector, Functor> vs(fct, ptr);
		DomainFunction::call(res, key, vs);
	}
	
	VectorAbstract* clone() const {
		CopyVectorFunctor Fct;
		void *v;
		launch(v, Fct);
		return new VectorEnvelope<Vector>(v, key, _info);
	}

	void * getPtr() const { return ptr;}

	const DomainKey& getDomainKey() const {return key;}
	
	const char* info() const {
		std::string msg= "[ LinBox Vector (storage = ";
		msg+= std::string(_info);
		msg+= std::string(", domain = [LinBox Domain (type = ");
		msg+= std::string(key.Type());
		msg+= std::string(", charact = ");
		msg+= std::string(key.Characteristic());
		msg+= std::string(")] )]\n");
		return msg.c_str();
	}

	void rebind(const DomainKey &k) {
		RebindVectorFunctor<Vector> Fct(ptr);		
		DomainFunction::call(k, key, Fct);	
		key = k;
		key.set_autogc();
	}	
};



/*********************************
 * Functors to construct Vectors *
 *********************************/

template<template<class, class> class Vector, template <class> class Alloc=std::allocator>
class CreateVectorFunctor{
	size_t &_dim;
public:
	CreateVectorFunctor( size_t &n) : _dim(n) {}

	template<class Domain>
	void operator()(void *&res, Domain *D) const {		
		typename Domain::Element zero;
		D->init(zero, 0UL);
		res = new Vector<typename Domain::Element, Alloc<typename Domain::Element> >(_dim, zero);	
	}
};

template<template<class,class> class Vector, template <class> class Alloc=std::allocator>
class CreateVectorFromStreamFunctor{	
	std::istream &in;
public:
	CreateVectorFromStreamFunctor(std::istream &i) : in(i) {}

	template<class Domain>
	void operator()(void *&res, Domain *D) const {
		size_t n;
		LinBox::integer tmp;
		in>>n;
		Vector<typename Domain::Element, Alloc<typename Domain::Element> > * v = new Vector<typename Domain::Element, Alloc<typename Domain::Element> >(n);
		typename Vector<typename Domain::Element,Alloc<typename Domain::Element> >::iterator it = v->begin();
		for (; it != v->end(); ++it){
			in>>tmp;
			D->init(*it, tmp);
		}
		res =v;
	}
};

/******************************************************
 * Vector construction function used in the Factory *
 ******************************************************/

template<template<class T, class Allocator=std::allocator<T> > class Vector>
VectorAbstract* constructVector_from_size(const DomainKey &k, size_t n, const char* info){
	CreateVectorFunctor<Vector> fct(n);
	void *bb;
	DomainFunction::call(bb, k, fct);
	VectorEnvelope<Vector> *bbe = new VectorEnvelope<Vector> (bb, k, info);
	return bbe;
}

template<template<class T, class Allocator=std::allocator<T> > class Vector>
VectorAbstract* constructVector_from_stream (const DomainKey &k, std::istream &in, const char *info){
	CreateVectorFromStreamFunctor<Vector> fct(in);
	void *v;
	DomainFunction::call(v, k, fct);
	VectorEnvelope<Vector> *ve = new VectorEnvelope<Vector> (v, k, info);
	return ve;
}


#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
