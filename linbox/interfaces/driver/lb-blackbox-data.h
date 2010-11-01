/* lb-blackbox-data.h
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

#ifndef __LINBOX_lb_blackbox_data_H
#define __LINBOX_lb_blackbox_data_H

#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/util/matrix-stream.h>
#include <map>
#include <utility>



#include <lb-utilities.h>
#include <lb-blackbox-collection.h>
#include <lb-blackbox-abstract.h>
#include <lb-domain-function.h>
#include <lb-garbage.h>


extern BlackboxTable blackbox_hashtable;

/******************************************************
 * Functor to determine domain in Abstract Blackboxes *
 ******************************************************/

template<template<class Domain> class Blackbox, class Functor>
class BlackboxSpecFunctor{
	const Functor & fct;
	void *ptr;
public:
	BlackboxSpecFunctor(const Functor &f, void *p) : fct(f), ptr(p) {}
	
	template<class Domain, class Result>
	void  operator()(Result &res, Domain *d) const{
		Blackbox<Domain> * b= static_cast<Blackbox<Domain>*> (ptr);
		fct(res, b);
	}
};


/********************************************
 * Factory to construct Abstract Blackboxes *
 ********************************************/

class Blackbox_Factory {
public:
	typedef BlackboxAbstract* (*createBlackbox_1_CallBack)(const DomainKey &, size_t, size_t, const char*);
	typedef BlackboxAbstract* (*createBlackbox_2_CallBack)(const DomainKey &, std::istream&, const char*);
	typedef std::map<const char*, std::pair<createBlackbox_1_CallBack, createBlackbox_2_CallBack > , ltstr> CallBackMap;
	
	bool add(const char *name, std::pair<createBlackbox_1_CallBack, createBlackbox_2_CallBack> createD){
		return _callback.insert(CallBackMap::value_type(name, createD)).second;
	}
	
	bool remove(const char *name){return _callback.erase(name) == 1;}
			
	BlackboxAbstract* create(const char *name, const DomainKey &k, size_t m, size_t n){
		CallBackMap::iterator it = _callback.find(name);
		if (it != _callback.end()){
			return it->second.first(k, m, m,name);
		}
		else {
			std::string mes("LinBox ERROR: you are trying to construct a non defined blackbox << ");
			mes+= std::string(name);
			mes+= std::string(" >>\n");
			throw lb_runtime_error(mes.c_str());// throw an exception
		}
	}

	BlackboxAbstract* create(const char *name, const DomainKey &k, std::istream &is){
		CallBackMap::iterator it = _callback.find(name);
		if (it != _callback.end()){ return it->second.second(k, is, name); }
		else {
			std::string mes("LinBox ERROR: you are trying to construct a non defined blackbox << ");
			mes+= std::string(name);
			mes+= std::string(" >>\n");
			throw lb_runtime_error(mes.c_str());// throw an exception
		}		
		
	}

	size_t size() { return _callback.size(); }

private:	
	CallBackMap _callback;
};





/******************************
 * Functor to copy Blackboxes *
 ******************************/
class CopyBlackboxFunctor {
public:
	template<class Blackbox>
	void operator()(void *&res, Blackbox *B) const {
		res= new Blackbox(static_cast<const Blackbox&>(*B));
	}
};


/*********************************
 * Functors to rebind Blackboxes *
 *********************************/

template<template<class T> class Blackbox>
class RebindBlackboxFunctor{
	void            *&ptr;
public:
	RebindBlackboxFunctor(void *&p) : ptr(p) {}

	template<class Domain>
	void operator()(const DomainKey &key, Domain *D) const {		
		RebindBlackboxFunctor fct(ptr);
		DomainFunction::call(*D, key, fct);
	}

	
	template<class DomainSource, class DomainTarget>
	void operator()(DomainSource &res, DomainTarget *D) const {
		Blackbox<DomainSource> *B_source= static_cast<Blackbox<DomainSource> * >  (ptr);				
		Blackbox<DomainTarget> *B_target;
		typename Blackbox<DomainSource>::template rebind<DomainTarget>()(B_target, *B_source, *D);					
		delete B_source;
		ptr = B_target;
	}
};


/************************************************************
 * Blackbox Envelope to be compliant with Blackbox Abstract *
 ************************************************************/
template<template<class Domain> class Blackbox>
class BlackboxEnvelope : public BlackboxAbstract{
protected:
	void         *ptr;
	DomainKey     key;
	const char* _info;
public:

	BlackboxEnvelope(void *p, const DomainKey &k, const char *info) : ptr(p), key(k, true), _info(info) {}

	~BlackboxEnvelope(){}

	BlackboxAbstract* clone() const {
		CopyBlackboxFunctor Fct;
		void *b;
		launch(b, Fct);
		return new BlackboxEnvelope<Blackbox>(b, key, _info);
	}

	void * getPtr() const { return ptr;}
	
	virtual const DomainKey& getDomainKey() const {return key;}

	LINBOX_VISITABLE();
	
	template<class Functor, class Result>
	void  launch (Result &res, const Functor &fct) const {
		BlackboxSpecFunctor<Blackbox, Functor> bbs(fct, ptr);
		DomainFunction::call(res, key, bbs);
	}

	const char* info() const {
		std::string msg= "[ LinBox Blackbox (storage = ";
		msg+= std::string(_info);
		msg+= std::string(", domain = [LinBox Domain (type = ");
		msg+= std::string(key.Type());
		msg+= std::string(", charact = ");
		msg+= std::string(key.Characteristic());
		msg+= std::string(")] )]\n");
		return msg.c_str();
	}

	void rebind(const DomainKey &k) {
		RebindBlackboxFunctor<Blackbox> Fct(ptr);		
		DomainFunction::call(k, key, Fct);	
		key = k;
		key.set_autogc();
	}	

};



/**********************************
 * Functors to construct Blackbox *
 **********************************/

template<template<class T> class Blackbox>
class CreateBlackboxFunctor{
	size_t &_row;
	size_t &_col;
public:
	CreateBlackboxFunctor(size_t &m, size_t &n) : _row(m), _col(n) {}

	template<class Domain>
	void operator()(void *&res, Domain *D) const {
		res = new Blackbox<Domain>(*D, _row, _col);	
	}
};

template<template<class T> class Blackbox>
class CreateBlackboxFromStreamFunctor{	
	std::istream &in;
public:
	CreateBlackboxFromStreamFunctor(std::istream &i) : in(i) {}

	template<class Domain>
	void operator()(void *&res, Domain *D) const {
		LinBox::MatrixStream<Domain> ms(*D, in);
		res = new Blackbox<Domain>(ms);
	}
};


/******************************************************
 * Blackbox construction function used in the Factory *
 ******************************************************/

template<template<class T> class Blackbox>
BlackboxAbstract* constructBlackbox_from_size(const DomainKey &k, size_t m, size_t n, const char* info){
	CreateBlackboxFunctor<Blackbox> fct(m,n);
	void *bb;
	DomainFunction::call(bb, k, fct);
	BlackboxEnvelope<Blackbox> *bbe = new BlackboxEnvelope<Blackbox> (bb, k, info);
	return bbe;
}


template<template<class T> class Blackbox>
BlackboxAbstract* constructBlackbox_from_stream (const DomainKey &k, std::istream &in, const char *info){
	CreateBlackboxFromStreamFunctor<Blackbox> fct(in);
	void *bb;
	DomainFunction::call(bb, k, fct);
	BlackboxEnvelope<Blackbox> *bbe = new BlackboxEnvelope<Blackbox> (bb, k, info);
	return bbe;
}


/************************************************************
 * Function to add an abstract blackbox in linbox hashtable *
 ************************************************************/
const BlackboxKey& addBlackbox(BlackboxAbstract * v){
	
	std::pair<BlackboxTable::const_iterator, bool> status;
	status = blackbox_hashtable.insert(std::pair<BlackboxKey, BlackboxAbstract*> (BlackboxKey(v), v));
	if (status.second)
		return status.first->first;
	else
		throw lb_runtime_error("LinBox ERROR: blackbox creation failed \n");// throw an exception
}








#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
