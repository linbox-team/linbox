/* lb-garbage.C
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

#ifndef __LINBOX_lb_garbage_C
#define __LINBOX_lb_garbage_C

#include <lb-garbage.h>
#include <lb-domain-function.h>
#include <lb-element-data.h>
#include <lb-blackbox-function.h>
#include <lb-vector-function.h>

#define __LINBOX_NO_GC_EXCEPTION


/***************************************
 * API to delete a domain from ist key *
 ***************************************/
class DeleteElementFunctor{
protected:
	EltAbstract  *elt;
public:	
	DeleteElementFunctor(EltAbstract *e) : elt(e) {}

	template<class Domain>
	void operator() (void *, Domain *D) const {
		if (EltEnvelope<typename Domain::Element> *ptr = dynamic_cast<EltEnvelope<typename Domain::Element>*>(elt))
			delete ptr;
		else
			throw lb_runtime_error("LinBox ERROR: incompatible domain and element type (freeing element impossible)");
	}
};

void deleteElement(const EltKey &key){
	EltTable::iterator it = element_hashtable.find(key);
	if ( it == element_hashtable.end()){
#ifndef __LINBOX_NO_GC_EXCEPTION
		throw lb_runtime_error("LinBox ERROR: invalid element (freeing impossible)");
#endif
	}
	else {
		DeleteElementFunctor Fct(it->second);
#ifdef __LINBOX_NO_GC_EXCEPTION
		try {
#endif
		DomainFunction::call(it->second->getDomainKey(), Fct);	
#ifdef __LINBOX_NO_GC_EXCEPTION
		} catch(lb_runtime_error &t){std::cout<<"LinBox exception catched\n"; exit(0);}
#endif
		element_hashtable.erase(it);
	}
}



/*****************************
 * API to collect all domain *
 *****************************/
void collectElement(){
	EltTable::iterator it = element_hashtable.begin();
	for (;it != element_hashtable.end(); ++it){
		delete it->second;
		element_hashtable.erase(it);
	}
}


/***************************************
 * API to delete a domain from its key *
 ***************************************/
void deleteDomain(const DomainKey &key) {
	if (key.free()){
		DomainTable::iterator it= domain_hashtable.find(key);
		delete it->second;
		domain_hashtable.erase(it);	
	}
	else{
		key.dispose();
	}
}

/*****************************
 * API to collect all domain *
 *****************************/
void collectDomain(){
	DomainTable::iterator it= domain_hashtable.begin();
	for (; it != domain_hashtable.end(); ++it){
		delete it->second;
		domain_hashtable.erase(it);
	}
}


/*****************************************
 * API to delete a blackbox from its key *
 *****************************************/
class DeleteBlackboxFunctor{
public:	
	template<class Blackbox>
	void operator() (void *, Blackbox *B) const {
		delete B;		
	}
};

void deleteBlackbox (const BlackboxKey &key){
	BlackboxTable::iterator it = blackbox_hashtable.find(key);
	
	if ( it == blackbox_hashtable.end()){
#ifndef __LINBOX_NO_GC_EXCEPTION
		throw lb_runtime_error("LinBox ERROR: invalid blackbox (freeing impossible)\n");
#else
		std::cout<<"LinBox ERROR: invalid  blackbox (freeing impossible)\n";
#endif
	}
	else {
		DeleteBlackboxFunctor Fct;
#ifdef __LINBOX_NO_GC_EXCEPTION
		try {
#endif
			BlackboxFunction::call(key, Fct);
#ifdef __LINBOX_NO_GC_EXCEPTION	
		} catch (lb_runtime_error &t) {std::cout<<"LinBox exception catched: "<<t<<"\n"; exit(0);}
#endif
		delete it->second;	
		blackbox_hashtable.erase(it);		
	}
}


/*******************************
 * API to collect all blackbox *
 *******************************/
void collectBlackbox(){
	DeleteBlackboxFunctor Fct;
	BlackboxTable::iterator it= blackbox_hashtable.begin();
	for (; it != blackbox_hashtable.end(); ++it){std::cout<<"bb to delete: "<<it->second->info()<<"\n";
#ifdef __LINBOX_NO_GC_EXCEPTION
		try{
#endif
		BlackboxFunction::call(*it, Fct);
#ifdef __LINBOX_NO_GC_EXCEPTION
		} catch (lb_runtime_error &t) {std::cout<<"LinBox exception catched\n"; exit(0);}
#endif
		delete it->second;
		blackbox_hashtable.erase(it);
	}
}


/***************************************
 * API to delete a vector from its key *
 ***************************************/

class DeleteVectorFunctor{
public:	
	template<class Vector>
	void operator() (void*, Vector *V) const {
		delete V;		
	}
};


void deleteVector (const VectorKey &key){
	VectorTable::iterator it = vector_hashtable.find(key);
	if ( it == vector_hashtable.end()){
#ifndef __LINBOX_NO_GC_EXCEPTION
		throw lb_runtime_error("LinBox ERROR: invalid vector (freeing impossible)");
#endif
	}
	else {
		DeleteVectorFunctor Fct;
#ifdef __LINBOX_NO_GC_EXCEPTION
		try{
#endif
			VectorFunction::call(key, Fct);	
#ifdef __LINBOX_NO_GC_EXCEPTION
		} catch(lb_runtime_error &t) {exit(0);}
#endif
		delete it->second;
		vector_hashtable.erase(it);	
	}
}


/*****************************
 * API to collect all vector *
 *****************************/
void collectVector(){
	DeleteVectorFunctor Fct;
	VectorTable::iterator it= vector_hashtable.begin();
	for (; it != vector_hashtable.end(); ++it){
#ifdef __LINBOX_NO_GC_EXCEPTION
		try{
#endif
			VectorFunction::call(*it, Fct);		
#ifdef __LINBOX_NO_GC_EXCEPTION
		} catch (lb_runtime_error &t) {std::cout<<"LinBox exception catched\n"; exit(0);}
#endif
		delete it->second;
		vector_hashtable.erase(it);
	}
}


/***********************************************
 * API to collect all data allocated by LinBox *
 **********************************************/
void LinBoxCollect(){
	collectVector();   
	collectBlackbox();  
	collectElement();   
	collectDomain();    
} 


#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
