/* lb-vector.C
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

#ifndef __LINBOX_lb_vector_C
#define __LINBOX_lb_vector_C

#include <linbox-config.h>

#include <lb-vector.h>
#include <lb-vector-function.h>

#include <lb-domain-collection.h>


/*\**************************
 * Allocate global variable *
 ****************************/

// global hash table for allocated vector 
VectorTable vector_hashtable;

// global variable for the factory
Vector_Factory linbox_vector;

// global variable for current vector type
const char* current_vector  = default_vector;



/**********************************************************
 * function to add an abstract vector in linbox hashtable *
 **********************************************************/
const VectorKey& addVector(VectorAbstract * v){
	
	std::pair<VectorTable::const_iterator, bool> status;
	status = vector_hashtable.insert(std::pair<VectorKey, VectorAbstract*> (VectorKey(v), v));
	if (status.second)
		return status.first->first;
	else
		throw lb_runtime_error("LinBox ERROR: vector creation failed \n");// throw an exception
}


/*********************************************************
 * API to contruct a n dimensional vector  over a domain *
 *********************************************************/
const VectorKey& createVector(const DomainKey &k, size_t n, const char *name){
	const char *type=name;
	if (type == NULL)
		type= current_vector;

	VectorAbstract *v = linbox_vector.create(type, k, n); 
	return addVector(v);
}



/******************************************
 * API to contruct a vector from a stream *
 ******************************************/
const VectorKey& createVector(const DomainKey &k, std::istream &in, const char *name){
	const char *type=name;
	if (type == NULL)
		type= current_vector;
	
	VectorAbstract *v = linbox_vector.create(type, k, in); 
	return addVector(v);
}


/**********************************
 * API to copy an existing vector *
 **********************************/
const VectorKey& copyVector(const VectorKey &k){
	
	VectorTable::iterator it = vector_hashtable.find(k);
	if (it == vector_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: vector does not exist (copying impossible)\n");
	
	VectorAbstract *v = it->second->clone();
	return addVector(v);       
}


/******************************************
 * API to get the dimensions of  a vector *
 ******************************************/
class VectorDimensionFunctor{
public:
	template<class Vector>
	void operator()(size_t &dim, Vector *V) const{
		dim = V->size();
	}
};

size_t getVectorDimension(const VectorKey &key){
	size_t dim;
	VectorDimensionFunctor Fct;
	VectorFunction::call(dim, key, Fct);
	return dim;
}

/*****************************************
 * API to set a vector with random value *
 *****************************************/
template<class VElement, class DElement>
class RandomVector{
public:
	template<class Vector, class Domain>
	void operator()(Vector *V, Domain *D) {
		throw lb_runtime_error("LinBox ERROR: incompatible type (vector writing impossible)\n");
	}
};


template<class Element>
class RandomVector<Element, Element>{
public:
	template<class Vector, class Domain>
	void operator()(Vector *V, Domain *D) {
		typename Domain::RandIter G(*D);
		for (size_t i=0; i< V->size();++i) 
			G.random((*V)[i]);
	}
};

template<class Vector>
class VectorAtRandomFunctorSpec{
protected:
	Vector *vect;
public:
	VectorAtRandomFunctorSpec(Vector *V) : vect(V) {}
	
	template<class Domain> 
	void operator()(void *, Domain *D) const {
		RandomVector<typename Vector::value_type, typename Domain::Element>()(vect, D);
	}
};

class VectorAtRandomFunctor{
protected:
	const VectorKey key;
public:
	VectorAtRandomFunctor(const VectorKey &k) : key(k) {}
	
	template<class Vector>
	void operator()(void*, Vector *V) const {
		VectorTable::iterator it = vector_hashtable.find(key);
		if ( it == vector_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: invalid vector (set random value impossible)");
		
		VectorAtRandomFunctorSpec<Vector> Fct(V);
		DomainFunction::call(it->second->getDomainKey(), Fct);
	}
};


void setVectorAtRandom(const VectorKey &k){
	VectorAtRandomFunctor Fct(k);
	VectorFunction::call(k, Fct);
}



/********************************************
 * API to rebind a vector over a new domain *
 ********************************************/
void rebindVector(const VectorKey &Vkey, const DomainKey &Dkey){
	VectorTable::iterator it = vector_hashtable.find(Vkey);
	if ( it == vector_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: invalid vector (rebinding to another domain impossible)");
	it->second->rebind(Dkey);
}


/*****************************************
 * API to write a vector over an ostream *
 *****************************************/
template<class Element, class VElement>
class PrintVector{
public:
	template< class Domain, class Poly, class Out>
	void operator()(Domain *D, Poly *poly, Out &os){
		throw lb_runtime_error("LinBox ERROR: incompatible type (vector writing impossible)\n");
	}
};

template<class Element>
class PrintVector<Element, Element> {
public:
	template< class Domain, class Vector, class Out>
	void operator()(Domain *D, Vector *vect, Out &os){
		os<<" [ "; 
		for (size_t i=0; i<vect->size()-1; ++i){
			D->write(os, (*vect)[i])<<" , ";
		}
		D->write(os, (*vect)[vect->size()-1]);
		os<<" ]\n";
	}
};

template<class Vector>
class WriteVectorSpecFunctor{
protected:
	std::ostream     &os;
	Vector         *vect;
public:
	WriteVectorSpecFunctor(std::ostream &o, Vector *V) : os(o), vect(V){}

	template<class Domain>
	void operator()(void *&, Domain *D) const {
		PrintVector<typename Domain::Element, typename Vector::value_type> ()(D, vect, os);
	}
};

class WriteVectorFunctor{
protected:
	std::ostream     &os;
	const VectorKey &key;
public:	
	WriteVectorFunctor(std::ostream &o, const VectorKey &k) : os(o), key(k) {}

	template<class Vector>
	void operator() (void*, Vector *V) const {
	
		VectorTable::iterator it = vector_hashtable.find(key);
		if ( it == vector_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: invalid vector (writing impossible)");

		WriteVectorSpecFunctor<Vector> Fct(os,V);
		DomainFunction::call(it->second->getDomainKey(), Fct);
	}
};


void writeVector (const VectorKey &key, std::ostream &os){
	WriteVectorFunctor Fct(os, key);
	VectorFunction::call(key, Fct);	
}





/*******************************************
 * API to modify the current vector type *
 *******************************************/
void setVector(const char* t){
	current_vector= t;
}

/*********************************
 * API to write info on a vector *
 *********************************/
void writeVectorInfo(const VectorKey &key, std::ostream& os){
	VectorTable::iterator it = vector_hashtable.find(key);
	if ( it == vector_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: invalid vector (writing vector info impossible)");
	os<<it->second->info();
}




/**********************************
 * API to serialize a  vector *
 **********************************/
template <class PElement, class DElement, class Category>
class ToSerialVector {
public:
	template<class Domain, class Vector>
	inline void operator()(SerialVector &s, Domain *D, Vector *poly){
		throw lb_runtime_error("LinBox ERROR: incompatible type (vector writing impossible)\n");
	}
};

template<class Element, class Category>
class ToSerialVector<Element, Element, Category>{
public:
	template<class Domain, class Vector>
	inline void operator()(SerialVector &s, Domain *D, Vector *poly){
		s.type = "integer";
		s.list.resize(poly->size());
		typename Vector::const_iterator    it_P= poly->begin();
		std::vector<LinBox::integer>::iterator it_s= s.list.begin();
		for (; it_P != poly->end(); ++it_P, ++it_s)
			D->convert(*it_s, *it_P);
	}
};

template<class Element>
class ToSerialVector<Element, Element, LinBox::RingCategories::RationalTag>{
public:
	template<class Domain, class Vector>
	inline void operator()(SerialVector &s, Domain *D, Vector *poly){
		s.type = "rational";
		s.list.resize(2*poly->size());
	
		typename Vector::const_iterator    it_P= poly->begin();
		std::vector<LinBox::integer>::iterator it_s= s.list.begin();
		for (; it_P != poly->end(); ++it_P, ++it_s){
			D->get_num(*it_s, *it_P);
			++it_s;
			D->get_den(*it_s, *it_P);
		}
	}
};


template<class Vector>
class SerializeVectorSpecFunctor{
protected:
	Vector     *poly;
public:
	SerializeVectorSpecFunctor(Vector *P) :  poly(P){}

	template<class Domain>
	void operator()(SerialVector &s, Domain *D) const {
		ToSerialVector<typename Vector::value_type, typename Domain::Element, typename LinBox::FieldTraits<Domain>::categoryTag>()(s, D, poly);
	}
};


class SerializeVectorFunctor{
protected:
	const VectorKey &key;
public:	
	SerializeVectorFunctor(const VectorKey &k) :  key(k) {}

	template<class Vector>
	void operator() (SerialVector& s, Vector *P) const {
	
		VectorTable::iterator it = vector_hashtable.find(key);
		if ( it == vector_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: invalid vector (serializing impossible)");

		SerializeVectorSpecFunctor<Vector> Fct(P);
		DomainFunction::call(s, it->second->getDomainKey(), Fct);
	}
};


void  SerializeVector (SerialVector &s, const VectorKey &key) {       
	SerializeVectorFunctor Fct(key);
	VectorFunction::call(s, key, Fct);	
}


#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
