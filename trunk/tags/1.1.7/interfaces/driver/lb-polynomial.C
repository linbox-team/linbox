/* lb-polynomial.C
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

#ifndef __LINBOX_lb_polynomial_C
#define __LINBOX_lb_polynomial_C

#include <lb-domain-collection.h>
#include <lb-domain-function.h>
#include <lb-vector-collection.h>
#include <lb-vector-function.h>
#include <lb-polynomial.h>

/***********************************************
 * Polynomial are handled through dense vector *
 ***********************************************/

// nothing to do at this stage


/*************************************************************
 * API to write a polynomial over the standard output stream *
 *************************************************************/

template<class Element, class VElement>
class PrintPoly{
public:
	template< class Domain, class Poly>
	void operator()(Domain *D, Poly *poly, std::ostream &os){
		throw lb_runtime_error("LinBox ERROR: incompatible type (polynomial writing impossible)\n");
	}
};

template<class Element>
class PrintPoly<Element, Element> {
public:
	template< class Domain, class Poly, class Out>
	void operator()(Domain *D, Poly *poly, Out &os){
		size_t deg=poly->size()-1;
		while  (D->isZero((*poly)[deg])) --deg;
		D->write(os,(*poly)[deg])<<".x^"<<deg;
		for (int i= (int)deg-1; i > 1 ; --i) {
			if (!D->isZero((*poly)[i])){os<<" + ";D->write(os, (*poly)[i])<<".x^"<<i;}
		} 
		if (deg > 1) 
			if (!D->isZero((*poly)[1])){os<<" + ";D->write(os,(*poly)[1])<<".x";}
		if (deg > 0)
			if (!D->isZero((*poly)[0])){os<<" + ";D->write(os,(*poly)[0]);}
		os<<"\n";
	}
};

template<class Polynomial>
class WritePolynomialSpecFunctor{
protected:
	std::ostream     &os;
	Polynomial     *poly;
public:
	WritePolynomialSpecFunctor(std::ostream &o, Polynomial *P) : os(o), poly(P){}

	template<class Domain>
	void operator()(void *&, Domain *D) const {
		PrintPoly<typename Domain::Element, typename Polynomial::value_type> ()(D, poly, os);
	}
};


class WritePolynomialFunctor{
protected:
	std::ostream     &os;
	const PolynomialKey &key;
public:	
	WritePolynomialFunctor(std::ostream &o, const PolynomialKey &k) : os(o), key(k) {}

	template<class Polynomial>
	void operator() (void*, Polynomial *P) const {
	
		VectorTable::iterator it = vector_hashtable.find(key);
		if ( it == vector_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: invalid polynomial (writing impossible)");

		WritePolynomialSpecFunctor<Polynomial> Fct(os,P);
		DomainFunction::call(it->second->getDomainKey(), Fct);
	}
};


void writePolynomial (const PolynomialKey &key, std::ostream &os){
	WritePolynomialFunctor Fct(os, key);
	VectorFunction::call(key, Fct);	
}


/**********************************
 * API to serialize a  polynomial *
 **********************************/
template <class PElement, class DElement>
class ToSerialPolynomial {
public:
	template<class Domain, class Polynomial>
	inline void operator()(SerialPolynomial &s, Domain *D, Polynomial *poly){
		throw lb_runtime_error("LinBox ERROR: incompatible type (polynomial writing impossible)\n");
	}
};

template<class Element>
class ToSerialPolynomial<Element, Element>{
public:
	template<class Domain, class Polynomial>
	inline void operator()(SerialPolynomial &s, Domain *D, Polynomial *poly){
		s.type = "integer";
		s.list.resize(poly->size());
		typename Polynomial::const_iterator    it_P= poly->begin();
		std::vector<LinBox::integer>::iterator it_s= s.list.begin();
		for (; it_P != poly->end(); ++it_P, ++it_s)
			D->convert(*it_s, *it_P);
	}
};


template<class Polynomial>
class SerializePolynomialSpecFunctor{
protected:
	Polynomial     *poly;
public:
	SerializePolynomialSpecFunctor(Polynomial *P) :  poly(P){}

	template<class Domain>
	void operator()(SerialPolynomial &s, Domain *D) const {
		ToSerialPolynomial<typename Polynomial::value_type, typename Domain::Element>()(s, D, poly);
	}
};


class SerializePolynomialFunctor{
protected:
	const PolynomialKey &key;
public:	
	SerializePolynomialFunctor(const PolynomialKey &k) :  key(k) {}

	template<class Polynomial>
	void operator() (SerialPolynomial& s, Polynomial *P) const {
	
		VectorTable::iterator it = vector_hashtable.find(key);
		if ( it == vector_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: invalid polynomial (serializing impossible)");

		SerializePolynomialSpecFunctor<Polynomial> Fct(P);
		DomainFunction::call(s, it->second->getDomainKey(), Fct);
	}
};


void  SerializePolynomial (SerialPolynomial &s, const PolynomialKey &key) {       
	SerializePolynomialFunctor Fct(key);
	VectorFunction::call(s, key, Fct);	
}


#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
