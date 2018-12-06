/* lb-solve.C
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#include "linbox/solutions/solve.h"

#include "lb-solve.h"
#include "lb-domain-function.h"
#include "lb-blackbox-function.h"
#include "lb-vector-function.h"
#include "lb-blackbox.h"
#include "lb-vector.h"
#include "lb-domain.h"
#include "lb-garbage.h"


extern BlackboxTable blackbox_hashtable;
extern VectorTable   vector_hashtable;
extern const char* current_rational_field;

/****************************************************************
 * Functor interface to launch LinBox solving                   *
 * deal with different API and avoid different type compilation *
 ****************************************************************/

// generic version to handle different type compilation
template<class ResultField, class BlackboxField, class VectorField, class InDomainCategory, class OutDomainCategory>
class LaunchSolveFunctor {
public:
	template<class Result, class Blackbox, class Vector, class Method>
	inline void  operator()(Result &s, Blackbox *B, Vector *v, const Method &m) const {
		throw lb_runtime_error("LinBox ERROR: incompatible domain in solving");// throw an exception for incompatible data type
	}
};

// specialization to launch LinBox solving using standard API
template<class SameField, class SameDomainCategory>
class LaunchSolveFunctor<SameField, SameField, SameField, SameDomainCategory, SameDomainCategory>{
public:
	template<class Result, class Blackbox, class Vector, class Method>
	inline void  operator()(Result &s, Blackbox *B, Vector *v, const Method &m) const {
		LinBox::solve(s, *B, *v, m);
	}
};


// specialization to launch LinBox solving using integer solving API (output is a GMPRational)
template<class SameField>
class LaunchSolveFunctor<SameField, SameField, SameField, LinBox::RingCategories::IntegerTag, LinBox::RingCategories::IntegerTag >{
public:
	template<class Result, class Blackbox, class Vector, class Method>
	inline void  operator()(Result &s, Blackbox *B, Vector *v, const Method &m) const {
		//LinBox::solve(s, *B, *v, m);
		//not yet handled
		throw lb_runtime_error("LinBox ERROR: integer system solving with diophantine result is not yet handled");
	}
};


// specialization to launch LinBox solving using integer solving API (output is a GMPRational)
template<class QField, class Field>
class LaunchSolveFunctor< QField, Field, Field, LinBox::RingCategories::IntegerTag, LinBox::RingCategories::RationalTag >{
public:
	template<class Result, class Blackbox, class Vector, class Method>
	inline void  operator()(Result &s, Blackbox *B, Vector *v, const Method &m) const {
		LinBox::solve(s, *B, *v, m);
	}
};


/**********************************************************************
 * Partial Solving linear system Functor according to a blackbox type *
 **********************************************************************/

template<class Blackbox, class Method>
class SolvePartialFunctor{
protected:
	Blackbox         *_BB;
	const Method   &meth;
public:

	SolvePartialFunctor (Blackbox *B, const Method &m) : _BB(B), meth(m) {}

	template<class Vector>
	void operator()(const VectorKey& Vkey, Vector *res) const {
		VectorTable::iterator it = vector_hashtable.find(Vkey);
		if (it == vector_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: right hand side vector does not exist (solving impossible)\n");

		SolvePartialFunctor<Blackbox, Method> fct(_BB, meth);
		VectorFunction::call(*res, Vkey, fct);
	}

      	template<class Vector, class Result>
	void operator() (Result &res, Vector *v) const {
		typedef typename Blackbox::Field  BField;
		typedef typename Vector::Field    VField;
		typedef typename Result::Field    RField;
		typedef typename LinBox::FieldTraits<RField>::categoryTag OutcategoryTag;
                typedef typename LinBox::FieldTraits<BField>::categoryTag IncategoryTag;
		LaunchSolveFunctor<RField, BField, VField, IncategoryTag, OutcategoryTag>()(res, _BB, v, meth);
	}
};



/**********************************************************
 * Modify the vector Result if it is integer computation  *
 * and vector result is not define over a rational domain *
 **********************************************************/
template<class Category>
class MutateVector{
public:
	void operator()(VectorAbstract *v){}
};

template<>
class MutateVector<LinBox::RingCategories::IntegerTag>{
public:
	void operator()(VectorAbstract *v){
		const DomainKey *k= &createDomain(0, current_rational_field);
		v->rebind(*k);
		deleteDomain(*k);
	}
};


class MutateVectorFunctor{
public:
	template<class Domain>
	void operator()(VectorAbstract *&v, Domain *D) const {
	  MutateVector<typename LinBox::FieldTraits<Domain>::categoryTag>()(v);
	}
};


void modifyResultVector(const VectorKey &key){
	VectorTable::iterator it = vector_hashtable.find(key);
	if (it == vector_hashtable.end())
			throw lb_runtime_error("LinBox ERROR: result vector does not exist (solving impossible)\n");

	const DomainKey *Dkey = &(it->second->getDomainKey());
	MutateVectorFunctor Fct;
	DomainFunction::call(it->second, *Dkey, Fct);
}


/*****************************************
 * Generic Solving linear system Functor *
 *****************************************/

template<typename Method = LinBox::Method::Hybrid >
class SolveFunctor{
protected:
	const VectorKey &_Vkey;
	Method meth;
public:
	SolveFunctor(const VectorKey &Vkey, Method m =Method()) : _Vkey(Vkey), meth(m) {}

	template<class Blackbox, class Result>
	inline void operator()(Result &res, Blackbox *B) const {
		SolvePartialFunctor<Blackbox, Method> fct(B, meth);
		VectorFunction::call(res, _Vkey, fct);
	}
};


/***********************************************************
 * API for solving linear systems                          *
 * vector solution  is returned through a given vector key *
 ***********************************************************/

void lb_solve(const VectorKey &res, const BlackboxKey &Bkey, const VectorKey &Vkey) {
	SolveFunctor<> fct(res);
	modifyResultVector(res);
	BlackboxFunction::call(Vkey, Bkey, fct);
}

/*****************************************************
 * API for solving linear systems                    *
 * vector solution  is returned through a vector key *
 *****************************************************/
 
const VectorKey&  lb_solve(const BlackboxKey &Bkey, const VectorKey &Vkey) {

	BlackboxTable::iterator it = blackbox_hashtable.find(Bkey);
	if (it == blackbox_hashtable.end())
		throw lb_runtime_error("LinBox ERROR: blackbox does not exist (solving impossible)\n");
	const DomainKey *Dkey = &it->second->getDomainKey();

	std::pair<size_t,size_t> dim;
	dim = getBlackboxDimension(Bkey);
	size_t coldim = dim.second;

	//const VectorKey *res = &copyVector(Vkey);
	const VectorKey *res = &createVector(*Dkey, coldim);

	lb_solve(*res, Bkey, Vkey);
	return *res;
}




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
