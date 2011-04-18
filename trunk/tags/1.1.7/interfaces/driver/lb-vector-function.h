/* lb-vector-function.h
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

#ifndef __LINBOX_lb_vector_function_H
#define __LINBOX_lb_vector_function_H


#include <lb-vector-abstract.h>
#include <lb-vector-functor.h> 


/*********************************************************
 * API to launch a generic function over a linbox vector *
 *********************************************************/

class VectorFunction {
public:
	// call a functor over a vector from the hashtable, result is given through 1st parameter
	template<class Functor, class Result>
	static void call(Result &res, const std::pair<const VectorKey, VectorAbstract*> &vector, const Functor &functor){
		ApplyVectorFunctor<Functor, Result> Ap(res, functor);
		(vector.second)->Accept(Ap);	
	}

	// call a functor over a vector from the hashtable, no result
	template<class Functor>
	static void call(const std::pair<const VectorKey, VectorAbstract*> &v, const Functor &f){
		void *dumbresult; 
		call(dumbresult,v,f);
	}
	
	// call a functor over a vector from its key, result is given through 1st parameter
	template<class Functor, class Result>
	static void call(Result &res, const VectorKey &key, const Functor &functor){
		VectorTable::const_iterator it = vector_hashtable.find(key);  
		if (it != vector_hashtable.end())
			VectorFunction::call(res, *it, functor);
		else
			throw lb_runtime_error("LinBox ERROR: use of a non allocated vector\n");// throw an exception
	}
	
	// call a functor over a vector from its key, no result
	template<class Functor>
	static void call(const VectorKey &k, const Functor &f) { 
		void *dumbresult; 
		call(dumbresult,k,f);
	}
};
	
#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
