/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/trace.h
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * See COPYING for license information.
 */

#ifndef __TRACE_H
#define __TRACE_H

#include <vector>
//#include <algorithm>

// must fix this list...
//#include "linbox/algorithms/wiedemann.h"
//#include "linbox/algorithms/lanczos.h"
//#include "linbox/algorithms/block-lanczos.h"
#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{


/* for trace we actually use only the blackbox method and local defs for sparse,
 dense, and a few other BB types.
*/
/** \brief sum of eigenvalues
 
Also sum of diagonal entries.
This is the generic one.

testing multiple doc comments.
*/
template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A)
{ return trace(t, A, Method::Hybrid()); }

// Any BBs that offer a local trace can specialize the BB class for the Hybrid method.
/** \brief our best guess

Hybrid method will choose based on matrix size and type
*/

template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A, 
		const Method::Hybrid& m)
{ return trace(t, A, Method::Blackbox(m.specifier())); }

/** \brief our elimination (a fake in this case)

Elimination method will go to blackbox.
*/
template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A, const Method::Elimination& m)
{ return trace(t, A, Method::Blackbox(m.specifier())); 
}

/*
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * See COPYING for license information.
 */

	/** Compute the trace of a linear operator A, represented as a black
	 * box. This class is parameterized by the black box type so that it can
	 * be specialized for different black boxes.
	 */

	template <class Blackbox>
	typename Blackbox::Field::Element &trace (typename Blackbox::Field::Element &res,
					const Blackbox          &A, 
					const Method::Blackbox& m)
	{

		typedef typename Blackbox::Field Field;
		typedef std::vector<typename Field::Element> Vector;
		Vector v, w;
		Field F = A.field();
		StandardBasisStream<Field, Vector> stream (F, A.coldim ());

		VectorWrapper::ensureDim (v, A.coldim ());
		VectorWrapper::ensureDim (w, A.rowdim ());

		F.init (res, 0);

		while (stream) {
			stream >> v;
			A.apply (w, v);
			F.addin (res, VectorWrapper::constRef<Field, Vector> (w, stream.pos () - 1));
		}

		return res;
	}
}

#endif // __TRACE_H
