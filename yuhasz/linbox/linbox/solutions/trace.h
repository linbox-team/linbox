/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/trace.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __TRACE_H
#define __TRACE_H

#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** Compute the trace of a linear operator A, represented as a black
	 * box. This class is parameterized by the black box type so that it can
	 * be specialized for different black boxes.
	 */

	template <class Vector, class Field, class Blackbox>
	typename Field::Element &trace (typename Field::Element &res,
					const Blackbox          &A,
					const Field             &F) 
	{

		Vector v, w;
		StandardBasisStream<Field, Vector> stream (F, A.coldim ());

		linbox_check (A.rowdim () == A.coldim ());

		VectorWrapper::ensureDim (v, A.coldim ());
		VectorWrapper::ensureDim (w, A.coldim ());

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
