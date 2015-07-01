/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/cra.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 *
 * Chinese remainder algorithm for integers and rationals
 */

#ifndef __CRA_H
#define __CRA_H

#include <algorithm>

#include "linbox/integer.h"
#include "linbox/rational.h"

namespace LinBox 
{

/** Chinese remainder algorithm for integers
 *
 * Given a vector of residues $x_i$ and a vector of moduli $n_i$, $1\leq i\leq
 * r, reconstruct the integer x that is congruent to $x_i$ modulo $n_i$ for
 * $1\leq i\leq r$.
 *
 * This code is generic with respect to the vector type, so it may be used,
 * e.g., with subvectors to reconstruct polynomial or integer vector
 * coefficients
 *
 * @param res Place to store result
 * @param residues Dense vector of residues
 * @param moduli Dense vector of moduli
 */

template <class Vector>
integer &cra (integer      &res,
	      const Vector &residues,
	      const Vector &moduli)
{
	linbox_check (residues.size () == moduli.size ());

	commentator.start ("Chinese remainder algorithm", "cra", residues.size ());

	integer pi = 1L, pi_m_j, s;
	typename Vector::const_iterator i, j;

	for (i = moduli.begin (); i != moduli.end (); ++i)
		integer::mulin (pi, (unsigned long) *i);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Product of moduli: " << pi << endl;

	res = 0L;

	for (i = residues.begin (), j = moduli.begin (); j != moduli.end (); ++i, ++j) {
		integer::div (pi_m_j, pi, (unsigned long) *j);
		integer::invmod (s, pi_m_j, integer (*j));
		integer::mulin (s, (unsigned long) *i);
		integer::modin (s, (unsigned long) *j);
		integer::mulin (s, pi_m_j);
		integer::addin (res, s);
		commentator.progress ();
	}

	integer::modin (res, pi);

	commentator.stop ("done", NULL, "cra");

	return res;
}								       		

} // namespace LinBox

#endif // __CRA_H
