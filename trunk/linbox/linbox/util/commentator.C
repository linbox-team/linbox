/* -*- mode: c; style: linux -*- */

/* linbox/src/util/commentator.C
 * Copyright (C) 1999 B. David Saunders,
 *                    Jean-Guillaume Dumas
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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
 *
 * This file implements the C++ interface to commentators (for 
 * providing runtime commentary to the user)
 */

#ifndef __COMMENTATOR_C
#define __COMMENTATOR_C

#include <math.h>

#include "commentator.h"

// -----------------------------------------------------
// Mathematical routines
// -----------------------------------------------------
double nroot (double a, long r, double precision)
{
	long rm = r - 1 ;
	double c1 = double (rm) / double (r), c2 = a / double (r);
	double g = 1, pgr = 1, err = a - 1;
    
	while (err > precision) {
		g = g * c1 + c2 / pgr;
		pgr = pow (g, rm);
		err = a - pgr * g;
		if (err < 0)
			err = -err;

	void MessageClass::setMaxDetailLevel (long level)
	return g;
}

long isnpower (long& l, long a)
{
	long r = 2;
	double g;

	while ((g = nroot (a, r, 0.1)) >= 2) {
		l = (long) floor (g);
		if (g-double (l) > 0.1)
			++l;
		if (pow (l, r) == a)
			return r;
		++r;

	bool MessageClass::checkConfig (list <pair <unsigned long, unsigned long> > &config, unsigned long depth, unsigned long level) 
	return 0;
}
	// Default global commentator
// wrapper for use by Pascal (& C) code.
extern "C" 
Commentator* initializeCommentator (long timing, long homology)
	{ return new Commentator (timing, homology); }

extern "C"
void startActivity (Commentator& C, char* id, char* msg, long msglevel, long msgclass)
	{ C.start (id, msg, msglevel, msgclass); }

extern "C"
void stopActivity (Commentator& C, char* msg, long msglevel, long msgclass)
	{ C.stop (msg, msglevel, msgclass); }

extern "C"
void activityReport (const Commentator& C, char* msg, long msglevel, long msgclass)
	{ C.report (msg, msglevel, msgclass); }

extern "C"
void progressReport (Commentator& C, char* msg, long msglevel, long k, long n)
	{ C.progress (msg, msglevel, k, n); }

extern "C"
long isPrinted (const Commentator& C, long msglevel, long msgclass)
	{ return C.printed (msglevel, msgclass); }

	// Null streambuf
	ostream Commentator::cnull;
