/* -*- mode: c; style: linux -*- */

/* linbox/tests/test-common.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __TEST_COMMON_H
#define __TEST_COMMON_H

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/archetype.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/integer.h"

enum ArgumentType {
	TYPE_NONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE
};

struct Argument 
{
	char             c;
	char            *example;
	char            *helpString;
	ArgumentType     type;
	void            *data;
};

template <class Field>
void printVector (Field &F, ostream &output, const vector <typename Field::element> &v) 
{
	int i;

	output << '(';
	for (i = 0; i < v.size (); i++) {
		F.write (output, v[i]);
		if (i < v.size () - 1)
			output << ", ";
	}
	output << ')' << endl;
}

template <class Field, class Polynomial>
void printPolynomial (Field &F, ostream &output, const Polynomial &v) 
{
	int i;
	int xpwr;

	for (i = 0, xpwr = 0; i < v.size () && F.isZero (v[i]); i++, xpwr++);

	for (i = v.size () - 1; i >= 0; i--) {
		if (F.isZero (v[i]))
			continue;

		if (!F.isOne (v[i]))
			F.write (output, v[i]);

		if (i > 0)
			output << " x^" << i;

		if (i > xpwr)
			output << " + ";
	}

	output << endl;
}

template <class Field, class Polynomial>
void applyPoly (const Field                                                         &F,
		vector <typename Field::element>                                    &w,
		const LinBox::Blackbox_archetype<vector <typename Field::element> > &A,
		const Polynomial                                                    &phi,
		const vector <typename Field::element>                              &v) 
{
	vector <typename Field::element> z(v.size ());
	int i, j;

	w.resize (v.size ());

	for (i = 0; i < v.size (); i++)
		F.mul (w[i], v[i], phi[phi.size () - 1]);

	for (i = phi.size () - 2; i >= 0; i--) {
		A.apply (z, w);

		for (j = 0; j < v.size (); j++)
			F.axpy (w[j], v[j], phi[i], z[j]);
	}
}

void parseArguments (int argc, char **argv, ofstream &report, Argument *args);

#endif // __TEST_COMMON_H
