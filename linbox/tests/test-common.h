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
#include "linbox/field/vector-domain.h"
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

template <class Field, class Vector>
void printVector (Field &F, ostream &output, const Vector &v) 
{ printVectorSpecialized(F, output, v, LinBox::VectorTraits<Vector>::VectorCategory()); }

template <class Field, class Vector>
void printVectorSpecialized(
		Field &F, 
		ostream &output, 
		const Vector &v, 
		LinBox::VectorCategories::DenseVectorTag tag
		)
{
	unsigned int i;

	output << '(';
	for (i = 0; i < v.size (); i++) {
		F.write (output, v[i]);
		if (i < v.size () - 1)
			output << ", ";
	}
	output << ')' << endl;
}

template <class Field, class Vector>
void printVectorSpecialized(
		Field &F, 
		ostream &output, 
		const Vector &v, 
		LinBox::VectorCategories::SparseSequenceVectorTag tag
		)
{
	typename Vector::const_iterator i;
	unsigned int j;

	output << '(';
	for (i = v.begin (), j = 0; i != v.end (); j++) {
		while (j < (*i).first) {
			output << "0, ";
			j++;
		}

		F.write (output, (*i).second);

		if (++i != v.end ())
			output << ", ";
	}
	output << ')' << endl;
}

template <class Field, class Vector>
void printVectorSpecialized(
		Field &F, 
		ostream &output, 
		const Vector &v, 
		LinBox::VectorCategories::SparseAssociativeVectorTag tag
		)
{
	typename Vector::const_iterator i;
	unsigned int j;

	output << '(';
	for (i = v.begin (), j = 0; i != v.end (); j++) {
		while (j < (*i).first) {
			output << "0, ";
			j++;
		}

		F.write (output, (*i).second);

		if (++i != v.end ())
			output << ", ";
	}
	output << ')' << endl;
}

template <class Field, class Vector>
bool areVectorsEqual (Field &F, const Vector &v, const Vector &w) 
{ return areVectorsEqualSpecialized(F, v, w, LinBox::VectorTraits<Vector>::VectorCategory()); }

template <class Field, class Vector>
bool areVectorsEqualSpecialized(
		Field &F, 
		const Vector &v, 
		const Vector &w, 
		LinBox::VectorCategories::DenseVectorTag tag
		)
{
	if (v.size() != w.size()) return false;

	for (size_t i = 0; i < v.size(); i++)
		if (!F.areEqual (w[i], v[i]))
			return false;
	
	return true;
}

template <class Field, class Vector>
bool areVectorsEqualSpecialized(
		Field &F, 
		const Vector &v, 
		const Vector &w, 
		LinBox::VectorCategories::SparseSequenceVectorTag tag
		)
{
	if (v.size() != w.size()) return false;

	typename Vector::const_iterator v_iter, w_iter;
	w_iter = w.begin();
	
	for ( v_iter = v.begin(); v_iter != v.end(); v_iter++, w_iter++)
		if ( (w_iter->first != v_iter->first) 
				|| (!F.areEqual (w_iter->second, v_iter->second)) )
			return false;
	
	return true;
}

template <class Field, class Vector>
bool areVectorsEqualSpecialized(
		Field &F, 
		const Vector &v, 
		const Vector &w, 
		LinBox::VectorCategories::SparseAssociativeVectorTag tag
		)
{
	if (v.size() != w.size()) return false;

	typename Vector::const_iterator v_iter, w_iter;
	w_iter = w.begin();
	
	for ( v_iter = v.begin(); v_iter != v.end(); v_iter++, w_iter++)
		if ( (w_iter->first != v_iter->first) 
				|| (!F.areEqual (w_iter->second, v_iter->second)) )
			return false;
	
	return true;
}

template <class Field, class Vector>
bool allZero (Field &F, const Vector &v) 
{ return allZeroSpecialized(F, v, LinBox::VectorTraits<Vector>::VectorCategory()); }

template <class Field, class Vector>
bool allZeroSpecialized(
		Field &F, 
		const Vector &v, 
		LinBox::VectorCategories::DenseVectorTag tag
		)
{
	for (size_t i = 0; i < v.size(); i++)
		if (!F.isZero (v[i]))
			return false;

	return true;
}
	
template <class Field, class Vector>
bool allZeroSpecialized(
		Field &F, 
		const Vector &v, 
		LinBox::VectorCategories::SparseSequenceVectorTag tag
		)
{
	if (0 != v.size()) 
		return false;
	else
		return true;
}

template <class Field, class Polynomial>
void printPolynomial (Field &F, ostream &output, const Polynomial &v) 
{
	int i;
	int val;

	for (val = 0; val < v.size () && F.isZero (v[val]); val++);

	if (v.size () == 0 || val == v.size ())
		output << "0";

	for (i = v.size () - 1; i >= 0; i--) {
		if (F.isZero (v[i]))
			continue;

		if (!F.isOne (v[i]) || i == 0)
			F.write (output, v[i]);

		if (i > 0)
			output << " x^" << i;

		if (i > val)
			output << " + ";
	}

	output << endl;
}

template <class Field, class Polynomial>
vector <typename Field::element> &
applyPoly (const Field                                                         &F,
	   vector <typename Field::element>                                    &w,
	   const LinBox::BlackboxArchetype<vector <typename Field::element> > &A,
	   const Polynomial                                                    &phi,
	   const vector <typename Field::element>                              &v) 
{
	typedef vector <typename Field::element> Vector;

	LinBox::VectorDomain <Field, Vector, Vector> VD (F);
	Vector z (v.size ());
	int i;

	VD.mul (w, v, phi[phi.size () - 1]);

	for (i = phi.size () - 2; i >= 0; i--) {
		A.apply (z, w);
		VD.axpy (w, v, phi[i], z);
	}

	return w;
}

/* Evaluate polynomial at a whole vector of points */

template <class Field, class Polynomial>
vector <typename Field::element> &
multiEvalPoly (const Field                            &F,
	       vector <typename Field::element>       &w,
	       const Polynomial                       &phi,
	       const vector <typename Field::element> &v) 
{
	typedef vector <typename Field::element> Vector;

	typename Field::element tmp;
	int i, j;

	w.resize (v.size ());

	for (i = 0; i < v.size (); i++)
		w[i] = phi[phi.size () - 1];

	for (i = phi.size () - 2; i >= 0; i--) {
		for (j = 0; j < v.size (); j++) {
			F.axpy (tmp, w[j], v[j], phi[i]);
			w[j] = tmp;
		}
	}

	return w;
}

/* Interpolate polynomial evaluated at a vector of points using Lagrange
 * interpolants */

template <class Field, class Polynomial>
Polynomial &
interpolatePoly (const Field                            &F,
		 Polynomial                             &f,
		 const vector <typename Field::element> &x,
		 const vector <typename Field::element> &y) 
{
	typedef vector <typename Field::element> Vector;

	int n = x.size ();

	// NB I leave one element in g always initialized to 0 as the ficticious
	// negative-first coefficient. This streamlines some of the code.
	static const int g_FUDGE = 1;
	Vector g(n + g_FUDGE);
	F.init (g[0], 0);

	typename Field::element gk, c1, c2;

	int i, j, k, d;

	f.resize (n);

	for (i = 0; i < n; i++)
		F.init (f[i], 0);

	for (j = 0; j < n; j++) {
		F.init (g[0 + g_FUDGE], 1);

		// d is the current degree of the Lagrange interpolant. i is the
		// current index in the array of x-coordonites
		for (d = 0, i = 0; d < n - 1; d++, i++) {
			if (i == j) i++;

			// Compute coefficients of this factor.
			F.sub (c1, x[j], x[i]);
			F.invin (c1);
			F.mul (c2, c1, x[i]);
			F.negin (c2);

			// Initialize the next element of the Lagrange interpolant
			F.init (g[d + 1 + g_FUDGE], 0);

			// Multiply this factor by the existing partial product
			for (k = d + 1 + g_FUDGE; k >= g_FUDGE; k--) {
				F.mul (gk, g[k - 1], c1);
				F.axpyin (gk, g[k], c2);
				g[k] = gk;
			}
		}

		for (i = 0; i < n; i++)
			F.axpyin (f[i], y[j], g[i + g_FUDGE]);
	}

	return f;
}

void parseArguments (int argc, char **argv, Argument *args);
bool isPower        (LinBox::integer n, LinBox::integer m);

#endif // __TEST_COMMON_H
