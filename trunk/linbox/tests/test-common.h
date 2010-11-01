/* linbox/tests/test-common.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Added parametrization to the VectorCategory tags to make them fit the
 * Rootbeer meeting design of VectorCategories being parametrized by
 * VectorTraits.
 * 
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __LINBOX_test_common_H
#define __LINBOX_test_common_H

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/archetype.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/integer.h"

using namespace std;

enum ArgumentType {
	TYPE_NONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE
};
#define TYPE_BOOL TYPE_NONE

struct Argument 
{
	char             c;
	const char            *example;
	const char            *helpString;
	ArgumentType     type;
	void            *data;
};
// example may be passed as null and will be generated intelligently
// eg "-b {YN+-}" for bools, "-v v" for all else

template <class Field, class Vector>
void printVector (Field &F, ostream &output, const Vector &v) 
{ printVectorSpecialized(F, output, v, typename LinBox::VectorTraits<Vector>::VectorCategory()); }

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
	size_t val;

	for (val = 0; val < v.size () && F.isZero (v[val]); val++) ;

	if (v.size () == 0 || val == v.size ())
		output << "0";

	for (i = v.size () - 1; i >= 0; i--) {
		if (F.isZero (v[i]))
			continue;

		if (!F.isOne (v[i]) || i == 0)
			F.write (output, v[i]);

		if (i > 0)
			output << " x^" << i;

		if (i > (int) val)
			output << " + ";
	}

	output << endl;
}

template <class Field, class Blackbox, class Polynomial, class Vector>
vector <typename Field::Element> &
applyPoly (const Field                             &F,
	   Vector                                  &w,
	   const Blackbox			   &A,
	   const Polynomial                        &phi,
	   const Vector                            &v) 
{
	LinBox::VectorDomain <Field> VD (F);
	Vector z;
	int i;

	LinBox::VectorWrapper::ensureDim (z, A.rowdim ());

	VD.mul (w, v, phi[phi.size () - 1]);

	for (i = phi.size () - 2; i >= 0; i--) {
		A.apply (z, w);
		VD.axpy (w, phi[i], v, z);
	}

	return w;
}

/* Evaluate polynomial at a whole vector of points */

template <class Field, class Polynomial>
vector <typename Field::Element> &
multiEvalPoly (const Field                            &F,
	       vector <typename Field::Element>       &w,
	       const Polynomial                       &phi,
	       const vector <typename Field::Element> &v) 
{
	typedef vector <typename Field::Element> Vector;

	typename Field::Element tmp;
	int i;
	size_t j;

	w.resize (v.size ());

	for (j = 0; j < v.size (); j++)
		w[j] = phi[phi.size () - 1];

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
		 const vector <typename Field::Element> &x,
		 const vector <typename Field::Element> &y) 
{
	typedef vector <typename Field::Element> Vector;

	int n = x.size ();

	// NB I leave one element in g always initialized to 0 as the ficticious
	// negative-first coefficient. This streamlines some of the code.
	static const int g_FUDGE = 1; 
	Vector g(n + g_FUDGE);
	F.init (g[0], 0);

	typename Field::Element gk, c1, c2;

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

void parseArguments (int argc, char **argv, Argument *args, bool printDefaults = true);

/** writes the values of all arguments, preceded by the programName */
std::ostream& writeCommandString (std::ostream& os, Argument *args, char* programName);

bool isPower        (LinBox::integer n, LinBox::integer m);

/* Give an approximation of the value of the incomplete gamma function at a, x,
 * to within the tolerance tol */

extern inline double incompleteGamma (double a, double x, double tol);

/* Give the value of the chi-squared cumulative density function for given
 * value of chi_sqr and the given degrees of freedom */

double chiSquaredCDF (double chi_sqr, double df);

#ifdef LinBoxSrcOnly
#include "test-common.C"
#endif
#endif // __LINBOX_test_common_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
