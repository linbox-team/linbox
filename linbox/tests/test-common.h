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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_test_common_H
#define __LINBOX_test_common_H

#include <iostream>
#include <fstream>
// #include <vector>
#include "linbox/vector/blas-vector.h"

#include "linbox/linbox-config.h"
#include "linbox/field/archetype.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/integer.h"

using namespace std;
#include "linbox/util/commentator.h"
#include "fflas-ffpack/utils/args-parser.h"


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
		F.write (output, v[(size_t)i]);
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
void printVector (Field &F, ostream &output, const Vector &v)
{
	printVectorSpecialized(F, output, v, typename LinBox::VectorTraits<Vector>::VectorCategory());
}

template <class Field, class Vector>
bool areVectorsEqual (Field &F, const Vector &v, const Vector &w)
{
	return areVectorsEqualSpecialized(F, v, w, LinBox::VectorTraits<Vector>::VectorCategory());
}

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
		if (!F.areEqual (w[(size_t)i], v[(size_t)i]))
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

	for ( v_iter = v.begin(); v_iter != v.end(); ++v_iter, ++w_iter)
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

	for ( v_iter = v.begin(); v_iter != v.end(); ++v_iter, ++w_iter)
		if ( (w_iter->first != v_iter->first)
		     || (!F.areEqual (w_iter->second, v_iter->second)) )
			return false;

	return true;
}

template <class Field, class Vector>
bool allZero (Field &F, const Vector &v)
{
       	return allZeroSpecialized(F, v, LinBox::VectorTraits<Vector>::VectorCategory());
}

template <class Field, class Vector>
bool allZeroSpecialized(
			Field &F,
			const Vector &v,
			LinBox::VectorCategories::DenseVectorTag tag
		       )
{
	for (size_t i = 0; i < v.size(); i++)
		if (!F.isZero (v[(size_t)i]))
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

	for (i = (int)v.size () - 1; i >= 0; i--) {
		if (F.isZero (v[(size_t)i]))
			continue;

		if (!F.isOne (v[(size_t)i]) || i == 0)
			F.write (output, v[(size_t)i]);

		if (i > 0)
			output << " x^" << i;

		if (i > (int) val)
			output << " + ";
	}

	output << endl;
}

template <class Field, class Blackbox, class Polynomial, class Vector>
LinBox::BlasVector <Field> &
applyPoly (const Field                             &F,
	   Vector                                  &w,
	   const Blackbox			   &A,
	   const Polynomial                        &phi,
	   const Vector                            &v)
{
	LinBox::VectorDomain <Field> VD (F);
	Vector z(F);
	int i;

	LinBox::VectorWrapper::ensureDim (z, A.rowdim ());

	VD.mul (w, v, phi[phi.size () - 1]);

	for (i = (int)phi.size () - 2; i >= 0; i--) {
		A.apply (z, w);
		VD.axpy (w, phi[(size_t)i], v, z);
	}

	return w;
}

/* Evaluate polynomial at a whole vector of points */

template <class Field, class Polynomial, class Vector>
Vector &
multiEvalPoly (const Field        &F,
	       Vector             &w,
	       const Polynomial   &phi,
	       const Vector       &v)
{

	typename Field::Element tmp;
	int i;
	size_t j;

	w.resize (v.size ());

	for (j = 0; j < v.size (); j++)
		w[(size_t)j] = phi[phi.size () - 1];

	for (i = (int)phi.size () - 2; i >= 0; i--) {
		for (j = 0; j < v.size (); j++) {
			F.axpy (tmp, w[(size_t)j], v[(size_t)j], phi[(size_t)i]);
			w[(size_t)j] = tmp;
		}
	}

	return w;
}

/* Interpolate polynomial evaluated at a vector of points using Lagrange
 * interpolants */

template <class Field, class Polynomial>
Polynomial &
interpolatePoly (const Field             &F,
		 Polynomial              &f,
		 const LinBox::BlasVector<Field> &x,
		 const LinBox::BlasVector<Field> &y)
{
	typedef LinBox::BlasVector<Field> Vector;

	int n = (int)x.size ();

	// NB I leave one element in g always initialized to 0 as the ficticious
	// negative-first coefficient. This streamlines some of the code.
	static const int g_FUDGE = 1;
	Vector g(F,(size_t)(n + g_FUDGE));
	F.init (g[0], 0);

	typename Field::Element gk, c1, c2;

	int i, j, k, d;

	f.resize ((size_t)n);

	for (i = 0; i < n; i++)
		F.init (f[(size_t)i], 0);

	for (j = 0; j < n; j++) {
		F.init (g[0 + g_FUDGE], 1);

		// d is the current degree of the Lagrange interpolant. i is the
		// current index in the array of x-coordonites
		for (d = 0, i = 0; d < n - 1; d++, i++) {
			if (i == j) i++;

			// Compute coefficients of this factor.
			F.sub (c1, x[(size_t)j], x[(size_t)i]);
			F.invin (c1);
			F.mul (c2, c1, x[(size_t)i]);
			F.negin (c2);

			// Initialize the next element of the Lagrange interpolant
			F.init (g[(size_t)(d + 1 + g_FUDGE)], 0);

			// Multiply this factor by the existing partial product
			for (k = d + 1 + g_FUDGE; k >= g_FUDGE; k--) {
				F.mul (gk, g[(size_t)k - 1], c1);
				F.axpyin (gk, g[(size_t)k], c2);
				g[(size_t)k] = gk;
			}
		}

		for (i = 0; i < n; i++)
			F.axpyin (f[(size_t)i], y[(size_t)j], g[(size_t)i + g_FUDGE]);
	}

	return f;
}


bool isPower        (LinBox::integer n, LinBox::integer m);

/* Give an approximation of the value of the incomplete gamma function at a, x,
 * to within the tolerance tol */

extern inline double incompleteGamma (double a, double x, double tol);

/* Give the value of the chi-squared cumulative density function for given
 * value of chi_sqr and the given degrees of freedom */

double chiSquaredCDF (double chi_sqr, double df);

#ifdef LinBoxTestOnly
#include "test-common.C"
#endif
#endif // __LINBOX_test_common_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
