/* File:	src/library/algorithms/BlackBox/wiedemann_det.h
 * Author:	William J. Turner for the LinBox group
 *
 * This file contains an algorithm to compute the determinant of a
 * generic black box matrix using the Wiedemann algorithm.
 */

#ifndef _WIEDEMANN_DET_
#define _WIEDEMANN_DET_

#include <vector>

#include "LinBox/blackbox_archetype.h"
#include "LinBox/diagonal.h"
#include "LinBox/compose.h"
#include "LinBox/vector_traits.h"
#include "LinBox/dotprod.h"
#include "LinBox/random_vector.h"
#include "LinBox/berlekamp_massey.h"

namespace LinBox
{
	
	/** Wiedemann black box determinant algorithm for a generic field.
	 * This function is templatized by the field used.
	 * @return	field element containing determinant of matrix
	 * @param	F	Field in which arithmetic is done
	 * @param	BB	Black box matrix
	 * @param	r	random field element generator
	 */
	template<class Field, class Vector> 
	typename Field::element wiedemann_det
	(
		Field& F,
		Blackbox_archetype<Vector>& BB,
		typename Field::randIter& r
	);
}

template<class Field, class Vector> 
typename Field::element 
LinBox::wiedemann_det(
	Field& F,
	Blackbox_archetype<Vector>& BB,
	typename Field::randIter& r
	)
{
	// Types

	
	typedef typename Field::element Element;
	typedef typename Field::randIter RandIter;
	
	Element zero, one;
	F.init(zero, 0);
	F.init(one, 1);

	// Dimensions

	integer n(BB.rowdim());

	// Precondition matrix with random diagonal matrix

	std::vector<Element> 
		d(random_vector<Field,std::vector<Element> >(F, n, r)); 
	
	Element prod(one);
	for (long i = 0; i < n; i++)
		F.mulin(prod, d[i]);
	
	for (long k = 0; (k < 64) && (F.isZero(prod)); k++)
	{
		F.assign(prod, one);
		for (long i = 0; i < n; i++) F.mulin(prod, d[i]);
	}

	Diagonal<Field, Vector> D(F, d); // preconditioning matrix
	Compose< Vector > DA(&D, &BB);	// preconditioned matrix

#ifdef TRACE
	cout << "The random elements of the diagonal matrix are:" << endl;
	for (long i = 0; i < n; i++)
	{
		cout << "i = " << i << ", \t D[i,i] = ";
		F.write(cout, d[i]);
		cout << endl;
	}
		
	cout << "Their product is ";
	F.write(cout, prod);
	cout << endl;
	
#endif // TRACE

	// Random vectors u and v

	//random_vector<Field, Vector>(F, n, r);

	Vector u(random_vector<Field, Vector>(F, n, r)); 
	Vector v(random_vector<Field, Vector>(F, n, r)); 
		
	integer N = 2*n;
	std::vector<Element> a(N);

	F.init(a[0], 0);
	F.assign(a[0], dotprod(F, u, v));

	// Compute sequence u^T A^i v and minimal polynomial

	for (long i = 1; i < N; i++)
	{
		DA.applyin(v);		// (DA)^i * v
		F.init(a[i], 0);
		F.assign(a[i], dotprod(F, u, v));
	} // for (long i = 1; i < 2*n; i++)

#ifdef TRACE
	cout << "The scalar sequence is:" << endl;
	for (unsigned long i = 0; i < a.size(); i++)
	{
		cout << "i = " << i << ", \t a[i] = ";
		F.write(cout, a[i]);
		cout << endl;
	}
#endif // TRACE

	std::vector<Element> minpoly(berlekamp_massey(F, a));

	if (F.isZero(minpoly[0])) return zero;
		
	if (n != integer(minpoly.size() - 1))
	{
		cerr << "Bad projection in Wiedemann algorithm" << endl;
		return zero;
	}

	Element det(minpoly[0]);
	F.divin(det, prod);
	if (0 != (n % 2)) F.negin(det);

#ifdef TRACE
	cout << "The determinant is ";
	F.write(cout, det);
	cout << endl;
#endif // TRACE

	return det;

} // wiedemann_det(F, BB, r)

#endif // _WIEDEMANN_DET_
