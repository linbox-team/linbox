/* File:	src/library/algorithms/BlackBox/wiedemann_minpoly.h
 * Author:	William J. Turner for the LinBox group
 *
 * This file contains an algorithm to compute the minimum polynomial of a 
 * generic (square) black box matrix using the Wiedemann algorithm.
 */

#ifndef _WIEDEMANN_MINPOLY_
#define _WIEDEMANN_MINPOLY_

#include <vector>

#include "LinBox/blackbox_archetype.h"
#include "LinBox/dotprod.h"
#include "LinBox/random_vector.h"
#include "LinBox/berlekamp_massey.h"

namespace LinBox
{
	
	/** Wiedemann black box determinant algorithm for a generic field.
	 * This function is templatized by the field used.
	 * @return	STL vector containing coefficients of minimum polynomial
	 * @param	F	Field in which arithmetic is done
	 * @param	BB	Black box matrix
	 * @param	r	random field element generator
	 */
	template<class Field, class Vector> 
	std::vector<typename Field::element>
	wiedemann_minpoly
	(
		Field& F,
		Blackbox_archetype<Vector>& BB,
		typename Field::randIter& r
	);
}

template<class Field, class Vector> 
std::vector<typename Field::element>
LinBox::wiedemann_minpoly(
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
		BB.applyin(v);		// (BB)^i * v
		F.init(a[i], 0);
		F.assign(a[i], dotprod(F, u, v));
	}

#ifdef TRACE
	cout << "The scalar sequence is:" << endl;
	for (size_t long i = 0; i < a.size(); i++)
	{
		cout << "\t a[" << i << "] = ";
		F.write(cout, a[i]);
		cout << endl;
	}
#endif // TRACE

	std::vector<Element> minpoly(berlekamp_massey(F, a));

#ifdef TRACE
	cout << "*** The degree of the minimal polynomial is " 
		<< (minpoly.size() - 1) << endl;

	cout << "*** The coefficients of the minimum polynomial are:" << endl;
	for (size_t i = 0; i < minpoly.size(); i++)
	{
		cout << "\tc[" << i << "] = ";
		F.write(cout, minpoly[i]);
		cout << endl;
	}
#endif // TRACE

	return minpoly;

} // wiedemann_minpoly(F, BB, r)

#endif // _WIEDEMANN_MINPOLY_
