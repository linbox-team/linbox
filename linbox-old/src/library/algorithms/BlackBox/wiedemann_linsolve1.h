/* File:	src/library/algorithms/BlackBox/wiedemann_linsolve1.h
 * Author:	William J. Turner for the LinBox group
 *
 * This file contains an algorithm to compute the determinant of a
 * generic black box matrix using the Wiedemann algorithm.
 */

#ifndef _WIEDEMANN_LINSOLVE1_
#define _WIEDEMANN_LINSOLVE1_

#include <vector>

#include "LinBox/blackbox_archetype.h"
#include "LinBox/scalarprod.h"
#include "LinBox/vaxpy.h"
#include "LinBox/wiedemann_minpoly.h"

namespace LinBox
{
	
	/** Wiedemann black box nonhomogeneous linear solver algorithm for a generic field.
	 * This function is templatized by the field used.
	 * @return	field element containing determinant of matrix
	 * @param	F	Field in which arithmetic is done
	 * @param	x	vector to contain output
	 * @param	BB	Black box matrix
	 * @param	b	right hand side vector: BB * x = b
	 * @param	r	random field element generator
	 */
	template<class Field, class Vector> 
	Vector& wiedemann_linsolve1
	(
		Field& F,
		Vector& x,
		Blackbox_archetype<Vector>& BB,
		const Vector& b,
		typename Field::randIter& r
	);
}

template<class Field, class Vector> 
Vector& LinBox::wiedemann_linsolve1(
	Field& F,
	Vector& x,
	Blackbox_archetype<Vector>& BB,
	const Vector& b,
	typename Field::randIter& r
	)
{
	// Types

	
	typedef typename Field::element Element;
	typedef typename Field::randIter RandIter;
	
	Element zero, one;
	F.init(zero, 0);
	F.init(one, 1);

        // Compute minimum polynomial of matrix
        
        std::vector<Element> minpoly(wiedemann_minpoly(F, BB, r));

	// Normalize minpoly so constant term is -1
	// Saves divide in linsolve step
	
	Element factor(minpoly[0]);
	F.negin(factor);
	F.invin(factor);

#ifdef TRACE
	cout << "*** - 1 / ";
	F.write(cout, minpoly[0]);
	cout << " = ";
	F.write(cout, factor);
	cout << endl;
#endif // TRACE
	
	for (std::vector<Element>::iterator iter = minpoly.begin();
			iter != minpoly.end(); iter++)
		F.mulin(*iter, factor);

	size_t N = minpoly.size();

#ifdef TRACE
	cout << "*** The coefficients of the normalized minimum polynomial are:" << endl;
	for (size_t i = 0; i < N; i++)
	{
		cout << "\t c[" << i << "] = ";
		F.write(cout, minpoly[i]);
		cout << endl;
	}
#endif // TRACE

	// Compute solution as linear combination of Krylov vectors
	
	Vector Ab(b);

	LinBox::scalarprod(F, x, minpoly[1], b); // x = minpoly[1] * b
	
	for (size_t k = 2; k < N; k++)
	{
		BB.applyin(Ab);
		vaxpyin(F, x, minpoly[k], Ab);
	}

	return x;

} // wiedemann_linsolve1(F, x, BB, b, r)

#endif // _WIEDEMANN_LINSOLVE1_
