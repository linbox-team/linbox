/* File:	src/library/algorithms/BlackBox/ks_rank.h
 * Author:	William J. Turner for the LinBox group
 *
 * This file contains an algorithm to compute the rank of a
 * generic square black box matrix using the Kaltofen-Saunders algorithm.
 */

#ifndef _KS_RANK_
#define _KS_RANK_

#include <vector>

#include "LinBox/blackbox_archetype.h"
#include "LinBox/diagonal.h"
#include "LinBox/cekstv_switch.h"
#include "LinBox/butterfly.h"
#include "LinBox/compose.h"
#include "LinBox/transpose.h"
#include "LinBox/wiedemann_minpoly.h"

namespace LinBox
{
	
	/** Kaltofen-Saunders black box rank algorithm for a generic field.
	 * This function is templatized by the field used.
	 * @return	integer containing rank of matrix
	 * @param	F	Field in which arithmetic is done
	 * @param	BB	Black box matrix
	 * @param	r	random field element generator
	 */
	template<class Field, class Vector> 
	integer ks_rank
	(
		Field& F,
		Blackbox_archetype<Vector>& BB,
		typename Field::randIter& r
	);
}

template<class Field, class Vector> 
LinBox::integer LinBox::ks_rank(
	Field& F,
	Blackbox_archetype<Vector>& BB,
	typename Field::randIter& r
	)
{
	// Types

	typedef typename Field::element Element;
	
	// Dimensions

	integer m(BB.rowdim());
	integer n(BB.coldim());

	if (n != m)
	{
#ifdef TRACE
		cerr << "Cannot compute rank; matrix not square." << endl;
#endif // TRACE
		return 0;
	}

	// Precondition matrix with random butterfly and diagonal matrices

	std::vector<Element> 
		d(random_vector<Field,std::vector<Element> >(F, n, r)); 

	Diagonal<Field, Vector> D(F, d); // preconditioning matrix

	// Butterfly matrices
	
	integer s = count_butterfly(n);

	std::vector<Element>
		s1(random_vector<Field,std::vector<Element> >(F, s, r)),
		s2(random_vector<Field,std::vector<Element> >(F, s, r)); 

		LinBox::cekstv_switch<Field> switch1(F, s1), switch2(F, s2);

	butterfly < Vector, LinBox::cekstv_switch<Field> > 
		B1(n, switch1), B2(n, switch2);


	// Precondition matrix
	
	Compose<Vector> B1BB (&B1, &BB);
	Transpose < Vector > B2T(&B2);
	Compose<Vector> B1BBB2T (&B1BB, &B2T);
	Compose<Vector> A (&B1BBB2T, &D);	// preconditioned matrix

        // Compute minimum polynomial of matrix
        
        std::vector<Element> minpoly(wiedemann_minpoly(F, A, r));

	integer rank;

	if (!F.isZero(minpoly[0])) 
		rank = n;
	else
		rank = integer(minpoly.size() - 2);

#ifdef TRACE
	cout << "The rank is " << rank << endl;
#endif // TRACE

	return rank;

} // ks_rank(F, BB, r)

#endif // _KS_RANK_
