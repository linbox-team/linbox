/* linbox/blackbox/mapleBB.h
 *
 * Written by David Saunders <saunders@cis.udel.edu>
 * See COPYING
 */

#ifndef __COMPANION_H
#define __COMPANION_H

#include "linbox/blackbox/triplesbb.h"

namespace LinBox {

template<class Field, class Vector, class Polynomial>
class Companion: public TriplesBB<Field, Vector> {

	Companion(const Field& F, const Polynomial& P)
        : TriplesBB(F, P.size()-1, P.size()-1)
	{	n = P.size() - 1;
		Field::Element one; F.init(one, 1);
	 	for (size_t i = 1; i < n; ++i) addEntry(one, i, i-1); 
	 	for (size_t i = 0; i < n; ++i) addEntry(one, P[i], n-1); 
	}// Companion cstor
 
	/** Companion cstor from random poly.  
	*Builds n by n matrix from degree n monic poly with other coefficients random.
	*/
	/*
	Companion(const Field& F, size_t n)
	: TriplesBB(F, n, n)
	{
		random source r;
		Field::Element one; F.init(one, 1);
		for (size_t i = 1; i < n; ++i) addEntry(one, i, i-1);
		for (size_t i = 0; i < n; ++i) addEntry(*r++, 1, n-1);
	}
	*/

// companion would be faster if built direct, using one axpy per entry: y_i = x_i-1 + p_i*x_n

}; //Companion class

} //LinBox
#endif //__COMPANION_H

