#ifndef __LINBOX_MIN_POL_H__
#define __LINBOX_MIN_POL_H__

#include "LinBox/lin_rand.h"                  // Random Iterator
#include "LinBox/lin_bbit.h"       // BB iterator
#include "LinBox/lin_massey.C"                // massey reccuring sequence solver

#include "LinBox/lin_methods.h"

template <class MT = MethodTrait::Wiedemann> class minpoly {
public:
	minpoly() {}

	/// A is supposed to be symmetric
	template<class BB, class Polynomial>
	void operator() (Polynomial& P, const BB& A, const MT& M = MT() ) {
		Random generator;
        	unsigned long deg;

	        BB_Container< SPBB > TF( & A, generator);
        	MasseyDom< SzCBB >  WD(&TF, M.Early_Term_Threshold() );

		WD.pseudo_minpoly(P, deg);
	}

};


#endif
