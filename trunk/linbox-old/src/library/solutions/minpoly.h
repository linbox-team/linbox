#ifndef __LINBOX_MINPOLY_H__
#define __LINBOX_MINPOLY_H__

#include "LinBox/lin_rand.h"       // Random Iterator
#include "LinBox/lin_bbit.h"       // BB iterator
#include "LinBox/lin_massey.C"     // massey recurring sequence solver

#include "LinBox/lin_methods.h"

template < class Polynomial, class BB, class Vector, class MT = MethodTrait::Wiedemann >
class Minpoly {
public:
	Polynomial & operator() (Polynomial & P,
				 const BB & A,
				 const MT & M = MT()) {

		Random          generator;
		unsigned long   deg;

		BB_Container < BB, Vector > TF(&A, generator);
		MasseyDom < BB_Container < BB > >WD(&TF, M.Early_Term_Threshold());

		WD.pseudo_minpoly(P, deg);
	};

};

#endif
