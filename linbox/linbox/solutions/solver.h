#ifndef __LINNONSING_H
#define __LINNONSING_H

#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h" 
#include "linbox/util/debug.h"
#include "linbox/field/vector-domain.h"
#include <vector>
#include <algorithm>

namespace LinBox 
{
	template <class Field, class BlackBox,class Vector>
	  Vector& solver(const BlackBox               &A,
			 Vector                       &x,		       
			 const Vector                 &y,
			 const Field                  &F)
	{
		linbox_check((x.size()==A.coldim())&&
			     (y.size()==A.rowdim()));

		typename Field::RandIter randiter(F);
		typedef std::vector<typename Field::Element> Poly;
		Poly P;
		unsigned long            deg;

		BlackboxContainer<Field, Vector> TF (&A, F, y);
		MasseyDomain< Field, BlackboxContainer<Field, Vector> > WD (&TF);
		
		WD.minpoly (P, deg);
		
		Poly::iterator p=P.begin();
		++p;
		for(;p!=P.end();++p)
		  {
		    F.divin(*p,P[0]);
		    F.negin(*p);
		  }
		
		std::copy(P.begin()+1,P.end(),P.begin());
		P.resize(P.size()-1);
		
		if(P.empty())
			return x=y;
		
		
		VectorDomain <Field> VD (F);		
		int i;		
		VD.mul (x, y, P[P.size () - 1]);
		
		for (i = P.size () - 2; i >= 0; i--)
		  {
		      A.applyIn(x);
		      VD.axpyin(x, P[i], y);
		  }
        

		return x;
	}								       		
}


#endif 
