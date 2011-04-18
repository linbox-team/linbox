/* Copyright (C) 2010 LinBox
 * 
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_smith_form_H
#define __LINBOX_smith_form_H

#include <list>
#include <vector>
#include <linbox/util/error.h>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/algorithms/smith-form-adaptive.h>
//#include <linbox/algorithms/smith-form.h>
//#include <linbox/algorithms/smith-form-local.h>

namespace LinBox
{
	template<class I1, class Lp>
	void distinct (I1 a, I1 b, Lp& c)
	{	typename I1::value_type e;
  		size_t count = 0;
  		if (a != b) {e = *a; ++a; count = 1;}
  		else return;
  		while (a != b)
  		{	if (*a == e) 
				++count;
     			else
     			{	c.push_back(typename Lp::value_type(e, count));
       				e = *a; count = 1;
     			}
     			++a;
  		}
  		c.push_back(typename Lp::value_type(e, count));
  		return;
	}


	/** Compute the Smith form of A
	 *
	 * The Smith form of a linear operator A, represented as a
	 * black box, is computed over a representation of Z or Z_m.
	 *
	 * @param Output S, a list of invariant/repcount pairs.
	 * @param A Matrix of which to compute the Smith form
	 * @param M may be a Method::Hybrid (default), which uses the 
	 algorithms/smith-form-adaptive.  
	 Other methods will be provided later.  For now see the examples/smith.C
	 for ways to call other smith form algorithms.
         \ingroup solutions
        */
    template <class Output, class Blackbox, class MyMethod>
    Output &smithForm(Output & S, 
		const Blackbox                              &A,
		const MyMethod                           &M) 
    {
	    smithForm(S, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	    return S;
    }

        // for specialization with respect to the DomainCategory
    template< class Output, class Blackbox, class SmithMethod, class DomainCategory>
    Output &smithForm(Output & S, 
                      const Blackbox        &A,
                      const DomainCategory  &tag,
                      const SmithMethod  &M)      
    {
        throw LinboxError( "Smith form solution implemented only for DenseMatrix<NTL_ZZ>");
    }

	// The smithForm with default Method 
    template<class Output, class Blackbox>
    Output &smithForm(Output& S, 
					  const Blackbox& A) {

        smithForm(S, A, Method::Hybrid());
		return S;
    }

#if 0
	// The smithForm for ModularTag 
    template<class Output, class Blackbox, class MyMethod>
    Output &smithForm(Output & S, 
        const Blackbox                            &A,
        const RingCategories::ModularTag          &tag,
		const MyMethod& M)
    {
		typename Blackbox::Field F = A.field();
		integer p, c; F.characteristic(p); F.cardinality(c);
		if (probab_prime(p) && p == c) 
		{	size_t r; rank(r, A);
			S.resize(0);
			size_t n = (A.rowdim() > A.coldim() ? A.coldim() : A.rowdim())-r;
			if (r > 0) S.push_back( std::pair<size_t, integer>(r, 1) );
			if (n > 0) S.push_back( std::pair<size_t, integer>(n, 0) );
		}
		else 
		{ 
			integr x; size_t c; 
			for(x = p, c = 0; divides(2, x); x /= 2, ++c);
		 
		 	if (x == 1 && c <= 32) // (a low power of 2) 
		  	{
		    	List L;
				LocalSmith<Local2_32> SmithForm;
				SmithForm( L, M, R );
				distinct(L.begin(), L.end(), S);

		  	}
		//  if (a odd prime power) call local-smith
		  	else 
			{
				IliopoulosElimination::smithIn (M);
				typedef std::list< PIR::Element > List;
				List L;
				for (size_t i = 0; i < M.rowdim(); ++i) L.push_back(M[i][i]);
				distinct(L.begin(), L.end(), S);
		  	}
		}
		  
		  return S;
    }
#endif

	// The smithForm with Hybrid Method 
    template<>
    std::list<std::pair<integer, size_t> > &smithForm(std::list<std::pair<integer, size_t> >& S, 
        const DenseMatrix<NTL_ZZ> 	&A,
        const RingCategories::IntegerTag          &tag,
		const Method::Hybrid& M)
    {
		std::vector<integer> v (A.rowdim() < A.coldim() ? A.rowdim() : A.coldim());
		SmithFormAdaptive::smithForm(v, A);
		distinct(v.begin(), v.end(), S);

		return S;
    }

#if 0
	// The smithForm with Elimination Method 
    template<class Output, class Ring>
    Output &smithForm(Output & S, 
		const DenseMatrix<Ring> &A,
		const RingCategories::IntegerTag          &tag,
		const Method::Elimination& M)
    {
		typename Ring::Element d;
		det(d, A, tag, M); // or just use default hybrid?  What does elim mean?
		integer D;
		A.field().convert(D, d);
		if (D < Modular<int>::MaxModulus)
		{  typedef Modular<int> Ring2;
			Ring2 R2(D);
			MatrixHom::map(B, A, R2);
	    	IliolopousElimination::smithIn(B);
        //return diagonal of B in Output object.
		}
		else
		{  typedef Modular<integer> Ring2;
			Ring2 R2(D);
			MatrixHom::map(B, A, R2);
	    	IliolopousElimination::smithIn(B);
        //return diagonal of B in Output object.
		}

    }
#endif

#if 0
	// The smithForm with BlackBox Method 
    template<class Output, class Blackbox>
    Output &smithForm(Output & S, 
		const Blackbox                      &A,
		const RingCategories::IntegerTag    &tag,
		const Method::Blackbox              &M)
    {
		// this will be binary search smith form (EGV')
    }
#endif
    

} // end of LinBox namespace
#endif // __LINBOX_smith_form_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
