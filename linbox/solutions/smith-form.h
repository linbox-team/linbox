/* Copyright (C) 2010 LinBox
 *
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_smith_form_H
#define __LINBOX_smith_form_H

#include <list>
#include <vector>
#include <iterator>
#include "linbox/util/error.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/algorithms/smith-form-adaptive.h"
#include "givaro/zring.h"

namespace LinBox
{
    template <typename PIR>
    using SmithPair=std::pair<typename PIR::Element,size_t>;
    template <typename PIR>
    using SmithList=std::list<SmithPair<PIR>>;

	template<class Ring, class Vector>
	SmithList<Ring> &
	compressedSmith(SmithList<Ring> & c,
                    const Vector& SmithDiagonal, const Ring& R,
                    const size_t m, const size_t n)
    {
        typename Ring::Element si(R.one);
        size_t num=0;
        for( auto dit : SmithDiagonal ) {
            if (R.areEqual(dit,si)) ++num;
            else {
                c.emplace_back(si,num);
                num=1;
                R.assign(si,dit);
            }
        }
        if (num>0) c.emplace_back(si,num);
        num = std::min(m,n) - SmithDiagonal.size();
        R.assign(si,R.zero);
        if (num>0) c.emplace_back(si,num);
        return c;
    }

	// Convert from vector of invariants (with repeats) to EC_LIST form.
	template<class Ring>
	SmithList<Ring> &
	compressedSmith(SmithList<Ring> & c, const BlasVector<Ring>& v)
	{
		typename Ring::Element e;
		size_t count = 0;
		const Ring& R = v.field();
		size_t n = v.size();
		if (n > 0) R.assign(e, v[0]); else return c;
		count = 1;
		for (size_t i = 1; i < v.size(); ++i)
		{	if (R.areEqual(v[i], e))
				++count;
			else
			{	c.emplace_back(e, count);
				R.assign(e, v[i]); count = 1;
			}
		}
		c.emplace_back(e, count);
		return c;
	}




	/** Compute the Smith form of A.
	 * \ingroup solutions
	 *
	 * The Smith form of a linear operator A, represented as a
	 * black box, is computed over a representation of \f$Z\f$ or \f$Z_m\f$.
	 *
	 * @param[out] S a list of invariant/repcount pairs.
	 * @param A Matrix of which to compute the Smith form
	 * @param M may be a \p Method::Auto (default), which uses the
	 algorithms/smith-form-adaptive.
	 @todo Other methods will be provided later.
	 For now see the examples/smith.C
	 for ways to call other smith form algorithms.
	 */
	/*
	BB has to be dense matrix
	PL means EC_list (list of value repcount pairs)
	VL means diag of smith form as a BlasVector.
	SNF function forms:
	template<BB> smithForm(PL, BB) -> add Auto
	template<BB> smithForm(VL, BB) -> add Auto
	template<BB,Meth> smithForm(PL, BB, Meth) -> add IntegerTag
	template<BB,Meth> smithForm(VL, BB, Meth) -> add IntegerTag
	smithForm(PL, BB, IntegerTag, Auto) -> call adaptive
	smithForm(VL, BB, IntegerTag, Auto) -> call adaptive
	*/

	template <class Blackbox, class Method>
	SmithList<typename Blackbox::Field> &
    smithForm(SmithList<typename Blackbox::Field> & S,
			  const Blackbox                     & A,
			  const Method                     & M)
	{
		smithForm(S, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
		return S;
	}
	template <class Blackbox, class Method>
	BlasVector<typename Blackbox::Field> &
	smithForm(BlasVector<typename Blackbox::Field> & V,
			  const Blackbox                     & A,
			  const Method                     & M)
	{
		smithForm(V, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
		return V;
	}

	// The smithForm with default Method
	template<class Blackbox>
	SmithList<typename Blackbox::Field> &
	smithForm(SmithList<typename Blackbox::Field> & S,
			  const Blackbox& A)
	{
		smithForm(S, A, Method::Auto());
		return S;
	}
	template<class Blackbox>
	BlasVector<typename Blackbox::Field> &
	smithForm(BlasVector<typename Blackbox::Field> & V,
			  const Blackbox& A)
	{
		smithForm(V, A, Method::Auto());
		return V;
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
			if (r > 0) S.emplace_back(r, 1);
			if (n > 0) S.emplace_back(n, 0);
		}
		else
		{
			integer x; size_t c;
			for(x = p, c = 0; divides(2, x); x /= 2, ++c);

			if (x == 1 && c <= 32) // (a low power of 2)
			{
				List L;
				LocalSmith<Local2_32> SmithForm;
				SmithForm( L, M, R );
				compressedSmith(L.begin(), L.end(), S);

			}
			//  if (a odd prime power) call local-smith
			else
			{
				IliopoulosElimination::smithIn (M);
				typedef std::list< PIR::Element > List;
				List L;
				for (size_t i = 0; i < M.rowdim(); ++i) L.push_back(M[i][i]);
				compressedSmith(L.begin(), L.end(), S);
			}
		}

		return S;
	}
#endif

	/*template<>
	std::list<std::pair<integer, size_t> > &
	smithForm(std::list<std::pair<integer, size_t> >& S,
	*/
	SmithList<Givaro::ZRing<Integer>> &
	smithForm(SmithList<Givaro::ZRing<Integer>> & S,
		  const BlasMatrix<Givaro::ZRing<Integer> > 	&A,
		  const RingCategories::IntegerTag      &tag,
		  const Method::Auto			& M)
	{
		Givaro::ZRing<Integer> Z;
		BlasVector<Givaro::ZRing<Integer> > v (Z,A.rowdim() < A.coldim() ? A.rowdim() : A.coldim());
		SmithFormAdaptive::smithForm(v, A);
		//distinct(v.begin(), v.end(), S);
		return compressedSmith(S,v);
	}
	BlasVector<typename Givaro::ZRing<Integer> > &
	smithForm(BlasVector<typename Givaro::ZRing<Integer> > & V,
		  const BlasMatrix<Givaro::ZRing<Integer> > 	&A,
		  const RingCategories::IntegerTag      &tag,
		  const Method::Auto			& M)
	{
		Givaro::ZRing<Integer> Z;
		SmithFormAdaptive::smithForm(V, A);
		return V;
	}

//#endif

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

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
