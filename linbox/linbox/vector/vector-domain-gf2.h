/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/vector/vector-domain.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * There are now two vector domain types: VectorDomainBase and
 * VectorDomain. VectorDomainBase, which in principle can be used independently,
 * contains all functions that require only one vector type (such as axpy, mul,
 * read, and write). VectorDomain inherits VectorDomainBase and implements
 * dotprod, which requires two vector types.
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Added the modifications for categories and vector traits that were designed
 * at the Rootbeer meeting. Added parametrization of VectorTags by VectorTraits.
 *
 * ------------------------------------
 * 2002-06-04 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Reverted change of 2002-04-10, reintegrating VectorDomain and
 * VectorDomainBase. Now using template specialization on the functions, now
 * that I know how to do it.
 *
 * ------------------------------------
 * 2002-06-21 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added methods add, addin, sub, subin, areEqual, isZero, and copy.
 *
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_field_vector_domain_gf2_H
#define __LINBOX_field_vector_domain_gf2_H

#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"

// DotProduct
namespace LinBox
{
	// Specialization of DotProductDomain for GF2


	template<>
	class DotProductDomain<GF2> : private virtual VectorDomainBase<GF2> {
	public:

		typedef bool Element;

		DotProductDomain (const GF2 &F) :
			VectorDomainBase<GF2> (F)
		{}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline BitVector::reference dotSpecializedDD (BitVector::reference res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline BitVector::reference dotSpecializedDSP (BitVector::reference res, const Vector1 &v1, const Vector2 &v2) const;
	};


}


// Specialization of vector domain
namespace LinBox
{

	template <>
	class VectorDomain<GF2> : private virtual VectorDomainBase<GF2>, private virtual DotProductDomain<GF2> {
	public:
		typedef bool Element;

		VectorDomain (const VectorDomain &VD) :
			VectorDomainBase<GF2> (VD._field), DotProductDomain<GF2> (VD._field)
		{}

		VectorDomain &operator = (const VectorDomain &) { return *this; }

		const GF2 &field () const
		{
			return _field;
		}

		template <class Vector>
		inline std::ostream &write (std::ostream &os, const Vector &x) const
		{
			return writeSpecialized (os, x, typename VectorTraits<Vector>::VectorCategory ());
		}

		template <class Vector>
		inline std::istream &read (std::istream &is, Vector &x) const
		{
			return readSpecialized (is, x, typename VectorTraits<Vector>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline Vector1 &copy (Vector1 &res, const Vector2 &v) const
		{
			return copySpecialized (res, v,
					  typename VectorTraits<Vector1>::VectorCategory (),
					  typename VectorTraits<Vector2>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline Vector1 &copy (Vector1 &res, const Vector2 &v, size_t i, size_t len = 0) const
		{
			return copySpecialized (res, v, i, len,
					  typename VectorTraits<Vector1>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline bool areEqual (const Vector1 &v1, const Vector2 &v2) const
		{
			return areEqualSpecialized (v1, v2,
					      typename VectorTraits<Vector1>::VectorCategory (),
					      typename VectorTraits<Vector2>::VectorCategory ());
		}

		template <class Vector>
		inline bool isZero (const Vector &v) const
		{
			return isZeroSpecialized (v, typename VectorTraits<Vector>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline Element &dot (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{
			return dotSpecialized (res, v1, v2,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline BitVector::reference dot (BitVector::reference res, const Vector1 &v1, const Vector2 &v2) const
		{
			return dotSpecialized (res, v1, v2,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{
			return dot (res, v1, v2);
		}

		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &add (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{
			return addSpecialized (res, y, x,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory (),
					 typename VectorTraits<Vector3>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline Vector1 &addin (Vector1 &y, const Vector2 &x) const
		{
			return addinSpecialized (y, x,
					   typename VectorTraits<Vector1>::VectorCategory (),
					   typename VectorTraits<Vector2>::VectorCategory ());
		}

		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &sub (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{
			return addSpecialized (res, y, x,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory (),
					 typename VectorTraits<Vector3>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline Vector1 &subin (Vector1 &y, const Vector2 &x) const
		{
			return addinSpecialized (y, x,
					   typename VectorTraits<Vector1>::VectorCategory (),
					   typename VectorTraits<Vector2>::VectorCategory ());
		}

		template <class Vector1, class Vector2>
		inline Vector1 &neg (Vector1 &res, const Vector2 &x) const
		{
			copy (res, x);
			return res;
		}

		template <class Vector>
		inline Vector &negin (Vector &y) const
		{
			return y;
		}

		template <class Vector1, class Vector2>
		inline Vector1 &mul (Vector1 &res, const Vector2 &x, const Element a) const
		{
			return mulSpecialized (res, x, a, typename VectorTraits<Vector1>::VectorCategory ());
		}

		template <class Vector>
		inline Vector &mulin (Vector &x, const Element a) const
		{
			return mulinSpecialized (x, a, typename VectorTraits<Vector>::VectorCategory ());
		}

		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &axpy (Vector1 &res, const Element a, const Vector2 &x, const Vector3 &y) const
		{
			if (a)
				add (res, x, y);
			else
				this->copy (res, y);
			return res;
		}

		template <class Vector1, class Vector2>
		inline Vector1 &axpyin (Vector1 &y, const Element a, const Vector2 &x) const
		{
			if (a)
				addin (y, x);
			return y;
		}

		VectorDomain (const GF2 &F) :
			VectorDomainBase<GF2> (F), DotProductDomain<GF2> (F)
		{}


		// Specialized function implementations
		template <class Vector>
		std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
						VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector>
		std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
						VectorCategories::SparseZeroOneVectorTag) const;

		template <class Vector>
		std::istream &readSpecialized (std::istream &is, const Vector &x,
					       VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector>
		std::istream &readSpecialized (std::istream &is, const Vector &x,
					       VectorCategories::SparseZeroOneVectorTag) const;

		template <class Vector1, class Vector2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::DenseZeroOneVectorTag,
					  VectorCategories::DenseZeroOneVectorTag) const;

		template <class Vector1, class Vector2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::DenseZeroOneVectorTag,
					  VectorCategories::SparseZeroOneVectorTag) const;
		template <class Vector1, class Vector2>
		inline bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						 VectorCategories::SparseZeroOneVectorTag,
						 VectorCategories::DenseZeroOneVectorTag) const
		{
			return areEqual (v2, v1);
		}

		template <class Vector1, class Vector2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseZeroOneVectorTag,
					  VectorCategories::SparseZeroOneVectorTag) const;


		template <class Vector>
		bool isZeroSpecialized (const Vector &v, VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector>
		inline bool isZeroSpecialized (const Vector &v,
					       VectorCategories::SparseZeroOneVectorTag) const
		{
			return v.empty ();
		}

		template <class Vector1, class Vector2>
		inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
						 VectorCategories::DenseZeroOneVectorTag,
						 VectorCategories::DenseZeroOneVectorTag) const
		{
			std::copy (v.wordBegin (), v.wordEnd (), res.wordBegin ()); return res;
		}

		template <class Vector1, class Vector2>
		inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len, VectorCategories::DenseZeroOneVectorTag) const
		{
			std::copy (v.begin (), (len == 0) ? v.end () : v.begin () + len, res.begin () + i);
			return res;
		}

		template <class Vector1, class Vector2>
		inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len, VectorCategories::DenseVectorTag) const
		{
			std::copy (v.begin (), (len == 0) ? v.end () : v.begin () + len, res.begin () + i);
			return res;
		}



		template <class Vector1, class Vector2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseZeroOneVectorTag,
					  VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector1, class Vector2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::DenseZeroOneVectorTag,
					  VectorCategories::SparseZeroOneVectorTag) const;
		template <class Vector1, class Vector2>
		inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
						 VectorCategories::SparseZeroOneVectorTag,
						 VectorCategories::SparseZeroOneVectorTag) const
		{ res = v; return res; }

		template <class Vector1, class Vector2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::DenseZeroOneVectorTag,
						VectorCategories::DenseZeroOneVectorTag) const
		{
			return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2);
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::DenseZeroOneVectorTag,
						VectorCategories::SparseZeroOneVectorTag) const
		{
			return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2);
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseZeroOneVectorTag,
						VectorCategories::DenseZeroOneVectorTag) const
		{
			return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1);
		}

		template <class Vector1, class Vector2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag) const;

		template <class Vector1, class Vector2>
		inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
							    VectorCategories::DenseZeroOneVectorTag,
							    VectorCategories::DenseZeroOneVectorTag) const
		{
			return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2);
		}

		template <class Vector1, class Vector2>
		inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
							    VectorCategories::DenseZeroOneVectorTag,
							    VectorCategories::SparseZeroOneVectorTag) const
		{
			return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2);
		}

		template <class Vector1, class Vector2>
		inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
							    VectorCategories::SparseZeroOneVectorTag,
							    VectorCategories::DenseZeroOneVectorTag) const
		{
			return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1);
		}

		template <class Vector1, class Vector2>
		BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
						     VectorCategories::SparseZeroOneVectorTag,
						     VectorCategories::SparseZeroOneVectorTag) const;

		template <class Vector1, class Vector2, class Vector3>
		Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::DenseZeroOneVectorTag,
					 VectorCategories::DenseZeroOneVectorTag,
					 VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector1, class Vector2, class Vector3>
		Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::DenseZeroOneVectorTag,
					 VectorCategories::DenseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag) const
		{
			copy (res, y);
			addin (res, x);
		}

		template <class Vector1, class Vector2, class Vector3>
		Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag) const;

		template <class Vector1, class Vector2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::DenseZeroOneVectorTag,
					   VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector1, class Vector2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::DenseZeroOneVectorTag,
					   VectorCategories::SparseZeroOneVectorTag) const;
		template <class Vector1, class Vector2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseZeroOneVectorTag,
					   VectorCategories::DenseZeroOneVectorTag) const
		{
			Vector1 xp, res;
			copy (xp, x);
			add (res, y, xp);
			copy (y, res);
			return y;
		}

		template <class Vector1, class Vector2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseZeroOneVectorTag,
					   VectorCategories::SparseZeroOneVectorTag) const
		{
			Vector1 res;
			add (res, y, x);
			this->copy (y, res);
			return y;
		}

		template <class Vector1, class Vector2>
		Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
					 VectorCategories::DenseZeroOneVectorTag ) const
		{
			if (a)
				this->copy (res, x);
			else
				std::fill (res.wordBegin (), res.wordEnd (), 0);
			return res;
		}

		template <class Vector1, class Vector2>
		Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
					 VectorCategories::SparseZeroOneVectorTag ) const
		{
			if (a)
				this->copy (res, x);
			else
				res.clear ();
			return res;
		}

		template <class Vector>
		inline Vector &mulinSpecialized (Vector &x, const Element a,
						 VectorCategories::DenseZeroOneVectorTag) const
		{
			if (!a)
				std::fill (x.wordBegin (), x.wordEnd (), 0);
			return x;
		}

		template <class Vector>
		inline Vector &mulinSpecialized (Vector &x, const Element a,
						 VectorCategories::SparseZeroOneVectorTag ) const
		{
			if (!a)
				x.clear ();
			return x;
		}

		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::GenericVectorTag,
						VectorCategories::GenericVectorTag,
						VectorCategories::GenericVectorTag) const
		{
			typename LinBox::Vector<GF2>::Sparse v;
			typename LinBox::Vector<GF2>::Sparse w;
			typename LinBox::Vector<GF2>::Sparse u;

			copy (v, x);
			copy (w, y);
			add (u, w, v);
			copy (res, u);

			return u;
		}

		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::GenericVectorTag,
						VectorCategories::GenericVectorTag,
						VectorCategories::GenericVectorTag) const
		{
			typename LinBox::Vector<GF2>::Sparse v;
			typename LinBox::Vector<GF2>::Sparse w;
			typename LinBox::Vector<GF2>::Sparse u;

			copy (v, x);
			copy (w, y);
			sub (u, w, v);
			copy (res, u);

			return u;
		}
	};

}



#include "linbox/vector/vector-domain-gf2.inl"

#endif // __LINBOX_field_vector_domain_H
