/* Copyright (C) 2010 LinBox
 * Adapted by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * (from other modular-balanced* files)
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

/*! @file field/modular/modular-int64.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c int64_t .
 */
#ifndef __LINBOX_modular_int64_H
#define __LINBOX_modular_int64_H


#include <cmath>
#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/ring/modular.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"

#include <givaro/modular-integral.h>

#ifndef LINBOX_MAX_INT64 /*  18446744073709551615L(L) is UINT64_MAX*/
#ifdef __x86_64__
#define LINBOX_MAX_INT64 INT64_MAX
#else
#define LINBOX_MAX_INT64 INT64_MAX
#endif
#endif

// Namespace in which all LinBox code resides
namespace LinBox
{

	template<class Field>
	class DotProductDomain;
	template<class Field>
	class FieldAXPY;
	template<class Field>
	class MVProductDomain;

	template <typename Compute_t>
	class FieldAXPY<Givaro::Modular<int64_t,Compute_t> > {
	public:

		typedef int64_t Element;
		typedef Givaro::Modular<int64_t,Compute_t> Field;

		FieldAXPY (const Field &F) : _field (&F), _y(0)
		{
			_two_64 = (uint64_t(1) << 32) % uint64_t(F.characteristic());
			_two_64 = (_two_64 * _two_64) % uint64_t(F.characteristic());
		}

		FieldAXPY (const FieldAXPY &faxpy) :
			_two_64 (faxpy._two_64), _field (faxpy._field), _y (0)
		{}

		FieldAXPY<Field> &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			_two_64 = faxpy._two_64;
			return *this;
		}

		inline const Field & field() const { return *_field; }

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) a * (uint64_t) x;
			_y += t;
			if (_y < t)
				return _y += _two_64;
			else
				return _y;
		}

		inline uint64_t& accumulate (const Element &t)
		{
			_y += (uint64_t)t;
			if (_y < (uint64_t)t)
				return _y += _two_64;
			else
				return _y;
		}

		inline Element& get (Element &y)
		{
			y =Element(_y % (uint64_t) field().characteristic());
			return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = (uint64_t)y;
			return *this;
		}

		inline void reset()
		{
			_y = 0;
		}

	public:
		uint64_t _two_64;

	protected:
		const Field *_field;
		uint64_t _y;
	};


	template <typename Compute_t>
	class DotProductDomain<Givaro::Modular<int64_t,Compute_t> > : public VectorDomainBase<Givaro::Modular<int64_t,Compute_t> > {

	public:
		typedef int64_t Element;
		typedef Givaro::Modular<int64_t,Compute_t> Field;
		using VectorDomainBase<Field>::faxpy;
		using VectorDomainBase<Field>::field;
		DotProductDomain(){}
		DotProductDomain (const Field &F) :
			VectorDomainBase<Field> (F)
		{}


	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;

			uint64_t y = 0;
			uint64_t t;

			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j)
			{
				t = ( (uint64_t) *i ) * ( (uint64_t) *j );
				y += t;

				if (y < t)
					y += faxpy()._two_64;
			}

			y %= (uint64_t) field().characteristic();
			return res = (Element)y;

		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;

			uint64_t y = 0;
			uint64_t t;

			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt)
			{
				t = ( (uint64_t) *i_elt ) * ( (uint64_t) v2[*i_idx] );
				y += t;

				if (y < t)
					y += faxpy()._two_64;
			}


			y %= (uint64_t) field().characteristic();

			return res = (Element) y;
		}
	};

	// Specialization of MVProductDomain for int64_t modular field

	template <typename Compute_t>
	class MVProductDomain<Givaro::Modular<int64_t,Compute_t> > {
	public:

		typedef int64_t Element;
		typedef Givaro::Modular<int64_t,Compute_t> Field;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized
			(VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const
                {

                        linbox_check (A.coldim () == v.size ());
                        linbox_check (A.rowdim () == w.size ());

                        typename Matrix::ConstColIterator i = A.colBegin ();
                        typename Vector2::const_iterator j;
                        typename Matrix::Column::const_iterator k;
                        std::vector<uint64_t>::iterator l;

                        uint64_t t;

                        if (_tmp.size () < w.size ())
                                _tmp.resize (w.size ());

                        std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

                        for (j = v.begin (); j != v.end (); ++j, ++i)
                                {
                                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                                                {
                                                        t = ((uint64_t) *k) * ((uint64_t) *j);

                                                        *l += t;

                                                        if (*l < t)
                                                                *l += VD.faxpy()._two_64;
                                                }
                                }

                        typename Vector1::iterator w_j;

                        for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                                *w_j = *l % VD.field ().characteristic();

                        return w;
                }
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const
                {
                        linbox_check (A.coldim () == v.size ());
                        linbox_check (A.rowdim () == w.size ());

                        typename Matrix::ConstColIterator i = A.colBegin ();
                        typename Vector2::const_iterator j;
                        typename Matrix::Column::const_iterator k;
                        std::vector<uint64_t>::iterator l;

                        uint64_t t;

                        if (_tmp.size () < w.size ())
                                _tmp.resize (w.size ());

                        std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

                        for (j = v.begin (); j != v.end (); ++j, ++i)
                                {
                                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                                                {
                                                        t = ((uint64_t) k->second) * ((uint64_t) *j);

                                                        _tmp[k->first] += t;

                                                        if (_tmp[k->first] < t)
                                                                _tmp[k->first] += VD.faxpy()._two_64;
                                                }
                                }

                        typename Vector1::iterator w_j;

                        for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                                *w_j = *l % VD.field ().characteristic();

                        return w;
                }

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const
                {

                        linbox_check (A.coldim () == v.size ());
                        linbox_check (A.rowdim () == w.size ());

                        typename Matrix::ConstColIterator i = A.colBegin ();
                        typename Vector2::const_iterator j;
                        typename Matrix::Column::const_iterator k;
                        std::vector<uint64_t>::iterator l;

                        uint64_t t;

                        if (_tmp.size () < w.size ())
                                _tmp.resize (w.size ());

                        std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

                        for (j = v.begin (); j != v.end (); ++j, ++i)
                                {
                                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                                                {
                                                        t = ((uint64_t) k->second) * ((uint64_t) *j);

                                                        _tmp[k->first] += t;

                                                        if (_tmp[k->first] < t)
                                                                _tmp[k->first] += VD.faxpy()._two_64;
                                                }
                                }

                        typename Vector1::iterator w_j;

                        for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                                *w_j = *l % VD.field ().characteristic();

                        return w;
                }

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const
                {

                        linbox_check (A.coldim () == v.size ());
                        linbox_check (A.rowdim () == w.size ());

                        typename Matrix::ConstColIterator i = A.colBegin ();
                        typename Vector2::const_iterator j;
                        typename Matrix::Column::first_type::const_iterator k_idx;
                        typename Matrix::Column::second_type::const_iterator k_elt;
                        std::vector<uint64_t>::iterator l;

                        uint64_t t;

                        if (_tmp.size () < w.size ())
                                _tmp.resize (w.size ());

                        std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

                        for (j = v.begin (); j != v.end (); ++j, ++i)
                                {
                                        for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
                                             k_idx != i->first.end ();
                                             ++k_idx, ++k_elt, ++l)
                                                {
                                                        t = ((uint64_t) *k_elt) * ((uint64_t) *j);

                                                        _tmp[*k_idx] += t;

                                                        if (_tmp[*k_idx] < t)
                                                                _tmp[*k_idx] += VD.faxpy()._two_64;
                                                }
                                }

                        typename Vector1::iterator w_j;

                        for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                                *w_j = *l % VD.field ().characteristic();

                        return w;
                }


		mutable std::vector<uint64_t> _tmp;
	};
}

#undef LINBOX_MAX_INT64



#endif //__LINBOX_modular_int64_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
