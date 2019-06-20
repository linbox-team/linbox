/* Copyright (C) 2010 LinBox
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

/*! @file field/modular/modular-int32.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c int32_t .
 */
#ifndef __LINBOX_modular_int32_H
#define __LINBOX_modular_int32_H


#include <cmath>
#include <givaro/modular-integral.h>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-traits.h"
#include "linbox/ring/modular.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/write-mm.h"

#ifndef LINBOX_MAX_INT /* 2147483647 */
#define LINBOX_MAX_INT INT32_MAX
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

        template<class Compute>
	class FieldAXPY<Givaro::Modular<int32_t,Compute> > {
	public:

		typedef int32_t Element;
		typedef int64_t Abnormal;
		typedef Givaro::Modular<int32_t,Compute> Field;

		FieldAXPY (const Field &F) : _field (&F), _y(0)
		{
			_two_64 = (uint64_t(1) << 32) % uint64_t(F.characteristic());
			_two_64 = (_two_64 * _two_64) % uint64_t(F.characteristic());
		}
		
		FieldAXPY (const FieldAXPY &faxpy) :
			_two_64 (faxpy._two_64), _field (faxpy._field), _y (0)
		{}

		FieldAXPY<Field > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			_two_64 = faxpy._two_64;
			return *this;
		}

		 uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) a * (uint64_t) x;
			_y += t;
			if (_y < t) {
				 _y += _two_64;
				 return _y ;
			}
			else
				return _y;
		}

		 uint64_t& accumulate (const Element &t)
		{
			_y += (uint64_t) t;
			if (_y < (uint64_t)t)
				return _y += _two_64;
			else
				return _y;
		}

		 Element& get (Element &y)
		{
			y = Element (_y % (uint64_t) field().characteristic());
			return y;
		}

		 FieldAXPY &assign (const Element y)
		{
			_y = (uint64_t) y;
			return *this;
		}

		 void reset()
		{
			_y = 0;
		}

		inline const Field & field() const { return *_field; }

	public:
		uint64_t _two_64;

	protected:
		const Field * _field;
		uint64_t _y;
	};


	template <class Compute>
	class DotProductDomain<Givaro::Modular<int32_t,Compute> > : public VectorDomainBase<Givaro::Modular<int32_t,Compute> > {

	public:
		typedef int32_t Element;
		typedef Givaro::Modular<int32_t,Compute> Field;
		DotProductDomain(){}
		DotProductDomain (const Field&F) :
			VectorDomainBase<Field> (F)
		{}

		using VectorDomainBase<Givaro::Modular<int32_t,Compute>>::faxpy;
		using VectorDomainBase<Givaro::Modular<int32_t,Compute>>::field;


	protected:
		template <class Vector1, class Vector2>
		 Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
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
					y += (uint64_t) faxpy()._two_64;
			}

			y %= (uint64_t) field().characteristic();
			return res = Element(y);

		}

		template <class Vector1, class Vector2>
		 Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
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
					y += (uint64_t) faxpy()._two_64;
			}


			y %= (uint64_t) field().characteristic();

			return res = (Element) y;
		}
	};

	// Specialization of MVProductDomain for int32_t modular field

	template <class Compute>
	class MVProductDomain<Givaro::Modular<int32_t,Compute> > {
	public:

		typedef int32_t Element;
		typedef Givaro::Modular<int32_t,Compute> Field;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		 Vector1 &mulColDense
		(const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized
			(VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
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

                        std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

                        for (j = v.begin (); j != v.end (); ++j, ++i)
                                {
                                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                                                {
                                                        t = ((uint64_t) *k) * ((uint64_t) *j);

                                                        *l += t;

                                                        if (*l < t)
                                                                *l += (uint64_t) VD.faxpy ()._two_64;
                                                }
                                }

                        typename Vector1::iterator w_j;
                        typedef typename Vector1::value_type elements ;

                        for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                                *w_j = elements(*l % VD.field ().characteristic());

                        return w;
                }

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const
                {
                        linbox_check (A.coldim () == v.size ());
                        linbox_check (A.rowdim () == w.size ());

                        typename Matrix::ConstColIterator       i = A.colBegin ();
                        typename Vector2::const_iterator        j;
                        typename Matrix::Column::const_iterator k;
                        std::vector<uint64_t>::iterator         l;

                        uint64_t t;

                        if (_tmp.size () < w.size ())
                                _tmp.resize (w.size ());

                        std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

                        for (j = v.begin (); j != v.end (); ++j, ++i) {
                                for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
                                        t = ((uint64_t) k->second) * ((uint64_t) *j);

                                        _tmp[k->first] += t;

                                        if (_tmp[k->first] < t)
                                                _tmp[k->first] += (uint64_t)VD.faxpy ()._two_64;
                                }
                        }

                        typename Vector1::iterator w_j;
                        typedef typename Vector1::value_type val_t;

                        for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                                *w_j = (val_t)( (int32_t)(*l) % VD.field ().characteristic() );

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

                        std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

                        for (j = v.begin (); j != v.end (); ++j, ++i)
                                {
                                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                                                {
                                                        t = ((uint64_t) k->second) * ((uint64_t) *j);

                                                        _tmp[k->first] += t;

                                                        if (_tmp[k->first] < t)
                                                                _tmp[k->first] += VD.faxpy ()._two_64;
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

                        std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

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
                                *w_j = *l % VD.field().characteristic();

                        return w;
                }


		mutable std::vector<uint64_t> _tmp;
	};
}

#endif //__LINBOX_modular_int32_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
