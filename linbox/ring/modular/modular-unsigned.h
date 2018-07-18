/* linbox/ring/modular.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 * Copyright (C) 2011 LinBox
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * LargeModular is now replace by a class Givaro::Modular parameterized on the element
 * type. So, the old LargeModular is equivalent to Givaro::Modular<integer>. All other
 * interface details are exactly the same.
 *
 * Renamed from large-modular.h to modular.h
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

#ifndef __LINBOX_field_modular_unsigned_H
#define __LINBOX_field_modular_unsigned_H

//Dan Roche 7-2-04
#ifndef __LINBOX_MIN
#define __LINBOX_MIN(a,b) ( (a) < (b) ? (a) : (b) )
#endif

namespace LinBox { /*  uint8_t */

            /*! Specialization of FieldAXPY for uint8_t modular field */
    
    template<class Compute_t>
	class FieldAXPY<Givaro::Modular<uint8_t,Compute_t > > {
	public:

		typedef uint8_t Element;
		typedef uint64_t Abnormal;
		typedef Givaro::Modular<uint8_t, Compute_t> Field;

		FieldAXPY (const Field &F) :
                _k (((uint64_t) -1LL) / ((F.characteristic() - 1) * (F.characteristic() - 1))),
                _field (&F),
                _y (0),
                i (_k)
            {
            }

		FieldAXPY (const FieldAXPY &faxpy) :
                _k (faxpy._k),
                _field (faxpy._field),
                _y (0),
                i (_k)
            {}

		FieldAXPY<Field> &operator = (const FieldAXPY &faxpy)
            {
                _field = faxpy._field;
                _y = faxpy._y;
                _k = faxpy._k;
                return *this;
            }
        
		inline uint64_t& mulacc (const Element &a, const Element &x)
            {
                uint32_t t = (uint32_t) a * (uint32_t) x;

                if (!i--) {
                    i = int(_k);
                    return _y = _y % (uint32_t) field().characteristic() + t;
                }
                else
                    return _y += t;
            }

		inline uint64_t& accumulate (const Element &t)
            {

                if (!i--) {
                    i = int(_k);
                    return _y = _y % (uint32_t) field().characteristic() + t;
                }
                else
                    return _y += t;
            }

		inline Element &get (Element &y) const
            {
                const_cast<FieldAXPY<Field>*>(this)->_y %= (uint32_t) field().characteristic();
                if ((int32_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
                y = (uint8_t) _y;
                const_cast<FieldAXPY<Field>*>(this)->i = int(_k);
                return y;
            }

		inline FieldAXPY &assign (const Element y)
            {
                _y = y;
                i = int(_k);
                return *this;
            }

		inline void reset()
            {
                _y = 0;
            }

		inline const Field & field() const { return *_field; }
		
	public:
		uint64_t _k; 
		
	private:
		const Field *_field;
		uint64_t _y;
		int64_t i;
	};

        //! Specialization of DotProductDomain for unsigned short modular field

	template <class Compute_t>
	class DotProductDomain<Givaro::Modular<uint8_t, Compute_t> > : public  VectorDomainBase<Givaro::Modular<uint8_t, Compute_t> > {
	public:

		typedef uint8_t Element;
		typedef Givaro::Modular<uint8_t, Compute_t> Field;

		DotProductDomain(){}
		DotProductDomain (const Field &F) :
                VectorDomainBase<Field> (F)
            {}
		using VectorDomainBase<Field>::field;
		using VectorDomainBase<Field>::faxpy;

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
            {
                typename Vector1::const_iterator i = v1.begin ();
                typename Vector2::const_iterator j = v2.begin ();

                typename Vector1::const_iterator iterend = v1.begin () + (ptrdiff_t)(v1.size() % faxpy()._k);

                uint64_t y = 0;

                for (; i != iterend; ++i, ++j)
                    y += (uint64_t) *i * (uint64_t) *j;

                y %= (uint64_t) field().characteristic();

                for (; iterend != v1.end (); j += (ptrdiff_t)faxpy()._k) {
                    typename Vector1::const_iterator iter_i = iterend;
                    typename Vector2::const_iterator iter_j;

                    iterend += (ptrdiff_t)faxpy()._k;

                    for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
                        y += (uint64_t) *iter_i * (uint64_t) *j;

                    y %= (uint64_t) field().characteristic();
                }

                return res = (uint8_t) y;
            }

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
            {
                typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
                typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();

                uint64_t y = 0;

                if (v1.first.size () < faxpy()._k) {
                    for (; i_idx != v1.first.end (); ++i_idx, ++i_elt)
                        y += (uint64_t) *i_elt * (uint64_t) v2[*i_idx];

                    return res = uint8_t (y % (uint64_t) field().characteristic());
                }
                else {
                    typename Vector1::first_type::const_iterator iterend = v1.first.begin () +(ptrdiff_t)( v1.first.size() % faxpy()._k);

                    for (; i_idx != iterend; ++i_idx, ++i_elt)
                        y += (uint64_t) *i_elt * (uint64_t) v2[*i_idx];

                    y %= (uint64_t) field().characteristic();

                    while (iterend != v1.first.end ()) {
                        typename Vector1::first_type::const_iterator iter_i_idx = iterend;
                        typename Vector1::second_type::const_iterator iter_i_elt = i_elt;

                        iterend += (ptrdiff_t)faxpy()._k;
                        i_elt += (ptrdiff_t)faxpy()._k;

                        for (; iter_i_idx != iterend; ++iter_i_idx, ++iter_i_elt)
                            y += (uint64_t) *iter_i_elt * (uint64_t) v2[*iter_i_idx];

                        y %= (uint64_t) field().characteristic();
                    }

                    return res = (uint8_t) y;
                }
            }
	};

        //! Specialization of MVProductDomain for uint8_t modular field

	template<class Compute_t>
	class MVProductDomain<Givaro::Modular<uint8_t,Compute_t> > {
	public:

		typedef uint8_t Element;
		typedef Givaro::Modular<uint8_t,Compute_t> Field;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
            {
                return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
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
                typename Vector2::const_iterator j, j_end;
                typename Matrix::Column::const_iterator k;
                std::vector<uint32_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () + (ptrdiff_t)w.size (), 0);

                l_end = _tmp.begin () +(ptrdiff_t) w.size ();

                do {
                    j = v.begin ();
                    j_end = j + __LINBOX_MIN (A->coldim (), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                            *l += *k * *j;

                    j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = *l;

                return w;
            }

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const
            {
                linbox_check (A.coldim () == v.size ());
                linbox_check (A.rowdim () == w.size ());

                typename Matrix::ConstColIterator i = A.colBegin ();
                typename Vector2::const_iterator j, j_end;
                typename Matrix::Column::const_iterator k;
                std::vector<uint32_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () + (ptrdiff_t)w.size (), 0);

                l_end = _tmp.begin () + (ptrdiff_t)w.size ();


                do {
                    j = v.begin ();
                    j_end = j + __LINBOX_MIN (A->coldim (), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                            _tmp[k->first] += k->second * *j;

                    j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = *l;

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
                typename Vector2::const_iterator j, j_end;
                typename Matrix::Column::const_iterator k;
                std::vector<uint32_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () + (ptrdiff_t)w.size (), 0);

                l_end = _tmp.begin () +(ptrdiff_t) w.size ();

                do {
                    j = v.begin ();
                    j_end = j + __LINBOX_MIN (A->coldim (), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                            _tmp[k->first] += k->second * *j;

                    j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = *l;

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
                typename Vector2::const_iterator j, j_end;
                typename Matrix::Column::first_type::const_iterator k_idx;
                typename Matrix::Column::second_type::const_iterator k_elt;
                std::vector<uint32_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () + (ptrdiff_t)w.size (), 0);

                l_end = _tmp.begin () + (ptrdiff_t)w.size ();

                do {
                    j = v.begin ();
                    j_end = j + (ptrdiff_t)__LINBOX_MIN (uint64_t (A.coldim ()), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
                             k_idx != i->first.end ();
                             ++k_idx, ++k_elt, ++l)
                            _tmp[*k_idx] += *k_elt * *j;

                    j_end += (ptrdiff_t) __LINBOX_MIN (uint64_t (A.coldim () - (size_t)(j_end - v.begin ())), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;
                typedef typename Vector1::value_type val_t ;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = (val_t) *l;

                return w;
            }


		mutable std::vector<uint32_t> _tmp;
	};

}

namespace LinBox { /*  uint16_t */

        /*! Specialization of FieldAXPY for uint16_t modular field */
	template<class Compute_t>
	class FieldAXPY<Givaro::Modular<uint16_t,Compute_t> > {
	public:

		typedef uint16_t Element;
		typedef Givaro::Modular<uint16_t,Compute_t> Field;

		FieldAXPY (const Field &F) :
                _k (((uint64_t) -1LL) / ((F.characteristic() - 1) * (F.characteristic() - 1))),
                _field (&F),
                _y (0),
                i (_k)
            {}
		
		FieldAXPY (const FieldAXPY &faxpy) :
                _k (faxpy._k), _field (faxpy._field), _y (0), i (_k)
            {}

		FieldAXPY<Field > &operator = (const FieldAXPY &faxpy)
            {
                _field = faxpy._field;
                _y = faxpy._y;
                _k = faxpy._k;
                return *this;
            }

		inline uint64_t& mulacc (const Element &a, const Element &x)
            {
                uint64_t t = (uint64_t) ((long long) a * (long long) x);

                if (!i--) {
                    i = (int)_k;
                    return _y = _y % (uint64_t) field().characteristic() + t;
                }
                else
                    return _y += t;
            }

		inline uint64_t& accumulate (const Element &t)
            {
                if (!i--) {
                    i = (int)_k;
                    return _y = _y % (uint64_t) field().characteristic() + t;
                }
                else
                    return _y += t;
            }

		inline Element &get (Element &y) const
            {
                const_cast<FieldAXPY<Field>*>(this)->_y %= (uint64_t) field().characteristic();
                if ((int64_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
                y = (uint16_t) _y;
                const_cast<FieldAXPY<Field>*>(this)->i = int(_k);
                return y;
            }

		inline FieldAXPY &assign (const Element y)
            {
                _y = y;
                i = (int)_k;
                return *this;
            }

		inline void reset()
            {
                _y = 0;
            }

		inline const Field & field() const {return *_field;}
		
	public:
		uint64_t _k;
		
	private:
		const Field *_field;
		uint64_t _y;
		int64_t i;
	};

        //! Specialization of DotProductDomain for unsigned short modular field

	template<class Compute_t>
	class DotProductDomain<Givaro::Modular<uint16_t,Compute_t> > : public VectorDomainBase<Givaro::Modular<uint16_t,Compute_t> > {
	public:

		typedef uint16_t Element;
        typedef Givaro::Modular<uint16_t,Compute_t> Field;

		DotProductDomain () {}
		DotProductDomain (const Field &F) :
                VectorDomainBase<Field > (F)
            {}
		using VectorDomainBase<Field>::field;
		using VectorDomainBase<Field>::faxpy;

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
            {
                typename Vector1::const_iterator i = v1.begin ();
                typename Vector2::const_iterator j = v2.begin ();

                typename Vector1::const_iterator iterend = v1.begin () + (ptrdiff_t)(v1.size() % faxpy()._k);

                uint64_t y = 0;

                for (; i != iterend; ++i, ++j)
                    y += (uint64_t) *i * (uint64_t) *j;

                y %= (uint64_t) field().characteristic();

                for (; iterend != v1.end (); j += faxpy()._k) {
                    typename Vector1::const_iterator iter_i = iterend;
                    typename Vector2::const_iterator iter_j;

                    iterend += faxpy()._k;

                    for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
                        y += (uint64_t) *iter_i * (uint64_t) *j;

                    y %= (uint64_t) field().characteristic();
                }

                return res = (uint16_t) y;
            }

        
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
            {
                typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
                typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();

                uint64_t y = 0;

                if (v1.first.size () < faxpy()._k) {
                    for (; i_idx != v1.first.end (); ++i_idx, ++i_elt)
                        y += (uint64_t) *i_elt * (uint64_t) v2[*i_idx];

                    return res = (uint16_t) (y % (uint64_t) field().characteristic());
                }
                else {
                    typename Vector1::first_type::const_iterator iterend = v1.first.begin () +(ptrdiff_t)( v1.first.size() % faxpy()._k );

                    for (; i_idx != iterend; ++i_idx, ++i_elt)
                        y += (uint64_t) *i_elt * (uint64_t) v2[*i_idx];

                    y %= (uint64_t) field().characteristic();

                    while (iterend != v1.first.end ()) {
                        typename Vector1::first_type::const_iterator iter_i_idx = iterend;
                        typename Vector1::second_type::const_iterator iter_i_elt = i_elt;

                        iterend += faxpy()._k;
                        i_elt += faxpy()._k;

                        for (; iter_i_idx != iterend; ++iter_i_idx, ++iter_i_elt)
                            y += (uint64_t) *iter_i_elt * (uint64_t) v2[*iter_i_idx];

                        y %= (uint64_t) field().characteristic();
                    }

                    return res = (Element) y;
                }
            }
        
	};

        //! Specialization of MVProductDomain for uint16_t modular field

	template<class Compute_t>
	class MVProductDomain<Givaro::Modular<uint16_t,Compute_t> > {
	public:

		typedef uint16_t Element;
        typedef Givaro::Modular<uint16_t,Compute_t> Field;
	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
            {
                return mulColDenseSpecialized (VD, w, A, v, VectorTraits<typename Matrix::Column>::VectorCategory ());
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
                typename Vector2::const_iterator j = v.begin (), j_end;
                typename Matrix::Column::const_iterator k;
                    // Dan Roche, 7-1-04
                    // std::vector<uint32_t>::iterator l, l_end;
                std::vector<uint64_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

                l_end = _tmp.begin () +(ptrdiff_t) w.size ();

                do {
                    j = v.begin ();
                    j_end = j + __LINBOX_MIN (A->coldim (), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                            *l += *k * *j;

                    j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = *l;

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
                typename Vector2::const_iterator j, j_end;
                typename Matrix::Column::const_iterator k;
                    // Dan Roche, 7-1-04
                    // std::vector<uint32_t>::iterator l, l_end;
                std::vector<uint64_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

                l_end = _tmp.begin () +(ptrdiff_t) w.size ();

                do {
                    j = v.begin ();
                    j_end = j + __LINBOX_MIN (A->coldim (), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                            _tmp[k->first] += k->second * *j;

                    j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = *l;

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
                typename Vector2::const_iterator j, j_end;
                typename Matrix::Column::const_iterator k;
                    // Dan Roche, 7-1-04
                    // std::vector<uint32_t>::iterator l, l_end;
                std::vector<uint64_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

                l_end = _tmp.begin () +(ptrdiff_t) w.size ();

                do {
                    j = v.begin ();
                    j_end = j + __LINBOX_MIN (A->coldim (), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
                            _tmp[k->first] += k->second * *j;

                    j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = *l;

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
                typename Vector2::const_iterator j, j_end;
                typename Matrix::Column::first_type::const_iterator k_idx;
                typename Matrix::Column::second_type::const_iterator k_elt;
                    // Dan Roche, 7-1-04
                    // std::vector<uint32_t>::iterator l, l_end;
                std::vector<uint64_t>::iterator l, l_end;

                if (_tmp.size () < w.size ())
                    _tmp.resize (w.size ());

                std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

                l_end = _tmp.begin () +(ptrdiff_t) w.size ();

                do {
                    j = v.begin ();
                        //Dan Roche, 7-2-04
                        //j_end = j + __LINBOX_MIN (A->coldim (), VD.faxpy()._k);
                    j_end = j + __LINBOX_MIN (A.coldim (), VD.faxpy()._k);

                    for (; j != j_end; ++j, ++i)
                        for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
                             k_idx != i->first.end ();
                             ++k_idx, ++k_elt, ++l)
                            _tmp[*k_idx] += *k_elt * *j;

                        //j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.faxpy()._k);
                    j_end += __LINBOX_MIN (A.coldim () - (j_end - v.begin ()), VD.faxpy()._k);

                    for (l =_tmp.begin (); l != l_end; ++l)
                        *l %= VD.field ().characteristic();

                } while (j_end != v.end ());

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = *l;

                return w;
            }

		mutable std::vector<uint64_t> _tmp;
	};

}

#include <givaro/modular-integral.h>

namespace LinBox { /*  uint32_t */

	template<class Field>
	class DotProductDomain;
	template<class Field>
	class FieldAXPY;
	template<class Field>
	class MVProductDomain;


        /*! Specialization of FieldAXPY for unsigned short modular field */

	template<class Compute_t>
	class FieldAXPY<Givaro::Modular<uint32_t, Compute_t> > {
	public:

		typedef uint32_t Element;
		typedef Givaro::Modular<uint32_t, Compute_t> Field;

		FieldAXPY (const Field &F) :
                _field (&F), _y(0)
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
                _y += t;

                if (_y < t)
                    return _y += _two_64;
                else
                    return _y;
            }

		inline uint64_t& accumulate_special (const Element &t)
            {
                return _y += t;
            }

		inline Element &get (Element &y) const
            {
                const_cast<FieldAXPY<Field>*>(this)->_y %= (uint64_t) field().characteristic();
                    //if ((int64_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
                return y = (uint32_t) _y;
            }

		inline FieldAXPY &assign (const Element y)
            {
                _y = y;
                return *this;
            }

		inline void reset() {
			_y = 0;
		}

		inline const Field & field() const { return *_field; }
		
	public:
	
		uint64_t _two_64;
		
	private:

		const Field *_field;
		uint64_t _y;
	};

        //! Specialization of DotProductDomain for uint32_t modular field

	template<class Compute_t>
	class DotProductDomain<Givaro::Modular<uint32_t,Compute_t> > : public VectorDomainBase<Givaro::Modular<uint32_t,Compute_t> > {
	public:

		typedef uint32_t Element;
		typedef Givaro::Modular<uint32_t, Compute_t> Field;

		DotProductDomain () {}
		DotProductDomain (const Field &F) :
                VectorDomainBase<Field > (F)
            {}
		using VectorDomainBase<Field >::field;
		using VectorDomainBase<Field >::faxpy;

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
            {
                typename Vector1::const_iterator i;
                typename Vector2::const_iterator j;

                uint64_t y = 0;
                uint64_t t;

                for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
                    t = (uint64_t) *i * (uint64_t) *j;
                    y += t;

                    if (y < t)
                        y += faxpy()._two_64;
                }

                y %= (uint64_t) field().characteristic();

                return res = (uint32_t) y;
            }
        

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
            {
                typename Vector1::first_type::const_iterator i_idx;
                typename Vector1::second_type::const_iterator i_elt;

                uint64_t y = 0;
                uint64_t t = 0;

                for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
                    t = (uint64_t) *i_elt * (uint64_t) v2[*i_idx];
                    y += t;
                    if (y < t)
                        y += faxpy()._two_64;
                }

                y %= (uint64_t) field().characteristic();

                return res = (uint32_t)y;
            }

        
	};

        //! Specialization of MVProductDomain for uint32_t modular field

	template <class Compute_t>
	class MVProductDomain<Givaro::Modular<uint32_t,Compute_t> > {
	public:

		typedef uint32_t Element;
		typedef Givaro::Modular<uint32_t,Compute_t> Field;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
            {
                return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
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

                for (j = v.begin (); j != v.end (); ++j, ++i) {
                    for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
                        t = ((uint64_t) *k) * ((uint64_t) *j);

                        *l += t;

                        if (*l < t)
                            *l += VD.faxpy()._two_64;
                    }
                }

                typename Vector1::iterator w_j;
                typedef typename Vector1::value_type element;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = (element)(*l % VD.field ().characteristic());

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

                std::fill (_tmp.begin (), _tmp.begin () + (ptrdiff_t) w.size (), 0);

                for (j = v.begin (); j != v.end (); ++j, ++i) {
                    for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
                        t = ((uint64_t) k->second) * ((uint64_t) *j);

                        _tmp[k->first] += t;

                        if (_tmp[k->first] < t)
                            _tmp[k->first] += VD.faxpy()._two_64;
                    }
                }

                typename Vector1::iterator             w_j;
                typedef typename Vector1::value_type val_t;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = val_t(*l % VD.field ().characteristic());

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

                for (j = v.begin (); j != v.end (); ++j, ++i) {
                    for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
                        t = ((uint64_t) k->second) * ((uint64_t) *j);

                        _tmp[k->first] += t;

                        if (_tmp[k->first] < t)
                            _tmp[k->first] += VD.faxpy()._two_64;
                    }
                }

                typename Vector1::iterator w_j;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = (uint32_t) (uint32_t)*l % VD.field ().characteristic();

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

                for (j = v.begin (); j != v.end (); ++j, ++i) {
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

                typename Vector1::iterator     w_j;
                typedef typename Vector1::value_type val_t;

                for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
                    *w_j = val_t(*l % VD.field ().characteristic());

                return w;
            }


		mutable std::vector<uint64_t> _tmp;
	};

}


namespace LinBox { /*  uint64_t */

	template<class Field>
	class DotProductDomain;
	template<class Field>
	class FieldAXPY;
	template<class Field>
	class MVProductDomain;

        /*! Specialization of FieldAXPY for unsigned short modular field */

	template<typename Compute_t>
	class FieldAXPY<Givaro::Modular<uint64_t,Compute_t> > {
	public:

		typedef uint64_t Element;
		typedef Givaro::Modular<uint64_t,Compute_t> Field;

		FieldAXPY (const Field &F) :
                _field (&F), _y(0)
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
                return *this;
            }

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
                _y += t;

                if (_y < t)
                    return _y += _two_64;
                else
                    return _y;
            }

		inline uint64_t& accumulate_special (const Element &t)
            {
                return _y += t;
            }

		inline Element &get (Element &y) const
            {
                const_cast<FieldAXPY<Field>*>(this)->_y %= (uint64_t) field().characteristic();
                    //if ((int64_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
                return y = (uint64_t) _y;
            }

		inline FieldAXPY &assign (const Element y)
            {
                _y = y;
                return *this;
            }

		inline void reset() {
			_y = 0;
		}

		inline const Field & field() const { return *_field; }
		
	public:
	
		uint64_t _two_64;
		
	private:

		const Field *_field;
		uint64_t _y;
	};

        //! Specialization of DotProductDomain for uint64_t modular field

	template <typename Compute_t>
	class DotProductDomain<Givaro::Modular<uint64_t,Compute_t>> : public VectorDomainBase<Givaro::Modular<uint64_t,Compute_t> > {
      public:

		typedef uint64_t Element;
		typedef Givaro::Modular<uint64_t,Compute_t> Field;

		DotProductDomain () {}
		DotProductDomain (const Field &F) :
			VectorDomainBase<Field > (F)
		{}
		using VectorDomainBase<Field >::field;
		using VectorDomainBase<Field >::faxpy;

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

        //! Specialization of MVProductDomain for uint64_t modular field

	template <typename Compute_t>
	class MVProductDomain<Givaro::Modular<uint64_t,Compute_t> > {
	public:

		typedef uint64_t Element;
		typedef Givaro::Modular<uint64_t,Compute_t> Field;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
            {
                return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
            }

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Field > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

}


#endif // __LINBOX_field_modular_unsigned_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
