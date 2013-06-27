/* linbox/matrix/blas-matrix-multimod.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
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

/*! @file matrix/blas-matrix-multimod.h
 * @ingroup matrix
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field::Element</code>s. It also has the BlasBlackbox interface.
 *
 */

/*! @file matrix/blas-matrix-multimod.h
 * @brief specialisation for mutlimod field.
 */

#ifndef __LINBOX_blas_matrix_multimod_H
#define __LINBOX_blas_matrix_multimod_H

#include "linbox/util/debug.h"

// #include "linbox/vector/subiterator.h"
// #include "linbox/vector/subvector.h"

#include "linbox/matrix/matrix-category.h"
#include "linbox/linbox-tags.h"
#include "linbox/matrix/blas-matrix.h"

THIS_CODE_MAY_NOT_COMPILE_AND_IS_NOT_TESTED

namespace LinBox
{ /*  Specialisation of BlasMatrix for MultiModDouble field */

	/*! No Doc.
	*/
	template<>
	class BlasMatrix<MultiModDouble> {

	public:

		typedef MultiModDouble         Field;
		typedef std::vector<double>  Element;
		typedef BlasMatrix<MultiModDouble> Self_t;

	protected:

		MultiModDouble                 _field;
		const std::vector<MatrixDomain<Modular<double> > >   _MD;
		size_t                  _row,_col;
		Element                _One,_Zero;
		std::vector<BlasMatrix<Modular<double> >* > _rep;
		std::vector<double>       _entry;
	public:


		//BlasMatrix () {}

		BlasMatrix (const MultiModDouble& F) :
			_field(F) , _rep(F.size()), _entry(F.size())
		{}

		BlasMatrix (const Field& F, size_t m, size_t n, bool alloc=true) :
			_field(F), _row(m) , _col(n) , _rep(F.size()),  _entry(F.size())
		{
			for (size_t i=0;i<_rep.size();++i)
				_rep[i] =  new BlasMatrix<Modular<double> > (F.getBase(i), m, n);
		}

		BlasMatrix (const BlasMatrix<MultiModDouble> & A):
			_field(A._field),_row(A._row), _col(A._col),
			_rep(A._rep.size()), _entry(A._entry)
		{

			for (size_t i=0;i<_rep.size();++i)
				_rep[i]= new  BlasMatrix<Modular<double> > (const_cast<BlasMatrix<Modular<double> >& >( *A._rep[i]));
		}


		const BlasMatrix<MultiModDouble>& operator=(const BlasMatrix<MultiModDouble> & A)
		{
			_field   = A._field;
			_row = A._row;
			_col = A._col;
			_rep = std::vector<BlasMatrix<Modular<double> >* >(A._rep.size());
			_entry = A._entry;
			for (size_t i=0;i<_rep.size();++i)
				_rep[i]= new  BlasMatrix<Modular<double> > (const_cast<BlasMatrix<Modular<double> >& >( *A._rep[i]));
			return *this;
		}


		~BlasMatrix() {for (size_t i=0; i< _rep.size();++i) {delete _rep[i];} }

		template <class Vector1, class Vector2>
		Vector1&  apply (Vector1& y, const Vector2& x) const
		{
			for (size_t i=0;i<_rep.size();++i) {
				std::vector<double> x_tmp(x.size()), y_tmp(y.size());
				for (size_t j=0;j<x.size();++j)
					x_tmp[j]= x[j][i];

				_rep[i]->apply(y_tmp, x_tmp);

				for (size_t j=0;j<y.size();++j){
					y[j][i]=y_tmp[j];

				}
			}

			return y;
		}

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const
		{
			for (size_t i=0;i<_rep.size();++i) {
				std::vector<double> x_tmp(x.size()), y_tmp(y.size());
				for (size_t j=0;j<x.size();++j)
					x_tmp[i]= x[j][i];

				_rep[i]->applyTranspose(y_tmp, x_tmp);

				for (size_t j=0;j<y.size();++j)
					y[j][i]=y_tmp[i];
			}

			return y;
		}

#if 0
		template<typename _Tp1>
		struct rebind
		{
			typedef BlasMatrix<_Tp1> other;

			void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
				Ap = new other(F, A.rowdim(), A.coldim());
				Hom<Field, _Tp1> hom(A. field(), F);

				hom.image (*Ap_p, *A_p);
			}
		};
#endif

		size_t rowdim() const {return _row;}

		size_t coldim() const {return _col;}


		const Field &field() const  {return _field;}


		std::ostream& write(std::ostream& os) const
		{
			for (size_t i=0;i<_rep.size();++i)
				_rep[i]->write(os);
			return os;
		}


		void setEntry (size_t , size_t j, const Element &a_ij)
		{
			for (size_t i=0; i< _rep.size();++i)
				_rep[i]->setEntry(i,j,a_ij[i]);
		}


		const Element& getEntry (size_t , size_t j)
		{
			for (size_t i=0; i< _rep.size();++i)
				_entry[i]=_rep[i]->getEntry(i,j);
			return _entry;
		}

		BlasMatrix<Modular<double> >*& getMatrix(size_t i) {return _rep[i];}

	};


	template <>
	class MatrixContainerTrait<BlasMatrix<MultiModDouble> > {
	public:
		typedef MatrixContainerCategory::Blackbox Type;
	};
} // LinBox

#endif // __LINBOX_blas_matrix_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
