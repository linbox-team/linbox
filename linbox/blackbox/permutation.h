/* linbox/blackbox/permutation.h
 * Copyright (C) 2001 LinBox group
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *            Jean-Guillaume Dumas <jean-guillaume.dumas@imag.fr>
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

#ifndef __LINBOX_bb_permutation_H
#define __LINBOX_bb_permutation_H

#include <utility>
#include <algorithm>

#include "linbox/linbox-config.h"
#include "linbox/linbox-tags.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/blackbox/fibb.h"

namespace LinBox
{

	/** \brief

	  \ingroup blackbox
	 */
	template<class _Field, class _Matrix=DenseMatrix<_Field>>
	class Permutation : public  FIBB<_Field>
	{
	public:
		typedef _Field					Field;
		typedef typename Field::Element	Element;
		typedef Permutation<Field>		Self_t;
		typedef LightContainer<long>	Storage;
		typedef _Matrix					Matrix;
	protected:
		Storage _indices; // Vector of indices
		const Field& _field;
	public:

		/** Constructor from a vector of indices.
		 * This constructor creates a permutation matrix based on a vector of indices
		 * @param F
		 * @param indices Vector of indices representing the permutation
		 * Permutation P has 1 in the P_{i, _indices[i]} positions.
		 */
		Permutation (Storage & indices, const Field& F) :
			_field(F), _indices (indices)
		{}

		Self_t& init(size_t* P, size_t n)
		{	_indices.resize(n);
			for (size_t i = 0; i < n; ++i)
				_indices[i] = P[i];
			return *this;
		}

		Permutation (size_t* P, size_t n, const Field& F)
		: _field(F) {
            init(P, n);
        }

		/** \brief n x n permutation matrix, initially the identity.
		 * @param n The dimension of the matrix
		 * @param F field or ring
		 */
		Permutation (int n, const Field& F) : _field(F)
		{
			identity(n);
		}

		Permutation (const Field& F, size_t n=0, size_t m = 0) : _field(F)
		{
			identity((int)n);
		}



		//!@bug should be size_t
		void identity(int n)
		{
			this->_indices.resize ((size_t)n);
			for (typename Storage::value_type i=0; i < n; ++i)
				_indices[(size_t)i] = i;
		}

		void cyclicShift(size_t n)
		{
			this->_indices.resize ((size_t)n);
			for (typename Storage::value_type i=0; i < n; ++i)
				_indices[(size_t)i] = i+1 % n;
		}

        void random(unsigned int seed= static_cast<unsigned int>(std::time(nullptr)))
		{
			size_t n = rowdim();
			identity((int)n);
			MersenneTwister r(seed);
			// Knuth construction
			for (size_t i = 0; i < n-1; ++i) {
				size_t j = i + r.randomInt()%(n-i);
				std::swap(_indices[(size_t)i], _indices[(size_t)j]);
			}
		}


		/* Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Permutation (const Permutation &Mat) :
			_field(Mat._field),_indices (Mat._indices)
		{}

		// Destructor
		~Permutation (void) {}

		/** Application of BlackBox permutation matrix.
		  \f$y \leftarrow Px\f$.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template<class OutVector, class InVector>
		inline OutVector &apply (OutVector &y, const InVector &x) const
		{
			size_t i;

			linbox_check (x.size () == _indices.size ());

			// resizing y is now forbidden - bds //y.resize (x.size ());
			linbox_check (y.size () == _indices.size ());

			for (i = 0; i < x.size(); ++i)
				field().assign(y[(size_t)i], x[(size_t)_indices[(size_t)i]]);

			return y;
		}

		/** Application of BlackBox permutation matrix transpose.
		 * <code>y= transpose(P)*x</code>, equivalently <code>y= P^-1*x</code>
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		/// \f$y^T \leftarrow x^T P\f$.
		template<class OutVector, class InVector>
		inline OutVector &applyTranspose (OutVector &y, const InVector &x) const
		{
			size_t i;

			linbox_check (x.size () == _indices.size ());

			// resizing y is now forbidden - bds //y.resize (x.size ());
			linbox_check (y.size () == _indices.size ());

			for (i = 0; i < _indices.size (); ++i)
				field().assign(y[(size_t)_indices[(size_t)i]], x[(size_t)i]);

			return y;
		}

            // Permutes the rows of X
		Matrix& applyRight(Matrix& Y, const Matrix& X) const
		{
			Element x; field().init(x);
			for (size_t i = 0; i < Y.rowdim(); ++i){
				size_t k = _indices[i];
				for (size_t j = 0; j < Y.coldim(); ++j)
					Y.setEntry(i,j, X.getEntry(x, k, j));
			}
		/* desired form
			for (size_t i = 0; i < rowdim(); ++i)
			{
				Matrix Yrow(Y, i, 0, 1, Y.coldim());
				Matrix Xrow(X, _indices[i], 0, 1, X.coldim());
				Yrow.copy(Xrow); // right kind of copy?
			}
		*/
			return Y;
		}

            // Permutes the columns of X
		Matrix& applyLeft(Matrix& Y, const Matrix& X) const
		{
			Element x; field().init(x);
			for (size_t j = 0; j < Y.coldim(); ++j){
				size_t k = _indices[j];
				for (size_t i = 0; i < Y.rowdim(); ++i)
					Y.setEntry(i,k, X.getEntry(x, i, j));
			}
		/* desired form
			for (size_t i = 0; i < coldim(); ++i)
			{
				Matrix Ycol(Y, 0, _indices[i], Y.rowdim(), 1);
				Matrix Xcol(X, 0, i, X.rowdim(), 1);
				Ycol.copy(Xcol);
			}
		*/
			return Y;
		}
		/* FIBB functions */

		BBType bbTag() const { return permutation; }

		size_t& rank(size_t& r) const
		{ return r = rowdim(); }

		Element& det(Element& d) const
		{	size_t b = 0, i, j, k;
			Storage marks(_indices.size());
			for (i = 0; i < _indices.size(); ++i)
			if (not marks[i])
			{	for (k = 1, j = _indices[i]; i != j; ++k, j = _indices[j])
				;
				b &= k;
			}
			return d = b&1 ? field().mOne : field().one;
		}

            // Inverse permutation on the rows
		Matrix& solveRight(Matrix& Y, const Matrix& X) const
		{	Element x; field().init(x);
			for (size_t i = 0; i < Y.rowdim(); ++i){
				size_t k = _indices[i];
				for (size_t j = 0; j < Y.coldim(); ++j)
					Y.setEntry(k,j, X.getEntry(x, i, j));
			}
		/* desired form
			for (size_t i = 0; i < rowdim(); ++i)
			{
				Matrix Yrow(Y, _indices[i], 0, 1, Y.coldim());
				Matrix Xrow(X, i, 0, 1, X.coldim());
				Yrow.copy(Xrow);
			}
		*/
			return Y;
		}

            // Inverse permutation on the columns
		Matrix& solveLeft(Matrix& Y, const Matrix& X) const
		{	Element x; field().init(x);
			for (size_t j = 0; j < Y.coldim(); ++j){
				size_t k = _indices[j];
				for (size_t i = 0; i < Y.rowdim(); ++i)
					Y.setEntry(i,j, X.getEntry(x, i, k));
			}
		/* desired form
			for (size_t i = 0; i < coldim(); ++i)
			{
				Matrix Ycol(Y, 0, i, Y.rowdim(), 1);
				Matrix Xcol(X, 0, _indices[i], X.rowdim(), 1);
				Ycol.copy(Xcol);
			}
		*/
			return Y;
		}
		Matrix& nullspaceRandomRight(Matrix& N) const
		{	N.zero(); return N; }
		Matrix& nullspaceRandomLeft(Matrix& N) const
		{	N.zero(); return N; }
		Matrix& nullspaceBasisRight(Matrix& N) const
		{	N.resize(rowdim(), 0); return N; }
		Matrix& nullspaceBasisLeft(Matrix& N) const
		{	N.resize(0, coldim()); return N; }
		/* end FIBB section */

		template<typename _Tp1>
		struct rebind {
			typedef Permutation<_Tp1> other;
			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
				Ap->setStorage( A.getStorage() );
			}
		};



		/* Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		/// rowdim
		size_t rowdim (void) const
		{
			return _indices.size ();
		}

		/* Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		/// coldim
		size_t coldim (void) const
		{
			return _indices.size ();
		}

		/**
		 * this <-- transposition(i,j)*this
		 * indices (i = 0, j = 1):
		 * 0 2 1 * 1 0 2 = 1 2 0
		 * matrices (corresponding):
		 * 1 0 0   0 1 0   0 1 0
		 * 0 0 1 * 1 0 0 = 0 0 1
		 * 0 1 0   0 0 1   1 0 0
		 */
		void permute (size_t i, size_t j)
		{
			linbox_check (/*  i >= 0 &&*/ i < _indices.size ());
			linbox_check (/*  j >= 0 &&*/ j < _indices.size ());
			std::swap (_indices[i], _indices[j]);
		}

        size_t operator[](size_t i) const {
            return _indices[i];
        }

		const Field& field() const { return _field; }

		//!@bug needs a read. (needed by test-blackbox.h)
		std::istream &read(std::istream &os)
		{ return read(os, Tag::FileFormat::Plain); }

		//!@bug needs a MM version
		std::ostream &write(std::ostream &os) const
		{ return write(os, Tag::FileFormat::Plain); }

		std::ostream &write(std::ostream &os, Tag::FileFormat format) const
		{

			// Avoid unneeded overhead in the case that this
			// printing is disabled
			if (not os)
				return os;

			switch (format) {
			case Tag::FileFormat::Maple:
				{
					os << '[';
					bool firstrow=true;
					long nmu = (long)_indices.size()-1;
					for (typename Storage::const_iterator it=_indices.begin(); it!=_indices.end(); ++it) {
						if (firstrow) {
							os << '[';
							firstrow =false;
						}
						else
							os << ", [";

						long i=0;
						for( ; i< *it ; ++i) {
							field().write(os, field().zero);
							if (i < nmu) os << ',';
						}
						field().write(os, field().one);
						if (i < nmu) os << ',';
						for(++i ; i< static_cast<long>(_indices.size()) ; ++i) {
							field().write(os, field().zero);
							if (i < nmu) os << ',';
						}
						os << ']';

					}
					os << ']';
					break;
				}
			case Tag::FileFormat::Pretty:
				{
					for (typename Storage::const_iterator it=_indices.begin(); it!=_indices.end(); ++it) {
						os << "  [";

						long i=0;
						for( ; i< *it ; ++i) {
							field().write(os << ' ', field().zero);
						}
						field().write(os << ' ', field().one);
						for(++i ; i< static_cast<long>(_indices.size()) ; ++i) {
							field().write(os << ' ', field().zero);
						}
						os << " ]" << std::endl;

					}
					break;
				}

			default:
				os << '{';
				for (typename Storage::const_iterator it=_indices.begin(); it!=_indices.end(); ++it)
					os << *it << ' ';
				os << '}';
				break;


			}



			return os;
		}

		//!@bug there is no read here. (needed by test-blackbox.h)
		std::istream &read(std::istream &is, Tag::FileFormat format)
		{
            switch (format) {
                case Tag::FileFormat::Plain:
                {
                    char t;
                    is >> t;
                    Storage::value_type val;
                    _indices.resize(0);
                    while( t != '}') {
                        is >> val;
                        _indices.push_back(val);
                        is >> t;
                        if (t!='}') is.putback (t);
                    }
                    break;

                }
                default:
                    throw NotImplementedYet();
            }
			return is ;
		}


	Storage& setStorage(const Storage& s) { return _indices=s; }
	const Storage& getStorage() const { return _indices; }
	Storage& getStorage() { return _indices; }

	/// Generate next permutation in lex order.
	void next()
	{
		int n = _indices.size();
		if (n == 1) return;
		int i, j;
		for (i = n-2; i >= 0 and _indices[(size_t)i] >= _indices[(size_t)i+1]; --i);
		if (i < 0) {identity(n); return; }
		for (j = i+2; j < n and _indices[(size_t)i] <= _indices[(size_t)j]; ++j);
		std::swap(_indices[(size_t)i], _indices[(size_t)j-1]);
		std::reverse(_indices.begin() + i + 1, _indices.end());
	}

}; // template <Vector> class Permutation

} // namespace LinBox

#endif // __LINBOX_bb_permutation_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
