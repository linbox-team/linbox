/* linbox/blackbox/permutation.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
#ifndef __LINBOX_PERMUTATION_STORAGE
// #include "linbox/vector/light_container.h"
// #define __LINBOX_PERMUTATION_STORAGE LightContainer< long >
#include <vector>
#define __LINBOX_PERMUTATION_STORAGE std::vector< long >
#endif

#include "linbox-config.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/randiter/mersenne-twister.h"


// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief size is n.

	  \ingroup blackbox
	 * @param Storage \ref LinBox dense or sparse vector of field elements
	 */
	template<class _Field, class _Storage = __LINBOX_PERMUTATION_STORAGE >
	class Permutation : public  BlackboxInterface {
		const _Field* _field;
	public:
		typedef Permutation<_Field, _Storage>	Self_t;
		typedef _Storage 			Storage;
		typedef _Field				Field;
		typedef typename Field::Element 	Element;

		/** Constructor from a vector of indices.
		 * This constructor creates a permutation matrix based on a vector of indices
		 * @param F
		 * @param indices Vector of indices representing the permutation
		 */
		Permutation (Storage & indices, const Field& F = Field()) :
			_field(F), _indices (indices)
		{}

		/** Constructor from a dimension.
		 * This constructor creates an n x n permutation matrix, initialized to be the identity
		 * @param n The dimension of the matrix to create
		 * @param F
		 */
		Permutation (int n, const Field& F = Field()) :
			_field(&F)
		{
			identity(n);
		}

		Permutation (const Field& F = Field(), size_t n=0, size_t m = 0) :
			_field(&F)
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

		void random(size_t n)
		{
			identity((int)n);
			MersenneTwister r((unsigned int)time(NULL));
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



		template<typename _Tp1>
		struct rebind {
			typedef Permutation<_Tp1, Storage> other;
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
		 * Add a transposition to the matrix
		*/
		void permute (size_t row1, size_t row2)
		{
			linbox_check (/*  row1 >= 0 &&*/ row1 < _indices.size ());
			linbox_check (/*  row2 >= 0 &&*/ row2 < _indices.size ());
			std::swap (_indices[row1], _indices[row2]);

		}

		const Field& field() const { return *_field; }

		std::ostream &write(std::ostream &os) const //, LINBOX_enum(Tag::FileFormat) format = Tag::FileFormat::Maple) const
		{
			// 		for (typename Storage::const_iterator it=_indices.begin(); it!=_indices.end(); ++it)
			//                     std::cerr << *it << ' ';
			typename Field::Element one, zero; field().init(one,1UL);field().init(zero,0UL);
			os << "[";
			bool firstrow=true;
			long nmu = (long)_indices.size()-1;
			for (typename Storage::const_iterator it=_indices.begin(); it!=_indices.end(); ++it) {
				if (firstrow) {
					os << "[";
					firstrow =false;
				}
				else
					os << ", [";

				long i=0;
				for( ; i< *it ; ++i) {
					field().write(os, zero);
					if (i < nmu) os << ',';
				}
				field().write(os, one);
				if (i < nmu) os << ',';
				for(++i ; i< static_cast<long>(_indices.size()) ; ++i) {
					field().write(os, zero);
					if (i < nmu) os << ',';
				}
				os << " ]";
			}

			return os << "]";
		}

		//!@bug there is no read here. (needed by test-blackbox.h)

		Storage& setStorage(const Storage& s) { return _indices=s; }
		const Storage& getStorage() const { return _indices; }

		/// Generate next permutation in lex order.
		void next() {
			int n = _indices.size();
			if (n == 1) return;
			int i, j;
			for (i = n-2; i >= 0 and _indices[(size_t)i] >= _indices[(size_t)i+1]; --i);
			if (i < 0) {identity(n); return; }
			for (j = i+2; j < n and _indices[(size_t)i] <= _indices[(size_t)j]; ++j);
			std::swap(_indices[(size_t)i], _indices[(size_t)j-1]);
			reverse(_indices.begin() + i + 1, _indices.end());
		}
	private:
		// Vector of indices
		Storage _indices;

	}; // template <Vector> class Permutation

} // namespace LinBox

#endif // __LINBOX_bb_permutation_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

