/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/permutation.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __PERMUTATION_H
#define __PERMUTATION_H

#include <vector>

#include "linbox/blackbox/archetype.h"

#include "linbox/util/debug.h"
#include "linbox-config.h"

#ifdef XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <string>

using std::istream;
using std::ostream;
using std::string;

#endif


// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Blackbox permutation matrix.
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Vector>
	class Permutation : public BlackboxArchetype<Vector>
	{
	    public:

		typedef BlackboxArchetype<Vector> Blackbox;

		/** Constructor from a vector of indices
		 * This constructor creates a permutation matrix based on a vector of indices
		 * @param indices Vector of indices representing the permutation
		 */
		Permutation (vector<int> &_indices)
			: _indices (indices)
		{}

		/** Constructor from a dimension
		 * This constructor creates an n x n permutation matrix, initialized to be the identity
		 * @param n The dimension of hte matrix to create
		 */
		Permutation (int n)
		{
			int i;

			_indices.resize (n);

			for (i = 0; i < n; i++)
				_indices[i] = i;
		}

		/* Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Permutation (const Permutation<Vector> &M)
			: _indices (M._indices)
		{}

#ifdef XMLENABLED
		Permutation(Reader &R)
		{
			if(!R.expectTagName("MatrixOver")) return;
                        if(!R.expectChildTag()) return;
                        R.traverseChild();

                        if(!R.expectTagName("permutation") || !R.expectTagNumVector(_indices)) return;

			R.upToParent();
                        return;
		}
#endif


		// Destructor
		~Permutation (void) {}

		/* Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Permutation *clone () const
			{ return new Permutation (*this); }

		/* Application of BlackBox permutation matrix.
		 * y= P*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		/// #y \leftarrow Px#.
		inline Vector &apply (Vector &y, const Vector &x) const
		{
			size_t i;

			linbox_check (x.size () == _indices.size ());

			// resizing y is now forbidden - bds //y.resize (x.size ());
			linbox_check (y.size () == _indices.size ());

			for (i = 0; i < _indices.size (); i++)
				y[_indices[i]] = x[i];

			return y;
		}

		/* Application of BlackBox permutation matrix transpose.
		 * y= transpose(P)*x, equivalently y= P^-1*x
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		/// #y^T \leftarrow x^T P#.
		inline Vector &applyTranspose (Vector &y, const Vector &x) const
		{
			size_t i;

			linbox_check (x.size () == _indices.size ());

			// resizing y is now forbidden - bds //y.resize (x.size ());
			linbox_check (y.size () == _indices.size ());

			for (i = 0; i < _indices.size (); i++)
				y[i] = x[_indices[i]];

			return y;
		}

		/* Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			return _indices.size ();
		}
    
		/* Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const 
		{
			return _indices.size ();
		}

		/** Add a transposition to the matrix
		 */
		void permute (size_t row1, size_t row2) 
		{
			linbox_check (row1 >= 0 && row1 < _indices.size ());
			linbox_check (row2 >= 0 && row2 < _indices.size ());

			_swap (_indices[row1], _indices[row2]);
		}

#ifdef XMLENABLED

		ostream &write(ostream &out) const
		{
			Writer W;
			if( toTag(W) ) 
				W.write(out);
		
			return out;
		}

		bool toTag(Writer &W) const
		{
			string s;
			W.setTagName("MatrixOver");
			W.setAttribute("rows", Writer::numToString(s, _indices.size()));
			W.setAttribute("cols", Writer::numToString(s, _indices.size()));
			W.setAttribute("implDetail", "permutation");
			
			W.addTagChild();
			W.setTagName("permutation");
			W.addNumericalList(_indices);
			W.upToParent();

			return true;
		}
#endif


	    private:

		void _swap(int &x, int &y) const
		{
			int temp = x;
			x = y;
			y = temp;
		}

		// Vector of indices
		vector<int> _indices;

	}; // template <Vector> class Permutation

} // namespace LinBox

#endif // __PERMUTATION_H
