/* -*- mode: c; style: linux -*- */

/* linbox/blackbox/dense-matrix.h
 * Copyright (C) 2001 B. David Saunders,
 *               2001 Bradford Hovinen
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __DENSE_MATRIX_H
#define __DENSE_MATRIX_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/blackbox/archetype.h"
#include "linbox/field/vector-domain.h"

namespace LinBox
{
	/** Blackbox sparse matrix template. This is a class of dense matrices
	 * templatized by the {@link Fields field} in which the elements
	 * reside. The matrix is stored as an STL vector of STL vectors of field
	 * elements.
	 *
	 * The class is conforms to the {@link Archetypes archetype} for
	 * \Ref{Blackbox Matrices}.
	 *
	 * The class inherits from Blackbox_archetype, which is an abstract base
	 * class which ensures it adheres to the common object interface.
	 *
	 * @param Field \Ref{LinBox} field
	 */

	template <class Field>
	class DenseMatrix
		: public Blackbox_archetype< std::vector<typename Field::element> >
	{
	    public:
		typedef typename Field::element        element;
		typedef std::vector<element>           Vector;
		typedef std::vector<element>::iterator pointer;
      
		/** Constructor.
		 * @param  F the field of entries; passed so that a possible paramter 
		 *           such as a modulus is known to the matrix.
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		DenseMatrix (Field &F, size_t m, size_t n)
			: _F (F), _VD (F)
		{
			int i;

			linbox_check (m > 0);
			linbox_check (n > 0);

			_rep.resize (0);

			for (; m; m--)
				_rep.push_back (Vector (n));
		}

		/** Copy constructor
		 */
		DenseMatrix (const DenseMatrix &M)
			: _F (M._F), _rep (M._rep), _VD (M._F)
		{}

		Blackbox_archetype<Vector> *clone () const 
		{
			return new DenseMatrix (*this);
		}

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& apply (Vector& y, const Vector& x) const
		{
			int i;

			for (i = 0; i < _rep.size (); i++)
				_VD.dotprod (y[i], _rep[i], x);

			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& applyTranspose (Vector& y, const Vector& x) const
			{ /* FIXME */ return y; }

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
			{ return _rep.size();}
      
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const
			{ return _rep[0].size(); }

		/** Set the entry at (i, j)
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry (size_t i, size_t j, element a_ij) 
			{ _rep[i][j] = a_ij; }

		/** Read the matrix from an input stream
		 * @param file Input stream from which to read
		 */
		void read (std::istream &file)
		{
			int i;
			pointer p;

			for (i = 0; i < rows; ++i) {
				for (p = _res[i].begin (); p != _res[i].end (); ++p) {
					file.ignore (1);
					_VD.read (file, *p);
				}
			}
		}

		/** Write the matrix to an output stream
		 * @param os Output stream to which to write
		 */
		std::ostream &write(std::ostream &os = std::cout)
		{
			int i;
			pointer p;

			for (i = 0; i < _res.size (); ++i) {
				for (p = _res[i].begin (); p != _res[i].end(); ++p) {
					_VD.write (os, *p);
					os << " ";
				}

				os << endl;
			}	
		}

	    private:

		std::vector<Vector>                  _rep;
		Field                               &_F;
		VectorDomain<Field, Vector, Vector>  _VD;
	};
}

#endif
