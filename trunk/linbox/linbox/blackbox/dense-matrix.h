/* -*- mode: c; style: linux -*- */

/* linbox/src/blackbox/dense-matrix.h
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

	template <class MatrixDomain>
	class DenseMatrix
		: public Blackbox_archetype< std::vector<typename MatrixDomain::Element> >
	{
	    public:
		typedef typename MatrixDomain::Element Element;
		typedef std::vector<Element>           Vector;
		typedef std::vector<Element>::iterator pointer;
      
		/** Constructor.
		 * @param  _MD the field of entries; passed so that a possible paramter 
		 *             such as a modulus is known to the matrix.
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		DenseMatrix (MatrixDomain &MD, size_t m, size_t n)
			: _MD (MD)
		{
			int i;

			linbox_check (m > 0);
			linbox_check (n > 0);

			_rep.resize (0);

			for (; m; m--)
				_rep.push_back (vector<Element> (n));
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
				_MD.dotprod (y[i], _rep[i], x);
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
			{ _MD.multiDotprod (y, _rep, x); }

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		int rowdim ()
			{ return _rep.size();}
      
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		int coldim () 
			{ return _rep[0].size(); }

		void read (std::istream &file)
		{
			int i;
			pointer p;

			for (i = 0; i < rows; ++i) {
				for (p = _res[i].begin (); p != _res.[i].end (); ++p) {
					file.ignore (1);
					_MD.read (file, *p);
				}
			}
		}
      
		std::ostream &write(std::ostream &os = std::cout)
		{
			int i;
			pointer p;

			for (i = 0; i < _res.size (); ++i) {
				for (p = _res[i].begin (); p != _res[i].end(); ++p) {
					_MD.write (os,*p);
					os << " ";
				}

				os << endl;
			}	
		}

	    private:

		std::vector< std::vector <Element> >  _rep;
		MatrixDomain                         &_MD;
	};
}

#endif
