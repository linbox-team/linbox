/* linbox/solutions/getentry.h
 * Copyright(C) LinBox
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *  -bds
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

#ifndef __LINBOX_getentry_H
#define __LINBOX_getentry_H

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/solution-tags.h"

namespace LinBox
{

	/** \brief
	 * Getting the i,j entry of the blackbox.
	 */
	template <class BB>
	typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j);

	// Some specializations

	// Compose< Diagonal, BB > specialization
	template <class Field, class Trait, class BB>
	typename Field::Element& getEntry(typename Field::Element& x, const Compose<Diagonal<Field, Trait>, BB>& A, const size_t i, const size_t j);

	// Compose< BB, Diagonal > specialization
	template <class BB, class Field, class Trait>
	typename Field::Element& getEntry(typename Field::Element& x, const Compose<BB, Diagonal<Field, Trait> >& A, const size_t i, const size_t j);

	// Compose< Diagonal, Diagonal > specialization
	template <class Field, class T1, class T2>
	typename Field::Element& getEntry(typename Field::Element& x, const Compose<Diagonal<Field,T1>, Diagonal<Field, T2> >& A, const size_t i, const size_t j);

	/// To ignore methods
	template <class BB, class Method>
	typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, Method & m);

	/** GetEntryCategory is specialized for BB classes that offer a local getEntry

	   This includes SparseMatrix, SparseMatrixBase, Diagonal, ScalarMatrix.
	   It could and should include many more.
	 */
	template<class BB> struct GetEntryCategory;

/************************** internal forms **********************/
	// For the general case apply() will be used.
	template <class BB>
	typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, SolutionTags::Generic t);

	// Some BBs have a local getEntry method
	template <class BB>
	typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, SolutionTags::Local t );

} // LinBox

#include "linbox/solutions/getentry.inl"

#endif // __LINBOX_getentry_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
