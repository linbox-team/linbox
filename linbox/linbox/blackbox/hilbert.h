/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/blackbox/hilbert.h
 * Copyright (C) 2006 John P. May, B. David Saunders
 *
 * Written by John P. May <jpmay@cis.udel.edu>,
 *            B. David Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * The original hilbert.h, providing one of the first blackbox examples, but not using the JIT feature,
 * was written by Will Turner.
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
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file blackbox/hilbert.h
 * @ingroup blackbox
 * @brief NO DESC
 */

#ifndef __LINBOX_hilbert_H
#define __LINBOX_hilbert_H

#include <vector>
#include "linbox/blackbox/jit-matrix.h"


namespace LinBox
{

	/// The object needed to build a Hilbert matrix as a JIT matrix
	template<typename _Field>
	class Hilbert_JIT_Entry {

	public:
		typedef _Field Field;
		typedef typename _Field::Element Element;

		/// set up vector of 1/(i+1)
		Hilbert_JIT_Entry(Field& F, size_t m, size_t n);

		/// return 1/(i+j+2), zero based indexing.
		Element& operator()(Element &entry, size_t i, size_t j) const
		{
			return entry = _vecH[i+j+1];
		}

	private:
		std::vector<Element> _vecH;

	}; // Hilbert_JIT_Entry

	/// constructor
	template<typename _Field>
	Hilbert_JIT_Entry<_Field>::Hilbert_JIT_Entry(_Field& F, size_t m, size_t n) {

		Element temp, one;
		F.init(one, 1);
		F.init(temp, 0);

		_vecH = std::vector<Element>(m+n, temp);

		typename std::vector<Element>::iterator iter;

		// the ith entry of _vecH = 1/(i+1)
		for (iter=_vecH.begin(); iter != _vecH.end(); iter++) {
			F.addin(temp, one);
			F.inv(*iter, temp);
		}

	}//constructor


	/** \brief Example of a blackbox that is space efficient, though not time efficient.
	 *
	 *  \ingroup blackbox
	 *
	 *  Blackbox for the matrix whose i,j entry is 1/(i+j), i in 1..n, j in 1..n.
	 *
	 */
	template<typename _Field>
	class Hilbert : public JIT_Matrix<_Field, Hilbert_JIT_Entry<_Field> > {

	public:
		/** Constructor from field and size.
		 * @param F the field.
		 * @param n the size : size_t integer number of rows and columns of matrix.
		 */
		Hilbert(_Field& F, size_t n) :
			JIT_Matrix<_Field, Hilbert_JIT_Entry<_Field> >(F, n, n, Hilbert_JIT_Entry<_Field>(F, n, n))
		{};

	};



}//LinBox Namespace

#endif //__LINBOX_hilbert_H

