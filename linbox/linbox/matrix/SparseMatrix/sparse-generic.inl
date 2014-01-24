/* linbox/matrix/sparse.inl
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Based on sparse-base.h by William J Turner <wjturner@math.ncsu.edu>
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/sparse-base.inl to matrix/sparse.inl
 * ------------------------------------
 * 2002-11-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 *   - Renamed ColOfRowsIterator to RowIterator
 *   - Named template argument _Row rather than Row; add a typedef to Row
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

#ifndef __LINBOX_matrix_sparse_INL
#define __LINBOX_matrix_sparse_INL

#include "linbox-config.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>

#include "linbox/matrix/sparse.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"

namespace LinBox { namespace Protected {

	template <class Field, class Row, class Trait>
	std::istream &SparseMatrixReadWriteHelper<Field, Row, Trait> ::readTurner (SparseMatrix<Field, Row> &A, std::istream &is
										   , char *buf)
	{
		size_t i, j;

		A._matA.clear ();
		A._matA.resize (A._m);

		do {
			std::istringstream str (buf);

			str >> i;

			if (i == (size_t) -1) break; // return also if row index is -1
			str >> j;
			A.field().read (str, A.refEntry (i, j));

			is.getline (buf, 80);
		} while (is);

		return is;
	}

	template <class Field, class Row, class Trait>
	std::istream &SparseMatrixReadWriteHelper<Field, Row, Trait> ::readGuillaume (SparseMatrix<Field, Row> &A, std::istream &is
										      , char *buf)
	{
		size_t i = 0, j = 0 ;

		std::istringstream str (buf);
		str >> A._m >> A._n;

		A._matA.clear ();
		A._matA.resize (A._m);//cerr<<A.coldim()<<" "<<A.rowdim()<<endl;

		Element x;
		while (is >> i) {
			if (i == 0 || i == (size_t) -1) {is >> j; A.field().read(is, x); break;}
			is >> j;
			if (i > A._m || j > A._n)
				throw Exceptions::InvalidMatrixInput ();
			A.field().read (is, x);
			if (! A.field().isZero(x)) A.setEntry (i - 1, j - 1, x);
		}

		return is;

	}

	template <class Field, class Row, class Trait>
	std::istream &SparseMatrixReadWriteHelper<Field, Row, Trait> ::readMatlab (SparseMatrix<Field, Row> &A, std::istream &is
										   , char *buf)
	{
		size_t i = 0, j = 0;
		char c;
		Element a_ij;

		while (1) {
			do is >> c; while (is && !std::isdigit (c));
			if (!is) break;

			is.putback (c);

			A.field().read (is, a_ij);
			A.setEntry (i, j++, a_ij);

			do is >> c; while (is && c != ',' && c != ';' && c != ']');
			if (!is) break;;

			if (c == ';') {
				++i;
				j = 0;
			}
			else if (c == ']') break;
		}

		return is;
	}

	template <class Field, class Row, class Trait>
	std::istream &SparseMatrixReadWriteHelper<Field, Row, Trait> ::readPretty (SparseMatrix<Field, Row> &A, std::istream &is
										   , char *buf)
	{
		size_t i;
		Element a_ij;

		A._m = 0;
		A._matA.clear ();

		i = 0;

		do {
			char c;
			size_t j;
			++(A._m);
			A._matA.push_back (Row ());

			std::istringstream str (buf);

			do {str >> c; std::cout << "x" << c << "x" << std::endl;} while (isspace (c));
			if (c != '[')
				throw Exceptions::InvalidMatrixInput ();

			j = 0;

			while (str) {
				do str >> c; while (isspace (c));
				if (!str || c == ']') break;
				A.field().read (str, a_ij);

				++j;
				if (j > A._n)
					++(A._n);

				if (!A.field().isZero(a_ij))
					A.setEntry (i, j, a_ij);
			}

			is.getline (buf, 80);

			++i;
		} while (is);

			return is;

	}

	template <class Field, class Row, class Trait>
	std::istream &SparseMatrixReadWriteHelper<Field, Row, Trait> ::readMagmaCpt (SparseMatrix<Field, Row> &A, std::istream &is
										     , char *buf)
	{
		size_t i, j;
		Element a_ij;
		char c;
		const char matrixstart = '[', matrixend = ']';
		const char rowstart = '[', rowend = ']';
		const char pairstart = '[', pairend = ']';

		A._m = A._n = 0;
		A._matA.clear ();

		do {is.get(c);} while (c != matrixstart ); // find matrix start
		i = 0;
		while (true)
		{
			do {is.get(c);} while (c != matrixend && c != rowstart);
			if (c == matrixend) return is;
			else
			{
				++(A._m);
				A._matA.push_back (Row ());
				//processrow(i)
				while (true)
				{
					do {is.get(c);} while (c != pairstart && c != rowend );
					if (c == rowend) break;
					else
					{  //processpair( j v for row i);
						is >> j;
						if (j > A._n) A._n = j;
						do {is.get(c);} while (!std::isdigit(c) && c != '-' && c != '+');
						is.unget();
						A.field().read(is, a_ij);
						if (!A.field().isZero(a_ij)) A.setEntry (i, j-1, a_ij);
						do {is.get(c);} while (c != pairend);
					}
				}
				++i;
			}
		}
		//return is; //BB: unreachable
	}

	template <class Field, class Row, class Trait>
	std::istream &SparseMatrixReadWriteHelper<Field, Row, Trait> ::read (SparseMatrix<Field, Row> &A, std::istream &is
									     , LINBOX_enum(Tag::FileFormat) format)
	{
		char buf[80];
		buf[0]=0;

		switch (format) {
		case Tag::FileFormat::Detect: {
					    char c;
					    is.getline (buf, 80);
					    std::istringstream str (buf);
					    do str >> c; while (isspace (c));

					    if (c == '[') {
						    if (strchr (buf, ';') != NULL)
							    readMatlab (A, is, buf);
						    else
							    readPretty (A, is, buf);
					    }
					    else if (std::isdigit (c)) {
						    do str >> c; while (str && (isspace (c) || std::isdigit (c)));

						    if (c == 'M')
							    return readGuillaume (A, is, buf);
						    else
							    return readTurner (A, is, buf);
					    }
					    else
						    throw Exceptions::InvalidMatrixInput ();
					    break;
				    }

		case Tag::FileFormat::Turner:
				    return readTurner (A, is, buf);
		case Tag::FileFormat::Guillaume:
				    return readGuillaume (A, is, buf);
		case Tag::FileFormat::Matlab:
				    return readMatlab (A, is, buf);
		case Tag::FileFormat::Pretty:
				    return readPretty (A, is, buf);
		case Tag::FileFormat::MagmaCpt:
				    return readMagmaCpt (A, is, buf);
		default:
				    throw Exceptions::InvalidMatrixInput();
		}

		return is;
	}

	template <class Field, class Row, class Trait>
	std::ostream &SparseMatrixWriteHelper<Field, Row, Trait> ::write (const SparseMatrix<Field, Row> &A, std::ostream &os
									  , LINBOX_enum(Tag::FileFormat) format)
	{
		typename SparseMatrix<Field, Row>::Rep::const_iterator i;
		typename Row::const_iterator j;
		size_t i_idx, j_idx;
		//	int col_width;
		integer c;
		bool firstrow;

		// Avoid massive unneeded overhead in the case that this
		// printing is disabled
		if (not os)
			return os;

		switch (format) {
		case Tag::FileFormat::Detect:
			throw PreconditionFailed (__func__, __LINE__, "format != Tag::FileFormat::Detect");
			//break;//BB: unreachable

		case Tag::FileFormat::Turner:
			// The i j v triples, with zero based indices.
			for (i = A._matA.begin (), i_idx = 0; i != A._matA.end (); ++i, ++i_idx) {
				for (j = i->begin (), j_idx = 0; j != i->end (); ++j, ++j_idx) {
					os << i_idx << ' ' << j->first << ' ';
					A.field().write (os, j->second);
					os << std::endl;
				}
			}
			break;

		case Tag::FileFormat::OneBased:
			// The i j v triples, with zero based indices.
			for (i = A._matA.begin (), i_idx = 0; i != A._matA.end (); ++i, ++i_idx) {
				for (j = i->begin (), j_idx = 0; j != i->end (); ++j, ++j_idx) {
					os << i_idx + 1 << ' ' << j->first + 1 << ' ';
					A.field().write (os, j->second);
					os << std::endl;
				}
			}
			break;

		case Tag::FileFormat::Guillaume:
			// row col 'M' header line followed by the i j v triples, one based,
			// followed by 0 0 0.
			os << A._m << ' ' << A._n << " M" << std::endl;

			for (i = A._matA.begin (), i_idx = 0; i != A._matA.end (); ++i, ++i_idx) {
				for (j = i->begin (), j_idx = 0; j != i->end (); ++j, ++j_idx) {
					os << i_idx + 1 << ' ' << j->first + 1 << ' ';
					A.field().write (os, j->second);
					os << std::endl;
				}
			}

			os << "0 0 0" << std::endl;

			break;

		case Tag::FileFormat::Matlab:

			os << "[";

			for (i = A._matA.begin (), i_idx = 0; i != A._matA.end (); ++i, ++i_idx) {
				j = i->begin ();

				for (j_idx = 0; j_idx < A._n; ++j_idx) {
					if (j == i->end () || j_idx != j->first)
						A.field().write (os, A.field().zero);
					else {
						A.field().write (os, j->second);
						++j;
					}

					if (j_idx < A._n - 1)
						os << ", ";
				}

				os << "; ";
			}

			os << "]";

			break;

		case Tag::FileFormat::Maple:

			os << "[";
			firstrow=true;

			for (i = A._matA.begin (), i_idx = 0; i != A._matA.end (); ++i, ++i_idx) {
				if (firstrow) {
					os << "[";
					firstrow =false;
				}
				else
					os << ", [";

				j = i->begin ();

				for (j_idx = 0; j_idx < A._n; ++j_idx) {
					if (j == i->end () || j_idx != j->first)
						A.field().write (os, A.field().zero);
					else {
						A.field().write (os, j->second);
						++j;
					}

					if (j_idx < A._n - 1)
						os << ", ";
				}

				os << " ]";
			}

			os << "]";

			break;

		case Tag::FileFormat::Pretty:
			//A.field().characteristic (c);
			//col_width = (int) ceil (log ((double) c) / M_LN10);

			for (i = A._matA.begin (), i_idx = 0; i != A._matA.end (); ++i, ++i_idx) {
				os << "  [ ";

				j = i->begin ();

				for (j_idx = 0; j_idx < A._n; ++j_idx) {
					//os.width (col_width);

					if (j == i->end () || j_idx != j->first)
						A.field().write (os, A.field().zero);
					else {
						A.field().write (os, j->second);
						++j;
					}

					os << ' ';
				}

				os << ']' << std::endl;
			}

			break;

		case Tag::FileFormat::MagmaCpt:
			os << "sparse matrix written in MagmaCpt form is not implemented" << std::endl;
			break;

		default:
			os << "sparse matrix written in Tag::FileFormat::" << (int)format << " is not implemented" << std::endl;
		}

		return os;
	}


	template <class Field, class Row, class Tag>
	SparseMatrix<Field,Row,Tag> ::SparseMatrix( MatrixStream<Field>& ms ) :
		_matA(0), _m(0), _n(0)
	{
		Element val;
		size_t i, j;
		while( ms.nextTriple(i,j,val) ) {
			if( i >= _m ) {
				_m = i + 1;
				_matA.resize( _m );
			}
			if( j >= _n ) _n = j + 1;
			setEntry(i,j,val);
		}
		if( ms.getError() > END_OF_MATRIX )
			throw ms.reportError(__func__,__LINE__);
		if( !ms.getDimensions( i, _n ) )
			throw ms.reportError(__func__,__LINE__);
		if( i > _m ) {
			_m = i;
			_matA.resize(_m);
		}
	}



#if 0 /*  weird */
	template <class Field,  class _Row, class BRow>
	SparseMatrix<Field,_Row> *SparseMatrixFactory<Field, _Row, BRow>::makeBlackbox (const Field &F)
	{
		SparseMatrix<Field, _Row> *A = new SparseMatrix<Field, _Row> (F, rowdim (), coldim ());

		typename SparseMatrix<Field, BRow>::ConstIterator i;
		typename SparseMatrix<Field, BRow>::ConstIndexedIterator j;

		for (i = _matA.Begin (), j = _matA.IndexedBegin (); i != _matA.End (); ++i, ++j)
			A.field().init (A->refEntry (j.rowIndex (), j.colIndex ()), *i);

		return A;
	}
#endif

} // namespace LinBox
} // namespace Protected

#endif // __LINBOX_matrix_sparse_INL


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
