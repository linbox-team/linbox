/* Copyright (C) 2010 LinBox
 *
 *
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

#ifndef __LINBOX_format_sparse_row_H
#define __LINBOX_format_sparse_row_H

/* sparse-row.h
 * MatrixStreamReader specialization for matrices in the sparse row format:
 * 1st line: "50 60 S" // #rows #cols letter-S
 * Next lines: (# of entries in row) (column index of 1st entry) (value of 1st entry) ...
 * One next line for each row in matrix.
 */

#include <cstdlib>

#if 0
namespace LinBox__FORMAT_SPARSE_ROW_H
	{ static const char* name = "Sparse Row Format";
	  static const char* shortname = "sparserow"; }
#endif

namespace LinBox {

template<class Field>
class SparseRowReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int _base;
	size_t currentRow, colsLeft;

    protected:

	MatrixStreamError initImpl(const char* firstLine) {
		char* restLine;
		int i = 0;

		// Read m
		this->_m = strtoul(firstLine,&restLine,0);
		if( this->_m == 0 && restLine == firstLine )
			return NO_FORMAT;
		i = restLine - firstLine;

		// Read n
		this->_n = strtoul(firstLine+i,&restLine,0);
		if( this->_n == 0 && restLine == firstLine+i )
			return NO_FORMAT;
		i = restLine - firstLine;

		// Read "S"
		while( firstLine[i] && isspace(firstLine[i]) )
			++i;
		if( !firstLine[i] || (firstLine[i] != 'S' &&
		                      firstLine[i] != 's'   ) )
			return NO_FORMAT;

		// Check whitespace for rest of line
		++i;
		while( firstLine[i] && isspace(firstLine[i]) )
			++i;
		if( firstLine[i] ) return BAD_FORMAT;

		this->knowM = this->knowN = true;

		currentRow = (size_t) -1;
		colsLeft = 0;
		return GOOD;
	}

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
	        while( colsLeft == 0 ) {
	          	if( ++currentRow == this->_m ) return END_OF_MATRIX;
			this->ms->readWhiteSpace();
	          	*(this->sin) >> colsLeft;
	          	if( this->sin->eof() ) return END_OF_FILE;
	          	if( !this->sin->good() ) return BAD_FORMAT;
	        }
       
		this->ms->readWhiteSpace();
	        *(this->sin) >> n;
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;
       
		this->ms->readWhiteSpace();
	        this->ms->getField().read(*(this->sin),v);
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;

		n -= _base;
		m = currentRow;
		--colsLeft;

		if(  m >= this->_m ||
		     n >= this->_n ) return BAD_FORMAT;

		return GOOD;
	}

    public:
    	SparseRowReader( int base = 0 ) {
		_base = base;
		currentRow = colsLeft = (size_t) -1;
	}

	const char* getName() const {return "Sparse Row Format"; }//LinBox__FORMAT_SPARSE_ROW_H::name;
	const char* shortName() const
	{ return "sparserow"; }//LinBox__FORMAT_SPARSE_ROW_H::shortname; 

	bool isSparse() const { return true; }
};

}

#endif // __LINBOX_format_sparse_row_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
