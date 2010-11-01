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


#ifndef __LINBOX_format_dense_H
#define __LINBOX_format_dense_H

/* dense.h
 * MatrixStreamReader specialization for matrices in the generic dense format:
 * 1st line: #rows #cols
 * Subsequent lines: every entry, row by row.
 */

#if 0
namespace LinBox__FORMAT_DENSE_H
	{ static const char* name = "Generic Dense Format";
	  static const char* shortname = "dense"; }
#endif

namespace LinBox 
{

template<class Field>
class DenseReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	size_t currentRow, currentCol;
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
		if( this->_n == 0 && restLine == (firstLine+i) )
			return NO_FORMAT;
		i = restLine - firstLine;

		// Check whitespace for rest of line
		++i;
		while( firstLine[i] && isspace(firstLine[i]) )
			++i;
		if( firstLine[i] ) return BAD_FORMAT;

		this->knowM = this->knowN = true;

		currentRow = currentCol = 0;
		return GOOD;
	}

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
	    if( currentRow == this->_m ) return END_OF_MATRIX;
	    m = currentRow;
	    n = currentCol;

	    this->ms->readWhiteSpace();
	    this->ms->getField().read(*(this->sin),v);
	    if( this->sin->eof() ) return END_OF_FILE;
	    if( !this->sin->good() ) return BAD_FORMAT;
	    
	    if( ++currentCol == this->_n ) {
	    	++currentRow;
		currentCol = 0;
	    }
	    return GOOD;
	}

    public:
    	DenseReader() {
		currentRow = currentCol = (size_t) -1;
	}

	bool isSparse() const { return false; }
    	
	const char* getName() const { return "Generic Dense Format"; }//LinBox__FORMAT_DENSE_H::name; }
	
	const char* shortName() const 
		{ return "dense"; }//LinBox__FORMAT_DENSE_H::shortname; 

};

}

#endif // __LINBOX_format_dense_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
