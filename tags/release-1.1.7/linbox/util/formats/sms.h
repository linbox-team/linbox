/* Copyright (C) 2005 LinBox
 * Written by Dan Roche
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

/*! @file util/formats/sms.h
 * @brief MatrixStreamReader for sms matrix format
 * 
 * 1st line: "50 60 X" // #rows #cols a letter (e.g M, I, R, P ...)
 * Subsequent lines: i j v // row index, col index, value
 * last line: 0 0 0
 */

#ifndef __LINBOX_sms_H
#define __LINBOX_sms_H

#include <cstdlib>

/*
namespace LinBox__SMS_H
	{ static const char* name = "SMS Sparse Integer Matrix Format";
	  static const char* shortname = "sms"; }
	  */

namespace LinBox {

template<class Field>
class SMSReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int _base;

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

		// Read "M" or "R" or "P" or "I"
		while( firstLine[i] && isspace(firstLine[i]) )
			++i;
		if( !firstLine[i] || (firstLine[i] != 'M' &&
		                      firstLine[i] != 'm' &&
		                      firstLine[i] != 'I' &&
		                      firstLine[i] != 'i' &&
		                      firstLine[i] != 'R' &&
		                      firstLine[i] != 'r' &&
		                      firstLine[i] != 'P' &&
		                      firstLine[i] != 'p'   ) )
			return NO_FORMAT;

		// Check whitespace for rest of line
		++i;
		while( firstLine[i] && isspace(firstLine[i]) )
			++i;
		if( firstLine[i] ) return BAD_FORMAT;

		this->knowM = this->knowN = true;

		return GOOD;
	}

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
		this->ms->readWhiteSpace();
	        *(this->sin) >> m;
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;
       
		this->ms->readWhiteSpace();
	        *(this->sin) >> n;
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;
       
		this->ms->readWhiteSpace();
		if( this->sin->eof() ) return END_OF_FILE;
	        this->ms->getField().read(*(this->sin),v);
	        if( this->sin->eof() ) this->atEnd = true;
	        else if( !this->sin->good() ) return BAD_FORMAT;

		if( m == 0 && n == 0 ) return END_OF_MATRIX;

		m -= _base;
		n -= _base;

		if( m >= this->_m ||
		    n >= this->_n ) return BAD_FORMAT;

		return GOOD;
	}

    public:
    	SMSReader( int base = 1 ) {
		_base = base;
	}

	const char* getName() const {return "SMS Sparse Integer Matrix Format"; }//LinBox__SMS_H::name;}
	const char* shortName() const
	{ return "sms"; }//LinBox__SMS_H::shortname; }

	bool isSparse() const { return true; }
};

}

#endif // __LINBOX_sms_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
