/* Copyright (C) 2005 LinBox
 * Written by  Dan Roche
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_format_matrix_market_H
#define __LINBOX_format_matrix_market_H

/* matrix-market-array.h
 * MatrixStreamReader specialization for matrices in the MatrixMarket coordinate
 * format.
 */

#include <string>
#include <sstream>
#include "linbox/util/matrix-stream.h"
#include "linbox/integer.h"

#if 0
namespace LinBox__FORMAT_MATRIX_MARKET_H
{
	static const char* name = "Matrix Market Format";
	static const char* shortname = "mm";
}
#endif

namespace LinBox
{

#if 0
#ifndef __LINBOX_integer_H
	class integer;
#endif
#endif

#if 1 /* Unused */
	inline bool equalCaseInsensitive(const std::string &s1, const char* s2)
	{
		int len = int(s1.size());
		int counter = 0;
		while( counter < len && s2[counter] != '\0' &&
		       toupper(s1[(size_t)counter]) == toupper(s2[counter]) ) {
			++counter;
		}
		return( counter == len && s2[counter] == '\0' );
	}
#endif

	template<class Field>
	class MatrixMarketReader : public MatrixStreamReader<Field> {
	public:
		typedef typename MatrixStreamReader<Field>::Element Element;
	private:
		int entriesLeft;
		size_t currentCol, currentRow;
		bool array;
		bool pattern;
		bool symmetric;

		MatrixStreamError readHeader()
		{
			//Skip comments
			this->ms->readWhiteSpace();
			while( !this->sin->eof() && this->sin->peek() == '%' ) {
				char c;
				while( this->sin->get(c) ) {
					if( c == '\n' || c == '\r' ) {
						this->sin->putback(c);
						break;
					}
				}
				this->ms->readWhiteSpace();
			}

			*(this->sin) >> this->_m;
			this->ms->readWhiteSpace();
			*(this->sin) >> this->_n;
			this->ms->readWhiteSpace();
			if( !array ) {
				*(this->sin) >> entriesLeft;
				this->ms->readWhiteSpace();
			}

                        if (entriesLeft != 0) {
                                if( this->sin->eof()) return END_OF_FILE;
                                if( !this->sin->good()) return BAD_FORMAT;
                        }

			this->knowM = this->knowN = true;
			currentCol = currentRow = 1;

			if( symmetric && (this->_m != this->_n) ) return BAD_FORMAT;
			if( this->_m < 1 || this->_n < 1 ) return BAD_FORMAT;
			if( !array && (entriesLeft < 0 ||
				       (size_t)entriesLeft > this->_m*this->_n ) )
				return BAD_FORMAT;

			return GOOD;
		}

	protected:
		MatrixStreamError nextTripleImpl( size_t& i, size_t& j, Element& v )
		{
			if( currentCol == 0 && currentRow == 0 ) {
				MatrixStreamError mse = readHeader();
				if( mse != GOOD ) return mse;
			}

			if( array ) {
				if( currentCol == this->_n+1 ) return END_OF_MATRIX;
				j = currentCol;
				i = currentRow;
				if( ++currentRow == this->_m+1 ) {
					++currentCol;
					currentRow = (symmetric ? currentCol : 1);
				}
			}
			else {
				if( --entriesLeft < 0 ) return END_OF_MATRIX;

				this->ms->readWhiteSpace();
				*(this->sin) >> i;
				if( this->sin->eof() ) return END_OF_FILE;
				if( !this->sin->good() ) return BAD_FORMAT;

				this->ms->readWhiteSpace();
				*(this->sin) >> j;
				if( this->sin->eof() ) return END_OF_FILE;
				if( !this->sin->good() ) return BAD_FORMAT;
			}

			if( pattern ) this->ms->getField().assign(v,this->ms->getField().one);
			else {
				this->ms->readWhiteSpace();
				this->ms->getField().read(*(this->sin),v);
				if( this->sin->eof() ) return END_OF_FILE;
				if( !this->sin->good() ) return BAD_FORMAT;
			}

			--i;
			--j;
			if(  i >= this->_m || j >= this->_n )
				return BAD_FORMAT;
			if( symmetric && (i != j) ) this->saveTriple(j,i,v);

			return GOOD;
		}

#if 1 /* Unused */
		MatrixStreamError initImpl( const char* firstLine )
		{
			std::string st(firstLine);
			std::stringstream stin(st);

			//bb: cppcheck complained.
			// if( stin.get() != '%' || stin.get() != '%' ) return NO_FORMAT;
			if( stin.get() != '%')
				return NO_FORMAT;
			else if (stin.get() != '%' )
				return NO_FORMAT;

			if( !stin.good() ) return NO_FORMAT;

			std::string s;
			stin >> s;
			if( !stin.good() ) return NO_FORMAT;
			if( !equalCaseInsensitive(s,"MatrixMarket") ) return NO_FORMAT;

			stin >> s;
			if( !stin.good() ) return BAD_FORMAT;
			if( !equalCaseInsensitive(s,"matrix") ) return BAD_FORMAT;

			stin >> s;
			if( !stin.good() ) return BAD_FORMAT;
			if( equalCaseInsensitive(s,"array") ) array = true;
			else if( equalCaseInsensitive(s,"coordinate") ) array = false;
			else return BAD_FORMAT;

			stin >> s;
			if( !stin.good() ) return BAD_FORMAT;
			pattern = equalCaseInsensitive(s,"pattern");

			stin >> s;
			if( !stin.eof() && !stin.good() ) return BAD_FORMAT;
			if( equalCaseInsensitive(s,"symmetric") ) symmetric = true;
			else if( equalCaseInsensitive(s,"general") ) symmetric = false;
			else return BAD_FORMAT;

			stin >> s;
			if( !stin.eof() ) return BAD_FORMAT;

			if( array && pattern ) return BAD_FORMAT;

			currentRow = currentCol = 0;

			return GOOD;
		}
#endif

	public:
		MatrixMarketReader()
		{
			entriesLeft = -1;
			currentCol = currentRow = 0;
			//bb : why not intialised ?
			array     = false ;
			pattern   = false ;
			symmetric = false ;
		}

		bool isSparse() const
		{
			return !array;
		}

		const char* getName() const
		{
			return "Matrix Market Format";
		}// LinBox__FORMAT_MATRIX_MARKET_H::name;

		const char* shortName() const
		{
			return "mm"; // LinBox__FORMAT_MATRIX_MARKET_H::shortname;
		}
	};

}

#endif //__LINBOX_format_matrix_market_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
