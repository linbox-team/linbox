/* Copyright (C) 2005 LinBox
 * Written by  Dan Roche
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

#ifndef __LINBOX_format_maple_H
#define __LINBOX_format_maple_H

/* matrix-market-array.h
 * MatrixStreamReader for Maple-style (text) matrices.
 * The general format is Matrix(r,c,<init> [,...]),
 * where <init> is either a row-column-value specification, i.e.
 * {(i,j)=v,(i,j)=v,...} or a dense specification, i.e.
 * [[v11,v12,...],[v21,v22,...],...].
 * In the second case, the rows/columns specifications are optional,
 * and even the outer Matrix(...) part is. In any case, either the word
 * "Matrix" or two opening square brackets "[[" must be present on the
 * first line so the format can be identified.
 */

#include <string>
#include <sstream>
#include <linbox/util/matrix-stream.h>

#if 0
namespace LinBox__FORMAT_MAPLE_H
	{ static const char* name = "Maple Text Format";
	  static const char* shortname = "maple"; }
#endif

namespace LinBox 
{

template<class Field>
class MapleReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
	/* These are used to keep track of where we are in the data for
	 * array-based data. They are also used to keep track of where we
	 * are in reading the frontmatter of the matrix.
	 * If currentCol is 0, then we have not reached the matrix data yet.
	 * Then currentRow can be 0,2,3,4,5, depending on how many of the
	 * tokens ( r , c , have been read.
	 */
	size_t currentCol, currentRow;
	bool openParen;
	bool array;
	std::stringstream *stin;

	MatrixStreamError processCandidate( const char* cand ) {
		if( strlen(cand) < 6 ||
		    tolower(cand[1]) != 'a' ||
		    tolower(cand[2]) != 't' ||
		    tolower(cand[3]) != 'r' ||
		    tolower(cand[4]) != 'i' ||
		    tolower(cand[5]) != 'x' )
			return NO_FORMAT;
		
		int i = 6;
		while( cand[i] && isspace(cand[i]) ) ++i;
		
		if( !cand[i] ) {
			openParen = true;
			return GOOD;
		}
		
		if( cand[i] != '(' ) return BAD_FORMAT;

		openParen = true;
		currentCol = 1;

		++i;
		while( cand[i] && isspace(cand[i]) ) ++i;
		if( !cand[i] ) return GOOD;

		char* pastNum;
		if( std::isdigit(cand[i]) ) {
			this->_m = strtoul( cand+i, &pastNum, 0 );
			if( this->_m == 0 && pastNum == cand+i )
				return BAD_FORMAT;
			this->knowM = true;
			currentCol = 2;
			i = pastNum - cand;
		}
		else {
			currentCol = i;
			return GOOD;
		}

		while( cand[i] && isspace(cand[i]) ) ++i;
		if( !cand[i] ) return GOOD;
		if( cand[i] != ',' ) return BAD_FORMAT;
		currentCol = 3;

		++i;
		while( cand[i] && isspace(cand[i]) ) ++i;
		if( !cand[i] ) return GOOD;

		if( std::isdigit(cand[i]) ) {
			this->_n = strtoul( cand+i, &pastNum, 0 );
			if( this->_n == 0 && pastNum == cand+i )
				return BAD_FORMAT;
			this->knowN = true;
			currentCol = 4;
			i = pastNum - cand;
		}
		else {
			currentCol = i;
			return GOOD;
		}

		while( cand[i] && isspace(cand[i]) ) ++i;
		if( !cand[i] ) return GOOD;
		if( cand[i] != ',' ) return BAD_FORMAT;
		currentCol = i;
		return GOOD;
	}

	MatrixStreamError readUntil(char end) {
		if( stin ) {
			while( !stin->eof() && stin->get() != end )
			if( stin->eof() ) {
				delete stin;
				stin = NULL;
			}
			else return GOOD;
		}

		this->ms->readWhiteSpace();
		while( this->sin->get() != end )
			this->ms->readWhiteSpace();
		if( this->sin->eof() ) return END_OF_FILE;
		return GOOD;
	}
	
	MatrixStreamError readWhite() {
		if( stin ) {
			int peekVal = stin->peek();
			while( stin->good() && isspace(peekVal) ) {
				stin->get();
				peekVal = stin->peek();
			}
			if( !stin->good() || peekVal < 0 ) {
				delete stin;
				stin = NULL;
			}
		}
		if( !stin ) { // Note this is NOT equivalent to "else"
			this->ms->readWhiteSpace();
			if( this->sin->eof() ) return END_OF_FILE;
		}
		return GOOD;
	}

	MatrixStreamError readCharacter(char &c) {
		MatrixStreamError mse = readWhite();
		if( mse > GOOD ) return mse;
		if( stin ) stin->get(c);
		else this->sin->get(c);
		return GOOD;
	}

	MatrixStreamError readElement(Element& ele) {
		MatrixStreamError mse = readWhite();
		if( mse > GOOD ) return mse;
		if( stin ) {
			this->ms->getField().read(*stin,ele);
			if( stin->eof() ) {
				delete stin;
				stin = NULL;
			}
			else if( !stin->good() )
				return BAD_FORMAT;
		}
		else {
			this->ms->getField().read(*(this->sin),ele);
			if( !this->sin->eof() && !this->sin->good() )
				return BAD_FORMAT;
		}
		return GOOD;
	}
	
	MatrixStreamError readNumber(size_t& num) {
		MatrixStreamError mse = readWhite();
		if( mse > GOOD ) return mse;
		if( stin ) {
			(*stin) >> num;
			if( stin->eof() ) {
				delete stin;
				stin = NULL;
			}
			else if( !stin->good() )
				return BAD_FORMAT;
		}
		else {
			*(this->sin) >> num;
			if( !this->sin->eof() && !this->sin->good() )
				return BAD_FORMAT;
		}
		return GOOD;
	}

    protected:
    	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
		if( currentRow == 0 ) {
			while( currentCol < 7 ) {
				this->ms->readWhiteSpace();
				if( this->sin->eof() )
					return END_OF_FILE;

				char ch;
				switch( currentCol ) {
				    case 0:
					if( this->sin->get() != '(' )
						return BAD_FORMAT;
					break;
				    case 1:
				    	if( std::isdigit(this->sin->peek()) ) {
						*(this->sin) >> this->_m;
						if( this->sin->eof() )
							return END_OF_FILE;
						if( !this->sin->good() )
							return BAD_FORMAT;
						this->knowM = true;
					}
					else currentCol = 4;
					break;
				    case 2:
				    case 4:
				    	if( this->sin->get() != ',' )
						return BAD_FORMAT;
					break;
				    case 3:
				    	if( std::isdigit(this->sin->peek()) ) {
						*(this->sin) >> this->_n;
						if( this->sin->eof() )
							return END_OF_FILE;
						if( !this->sin->good() )
							return BAD_FORMAT;
						this->knowN = true;
					}
					else currentCol = 4;
					break;
				    case 5:
				    	ch = this->sin->get();
					if( ch == '{' ) {
						array = false;
						currentCol = 6;
					}
					else if( ch == '[' )
						array = true;
					else return BAD_FORMAT;
					break;
				    case 6:
				    	if( this->sin->get() != '[' )
						return BAD_FORMAT;
					break;
				}
				++currentCol;
			}

			currentRow = currentCol = 1;
		}

		if( array ) {
			MatrixStreamError mse = readElement(v);
			if( mse > GOOD ) return mse;
			m = currentRow - 1;
			n = currentCol - 1;
			
			char c;
			mse = readCharacter(c);
			if( mse > GOOD ) return mse;
			if( c == ',' ) {
				++currentCol;
				if( this->knowN && currentCol > this->_n )
					return BAD_FORMAT;
			}
			else if( c == ']' ) {
				if( !this->knowN ) {
					this->knowN = true;
					this->_n = currentCol;
				}
				currentCol = 1;
				mse = readCharacter(c);
				if( mse > GOOD ) return mse;
				if( c == ',' ) {
					++currentRow;
					if( this->knowM && 
					    currentRow > this->_m )
					    	return BAD_FORMAT;
					mse = readCharacter(c);
					if( mse > GOOD ) return mse;
					if( c != '[' ) return BAD_FORMAT;
				}
				else if( c == ']' ) {
					if( !this->knowM ) {
						this->knowM = true;
						this->_m = currentRow;
					}
					if( openParen ) {
						mse = readUntil(')');
						if( mse > GOOD ) return mse;
					}
					this->atEnd = true;
				}
				else return BAD_FORMAT;
			}
		}
		else {
			char c;
			MatrixStreamError mse = readCharacter(c);
			if( mse > GOOD ) return mse;
			if( c != '(' ) return BAD_FORMAT;
			mse = readNumber(m);
			if( mse > GOOD ) return mse;
			mse = readCharacter(c);
			if( mse > GOOD ) return mse;
			if( c != ',' ) return BAD_FORMAT;
			mse = readNumber(n);
			if( mse > GOOD ) return mse;
			mse = readCharacter(c);
			if( mse > GOOD ) return mse;
			if( c != ')' ) return BAD_FORMAT;
			mse = readCharacter(c);
			if( mse > GOOD ) return mse;
			if( c != '=' ) return BAD_FORMAT;
			mse = readElement(v);
			if( mse > GOOD ) return mse;
			mse = readCharacter(c);
			if( mse > GOOD ) return mse;
			if( c == '}' ) {
				if( openParen ) {
					mse = readUntil(')');
					if( mse > GOOD ) return mse;
				}
				this->atEnd = true;
			}
			--m;
			--n;
		}
		return GOOD;
	}

	MatrixStreamError initImpl( const char* firstLine ) {
		// First look for the word "Matrix"
		const char* candidate = strpbrk(firstLine,"mM");
		MatrixStreamError procErr;
		while( candidate ) {
			procErr = processCandidate(candidate);
			if( procErr < NO_FORMAT ) {
				if( procErr > GOOD ) return procErr;
				else break;
			}
			candidate = strpbrk(candidate+1,"mM");
		}

		bool lineend = false;
		int i = 0;

		if( candidate ) {
			lineend = currentCol <= 5;
			i = currentCol;
			if( !lineend ) currentCol = 5;
			while( !lineend && currentCol < 7 ) {
				while( candidate[i] && isspace(candidate[i]) )
					++i;
				switch(candidate[i]) {
				    case '\0':
				    	lineend = true;
					break;
				    case '[':
				    	array = true;
					++currentCol;
					break;
				    case '{':
				    	if( currentCol == 5 ) {
						array = false;
						currentCol = 7;
						break;
					}
				    default:
				    	return BAD_FORMAT;
				}
				++i;
			}
		}
		else { // Look for [[
			candidate = strchr(firstLine,'[');
			while( candidate ) {
				i = 1;
				while( candidate[i] && isspace(candidate[i]) )
					++i;
				if( candidate[i] == '[' ) {
					openParen = false;
					currentCol = 7;
					++i;
					break;
				}
				candidate = strchr(candidate+1,'[');
			}
		}

		if( !candidate ) return NO_FORMAT;
		else if( currentCol < 7 ) return GOOD;

		currentRow = currentCol = 1;
		if( this->knowM && !this->knowN ) {
			this->knowN = true;
			this->_n = this->_m;
		}
		if( !array && !this->knowM ) return BAD_FORMAT;

		// Now candidate+i is the beginning of the actual data.
		while( candidate[i] && isspace(candidate[i]) ) ++i;
		if( !candidate[i] ) return GOOD;

		std::string st(candidate+i);
		stin = new std::stringstream(st);
		
		return GOOD;
	}

    public:
    	MapleReader() {
		currentCol = currentRow = 0;
		stin = NULL;
	}

	~MapleReader() {
		if( stin ) delete stin;
	}

	bool isSparse() const { return !array; }
	
	const char* getName() const 
		{ return "Maple Text Format"; }// LinBox__FORMAT_MAPLE_H::name; }
	
	const char* shortName() const 
		{ return "maple"; }// LinBox__FORMAT_MAPLE_H::shortname; 
};


}

#endif // __LINBOX_format_maple_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
