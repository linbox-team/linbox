/* matrix-market-array.h
 * MatrixStreamReader specialization for matrices in the MatrixMarket coordinate
 * format.
 * Dan Roche, 1-25-05
 */

#ifndef __FORMAT_MATRIX_MARKET_H
#define __FORMAT_MATRIX_MARKET_H

#include <string>
#include <cctype>
#include <linbox/util/matrix-stream.h>

namespace LinBox__FORMAT_MATRIX_MARKET_H
	{ const char* name = "Matrix Market Format"; }

namespace LinBox {

#ifndef __INTEGER_H
class integer;
#endif

bool equalCaseInsensitive(const std::string s1, const char* s2) {
	int len = s1.size();
	int counter = 0;
	while( counter < len && s2[counter] != '\0' &&
	       toupper(s1[counter]) == toupper(s2[counter]) ) ++counter;
	return( counter == len && s2[counter] == '\0' );
}

template<class Field>
class MatrixMarketReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int entriesLeft;
	size_t currentCol, currentRow;
	bool array;
	bool pattern;
	bool symmetric;
    protected:
    	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
		if( array ) {
			if( currentCol == this->_n+1 ) return END_OF_MATRIX;
			n = currentCol;
			m = currentRow;
			if( ++currentRow == this->_m+1 ) {
				++currentCol;
				currentRow = (symmetric ? currentCol : 1);
			}
		}
		else if( --entriesLeft < 0 ) return END_OF_MATRIX;

		try {
		    if( !this->readBreaks() ||
		     (  !array &&
		      ( !this->readObject(m) ||
			!this->readWhiteSpace() ||
			!this->readObject(n) ||
		       (!pattern &&
		        !this->readWhiteSpace()))) ||
		     (  !pattern &&
		        !readElement(v)  ) ) return BAD_FORMAT;
		}
		catch( MatrixStreamError e ) {return e;}
		if( pattern ) this->ms->getField().init(v,(integer)1);
		--m;
		--n;
		if( m < 0 || m >= this->_m || n < 0 || n >= this->_n )
			return BAD_FORMAT;
		if( symmetric && (m != n) ) saveTriple(n,m,v);
		return GOOD;
	}

	MatrixStreamError initImpl() {
	    std::string s;
	    try {
		if( !this->readObject(s) ||
		    (s != "%%MatrixMarket") ||
		    !this->readWhiteSpace() ||
		    !this->readObject(s) ||
		    !equalCaseInsensitive(s,"matrix") ||
		    !this->readWhiteSpace() ||
		    !this->readObject(s) ) return NO_FORMAT;
		if( equalCaseInsensitive(s,"array") ) array = true;
		else if( equalCaseInsensitive(s,"coordinate") ) array = false;
		else return NO_FORMAT;
		if( !this->readWhiteSpace() ||
		    !this->readObject(s) ) return NO_FORMAT;
		pattern = equalCaseInsensitive(s,"pattern");
		if( !this->readWhiteSpace() ||
		    !this->readObject(s) ) return NO_FORMAT;
		if( equalCaseInsensitive(s,"symmetric") ) symmetric = true;
		else if( equalCaseInsensitive(s,"general") ) symmetric = false;
		else return NO_FORMAT;
		if( !this->readBreaks() ) return NO_FORMAT;
		char c;
		this->sin->get(c);
		while( c == '%' ) {
		    if( !this->readUntil('\n') ) return NO_FORMAT;
		    this->sin->get(c);
		    if( this->sin->eof() ) {
		    	if( !this->moreData() ) return END_OF_FILE;
			this->sin->get(c);
		    }
		}
		this->sin->putback(c);
		if( !this->readSomeWhiteSpace(true) ||
		    !this->readObject(this->_m) ||
		    !this->readWhiteSpace() ||
		    !this->readObject(this->_n) ) return NO_FORMAT;
		this->knowM = true;
		this->knowN = true;
		if( !array && !( this->readWhiteSpace() && this->readObject(entriesLeft) ) )
			return NO_FORMAT;
	    } catch( MatrixStreamError e ) { return e; }
	    if( array && pattern ) return BAD_FORMAT;
	    if( symmetric && (this->_m != this->_n) ) return BAD_FORMAT;
	    if( this->_m < 1 || this->_n < 1 ) return BAD_FORMAT;
	    if( !array && (entriesLeft < 0 || (size_t)entriesLeft > this->_m*this->_n ) )
		return BAD_FORMAT;
	    currentRow = currentCol = 1;
	    return GOOD;
	}

    public:
    	MatrixMarketReader() {
		entriesLeft = -1;
		currentCol = currentRow = 0;
	}

	bool isSparse() const { return !array; }
	
	const char* getName() const 
		{ return LinBox__FORMAT_MATRIX_MARKET_H::name; }
};

}

#endif
