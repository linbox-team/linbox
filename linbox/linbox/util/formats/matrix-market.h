/* matrix-market-array.h
 * MatrixStreamReader specialization for matrices in the MatrixMarket coordinate
 * format.
 * Dan Roche, 1-25-05
 */

#ifndef __FORMAT_MATRIX_MARKET_H
#define __FORMAT_MATRIX_MARKET_H

#include <string>
#include <cctype>

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
	using MatrixStreamReader<Field>:: readSomeWhiteSpace; 
	using MatrixStreamReader<Field>:: readWhiteSpace; 
	using MatrixStreamReader<Field>:: readObject;
	using MatrixStreamReader<Field>:: readBreaks;
	using MatrixStreamReader<Field>:: readUntil;
	using MatrixStreamReader<Field>:: atEnd;
	using MatrixStreamReader<Field>:: moreData;
	using MatrixStreamReader<Field>:: ms;
	using MatrixStreamReader<Field>:: sin;
	using MatrixStreamReader<Field>:: _m;
	using MatrixStreamReader<Field>:: _n;
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int entriesLeft, currentCol, currentRow;
	bool array;
	bool pattern;
	bool symmetric;
    protected:
    	MatrixStreamError nextTripleImpl( int& m, int& n, Element& v ) {
		if( array ) {
			if( currentCol == _n+1 ) return END_OF_MATRIX;
			n = currentCol;
			m = currentRow;
			if( ++currentRow == _m+1 ) {
				++currentCol;
				currentRow = (symmetric ? currentCol : 1);
			}
		}
		else if( --entriesLeft < 0 ) return END_OF_MATRIX;

		try {
		    if( !readBreaks() ||
		     (  !array &&
		      ( !readObject(m) ||
			!readWhiteSpace() ||
			!readObject(n) ||
		       (!pattern &&
		        !readWhiteSpace()))) ||
		     (  !pattern &&
		        !readElement(v)  ) ) return BAD_FORMAT;
		}
		catch( MatrixStreamError e ) {return e;}
		if( pattern ) ms->getField().init(v,(integer)1);
		--m;
		--n;
		if( m < 0 || m >= _m || n < 0 || n >= _n )
			return BAD_FORMAT;
		if( symmetric && (m != n) ) saveTriple(n,m,v);
		return GOOD;
	}

	MatrixStreamError initImpl() {
	    std::string s;
	    try {
		if( !readObject(s) ||
		    (s != "%%MatrixMarket") ||
		    !readWhiteSpace() ||
		    !readObject(s) ||
		    !equalCaseInsensitive(s,"matrix") ||
		    !readWhiteSpace() ||
		    !readObject(s) ) return NO_FORMAT;
		if( equalCaseInsensitive(s,"array") ) array = true;
		else if( equalCaseInsensitive(s,"coordinate") ) array = false;
		else return NO_FORMAT;
		if( !readWhiteSpace() ||
		    !readObject(s) ) return NO_FORMAT;
		pattern = equalCaseInsensitive(s,"pattern");
		if( !readWhiteSpace() ||
		    !readObject(s) ) return NO_FORMAT;
		if( equalCaseInsensitive(s,"symmetric") ) symmetric = true;
		else if( equalCaseInsensitive(s,"general") ) symmetric = false;
		else return NO_FORMAT;
		if( !readBreaks() ) return NO_FORMAT;
		char c;
		sin->get(c);
		while( c == '%' ) {
		    if( !readUntil('\n') ) return NO_FORMAT;
		    sin->get(c);
		    if( sin->eof() ) {
		    	if( !moreData() ) return END_OF_FILE;
			sin->get(c);
		    }
		}
		sin->putback(c);
		if( !readSomeWhiteSpace(true) ||
		    !readObject(_m) ||
		    !readWhiteSpace() ||
		    !readObject(_n) ) return NO_FORMAT;
		if( !array && !( readWhiteSpace() && readObject(entriesLeft) ) )
			return NO_FORMAT;
	    } catch( MatrixStreamError e ) { return e; }
	    if( array && pattern ) return BAD_FORMAT;
	    if( symmetric && (_m != _n) ) return BAD_FORMAT;
	    if( _m < 1 || _n < 1 ) return BAD_FORMAT;
	    if( !array && (entriesLeft < 0 || entriesLeft > _m*_n ) )
		return BAD_FORMAT;
	    currentRow = currentCol = 1;
	    return GOOD;
	}

    public:
    	MatrixMarketReader() {
		entriesLeft = currentCol = currentRow = -1;
	}

	bool isSparse() const { return !array; }
	
	const char* getName() const 
		{ return LinBox__FORMAT_MATRIX_MARKET_H::name; }
};

}

#endif
