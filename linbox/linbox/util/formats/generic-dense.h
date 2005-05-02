/* dense.h
 * MatrixStreamReader specialization for matrices in the generic dense format:
 * 1st line: #rows #cols
 * Subsequent lines: every entry, row by row.
 */

#ifndef __FORMAT_DENSE_H
#define __FORMAT_DENSE_H

namespace LinBox__FORMAT_DENSE_H
	{ const char* name = "Generic Dense Format"; }

namespace LinBox {

template<class Field>
class DenseReader :public MatrixStreamReader<Field> {
    public:
	using MatrixStreamReader<Field>:: readSomeWhiteSpace; 
	using MatrixStreamReader<Field>:: readWhiteSpace; 
	using MatrixStreamReader<Field>:: readObject;
	using MatrixStreamReader<Field>:: readBreaks;
	using MatrixStreamReader<Field>:: readUntil;
	using MatrixStreamReader<Field>:: atEnd;
	using MatrixStreamReader<Field>:: ms;
	using MatrixStreamReader<Field>:: sin;
	using MatrixStreamReader<Field>:: _m;
	using MatrixStreamReader<Field>:: _n;
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int currentRow, currentCol;
    protected:

	MatrixStreamError nextTripleImpl( int& m, int& n, Element& v ) {
	    if( currentRow == _m ) return END_OF_MATRIX;
	    m = currentRow;
	    n = currentCol;
	    try {
		if( (currentCol || currentRow) && !readWhiteSpace(true) )
			return BAD_FORMAT;
		if( !readElement(v) ) return BAD_FORMAT;
	    } catch( MatrixStreamError e ) { return e; }
	    if( ++currentCol == _n ) {
	    	++currentRow;
		currentCol = 0;
	    }
	    return GOOD;
	}

	MatrixStreamError initImpl() {
		try {
			int temp;
			if( !readSomeWhiteSpace() ||
			    !readObject( _m ) ||
			    !readWhiteSpace() ||
			    !readObject( _n ) ||
			    !readWhiteSpace(temp) ||
			    temp == 0 ) return NO_FORMAT;
		} catch( MatrixStreamError e ) { return e; }
		if( _m < 1 || _n < 1 ) return BAD_FORMAT;
		currentRow = currentCol = 0;
		return GOOD;
	}

    public:
    	DenseReader() {
		currentRow = currentCol = -1;
	}

	bool isSparse() const { return false; }
    	
	const char* getName() const { return LinBox__FORMAT_DENSE_H::name; }

};

}

#endif // __FORMAT_DENSE_H
