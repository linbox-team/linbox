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
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	size_t currentRow, currentCol;
    protected:

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
	    if( currentRow == this->_m ) return END_OF_MATRIX;
	    m = currentRow;
	    n = currentCol;
	    try {
		if( (currentCol || currentRow) && !this->readWhiteSpace(true) )
			return BAD_FORMAT;
		if( !readElement(v) ) return BAD_FORMAT;
	    } catch( MatrixStreamError e ) { return e; }
	    if( ++currentCol == this->_n ) {
	    	++currentRow;
		currentCol = 0;
	    }
	    return GOOD;
	}

	MatrixStreamError initImpl() {
		try {
			int temp=0;
			if( !this->readSomeWhiteSpace() ||
			    !this->readObject( this->_m ) ||
			    !this->readWhiteSpace() ||
			    !this->readObject( this->_n ) ||
			    !this->readWhiteSpace(temp) ||
			    temp == 0 ) return NO_FORMAT;
		} catch( MatrixStreamError e ) { return e; }
		if( this->_m < 1 || this->_n < 1 ) return BAD_FORMAT;
		this->knowM = this->knowN = true;
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
