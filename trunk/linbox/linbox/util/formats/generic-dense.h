/* dense.h
 * MatrixStreamReader specialization for matrices in the generic dense format:
 * 1st line: #rows #cols
 * Subsequent lines: every entry, row by row.
 */

#ifndef __FORMAT_DENSE_H
#define __FORMAT_DENSE_H

namespace LinBox__FORMAT_DENSE_H
	{ const char* name = "Generic Dense Format";
	  const char* shortname = "dense"; }

namespace LinBox {

template<class Field>
class DenseReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	size_t currentRow, currentCol;
    protected:

    	MatrixStreamError initImpl(char* firstLine) {
		char* restLine;

		// Read m
		this->_m = strtoul(firstLine,&restLine,0);
		if( this->_m == 0 && restLine == firstLine )
			return NO_FORMAT;
		firstLine = restLine;

		// Read n
		this->_n = strtoul(firstLine,&restLine,0);
		if( this->_n == 0 && restLine == firstLine )
			return NO_FORMAT;
		firstLine = restLine;

		// Check whitespace for rest of line
		++firstLine;
		while( *firstLine && isspace(*firstLine) )
			++firstLine;
		if( *firstLine ) return BAD_FORMAT;

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
		currentRow = currentCol = -1;
	}

	bool isSparse() const { return false; }
    	
	const char* getName() const { return LinBox__FORMAT_DENSE_H::name; }
	
	const char* shortName() const 
		{ return LinBox__FORMAT_DENSE_H::shortname; }

};

}

#endif // __FORMAT_DENSE_H
