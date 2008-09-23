/* dense.h
 * MatrixStreamReader specialization for matrices in the generic dense format:
 * 1st line: #rows #cols
 * Subsequent lines: every entry, row by row.
 */

#ifndef __FORMAT_DENSE_H
#define __FORMAT_DENSE_H

/*
namespace LinBox__FORMAT_DENSE_H
	{ static const char* name = "Generic Dense Format";
	  static const char* shortname = "dense"; }
	  */

namespace LinBox {

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
		currentRow = currentCol = -1;
	}

	bool isSparse() const { return false; }
    	
	const char* getName() const { return "Generic Dense Format"; }//LinBox__FORMAT_DENSE_H::name; }
	
	const char* shortName() const 
		{ return "dense"; }//LinBox__FORMAT_DENSE_H::shortname; }

};

}

#endif // __FORMAT_DENSE_H
