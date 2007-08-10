/* sparse-row.h
 * MatrixStreamReader specialization for matrices in the sparse row format:
 * 1st line: "50 60 S" // #rows #cols letter-S
 * Next lines: (# of entries in row) (column index of 1st entry) (value of 1st entry) ...
 * One next line for each row in matrix.
 */

#ifndef __FORMAT_SPARSE_ROW_H
#define __FORMAT_SPARSE_ROW_H

#include <cstdlib>

namespace LinBox__FORMAT_SPARSE_ROW_H
	{ const char* name = "Sparse Row Format";
	  const char* shortname = "sparserow"; }

namespace LinBox {

template<class Field>
class SparseRowReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int _base;
	size_t currentRow, colsLeft;

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

		// Read "S"
		while( firstLine[i] && isspace(firstLine[i]) )
			++i;
		if( !firstLine[i] || (firstLine[i] != 'S' &&
		                      firstLine[i] != 's'   ) )
			return NO_FORMAT;

		// Check whitespace for rest of line
		++i;
		while( firstLine[i] && isspace(firstLine[i]) )
			++i;
		if( firstLine[i] ) return BAD_FORMAT;

		this->knowM = this->knowN = true;

		currentRow = -1;
		colsLeft = 0;
		return GOOD;
	}

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
	        while( colsLeft == 0 ) {
	          	if( ++currentRow == this->_m ) return END_OF_MATRIX;
			this->ms->readWhiteSpace();
	          	*(this->sin) >> colsLeft;
	          	if( this->sin->eof() ) return END_OF_FILE;
	          	if( !this->sin->good() ) return BAD_FORMAT;
	        }
       
		this->ms->readWhiteSpace();
	        *(this->sin) >> n;
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;
       
		this->ms->readWhiteSpace();
	        this->ms->getField().read(*(this->sin),v);
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;

		n -= _base;
		m = currentRow;
		--colsLeft;

		if( m < 0 || m >= this->_m ||
		    n < 0 || n >= this->_n ) return BAD_FORMAT;

		return GOOD;
	}

    public:
    	SparseRowReader( int base = 0 ) {
		_base = base;
		currentRow = colsLeft = -1;
	}

	const char* getName() const {return LinBox__FORMAT_SPARSE_ROW_H::name;}
	const char* shortName() const
	{ return LinBox__FORMAT_SPARSE_ROW_H::shortname; }

	bool isSparse() const { return true; }
};

}

#endif // __FORMAT_SPARSE_ROW_H
