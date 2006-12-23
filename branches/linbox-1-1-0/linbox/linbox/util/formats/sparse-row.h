/* sparse-row.h
 * MatrixStreamReader specialization for matrices in the sparse row format:
 * 1st line: rows columns
 * next lines: (# of entries in row) (column index of 1st entry) (value of 1st entry) ...
 * # of lines of data = # of lines in matrix
 */

#ifndef __FORMAT_SPARSE_ROW_H
#define __FORMAT_SPARSE_ROW_H

namespace LinBox__FORMAT_SPARSE_ROW_H
	{ const char* name = "Sparse Row Format"; }

namespace LinBox {

template<class Field>
class SparseRowReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int _base;
	size_t currentRow, colsLeft;
	bool begin;

    protected:

	MatrixStreamError initImpl() {
		try {
			if( !this->readSomeWhiteSpace() ||
			    !this->readObject( this->_m ) ||
			    !this->readWhiteSpace() ||
			    !this->readObject( this->_n ) ||
			    !this->readBreaks() ) return NO_FORMAT;
			    this->knowM = this->knowN = true;
		} catch( MatrixStreamError e ) { return e; }
		if( this->_m < 1 || this->_n < 1 ) return BAD_FORMAT;
		currentRow = -1;
		colsLeft = 0;
		begin = true;
		return GOOD;
	}

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
		try {
			while( colsLeft == 0 ) {
				if( ++currentRow == this->_m ) return END_OF_MATRIX;
				if( begin ) begin = false;
				else if( !this->readBreaks() ) return BAD_FORMAT;
				if( !this->readObject( colsLeft ) ) return BAD_FORMAT;
			}
			if( !this->readWhiteSpace() ||
			    !this->readObject(n) ||
			    !this->readWhiteSpace() ||
			    !readElement(v) ) return BAD_FORMAT;
			n -= _base;
			m = currentRow;
			--colsLeft;
		} catch( MatrixStreamError e ) { return e; }
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

	bool isSparse() const { return true; }

};

}

#endif // __FORMAT_SPARSE_ROW_H
