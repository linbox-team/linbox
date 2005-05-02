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
    	int _base, currentRow, colsLeft;
	bool begin;

    protected:

	MatrixStreamError initImpl() {
		try {
			if( !readSomeWhiteSpace() ||
			    !readObject( _m ) ||
			    !readWhiteSpace() ||
			    !readObject( _n ) ||
			    !readBreaks() ) return NO_FORMAT;
		} catch( MatrixStreamError e ) { return e; }
		if( _m < 1 || _n < 1 ) return BAD_FORMAT;
		currentRow = -1;
		colsLeft = 0;
		begin = true;
		return GOOD;
	}

	MatrixStreamError nextTripleImpl( int& m, int& n, Element& v ) {
		try {
			while( colsLeft == 0 ) {
				if( ++currentRow == _m ) return END_OF_MATRIX;
				if( begin ) begin = false;
				else if( !readBreaks() ) return BAD_FORMAT;
				if( !readObject( colsLeft ) ) return BAD_FORMAT;
			}
			if( !readWhiteSpace() ||
			    !readObject(n) ||
			    !readWhiteSpace() ||
			    !readElement(v) ) return BAD_FORMAT;
			n -= _base;
			m = currentRow;
			--colsLeft;
		} catch( MatrixStreamError e ) { return e; }
		if( m < 0 || m >= _m ||
		    n < 0 || n >= _n ) return BAD_FORMAT;
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
