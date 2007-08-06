/* MatrixStreamReader for sms matrix format
 * Dan Roche, 1-12-05
 * 1st line: "50 60 M" // #rows #cols letter-M
 * Subsequent lines: i j v // row index, col index, value
 * last line: 0 0 0
 */

#ifndef __SMS_H
#define __SMS_H

#include <cstdlib>

namespace LinBox__SMS_H
	{ const char* name = "SMS Sparse Integer Matrix Format";
	  const char* shortname = "sms"; }

namespace LinBox {

template<class Field>
class SMSReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
    	int _base;

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

		// Read "M"
		while( *firstLine && isspace(*firstLine) )
			++firstLine;
		if( !(*firstLine) || ((*firstLine) != 'M' &&
		                      (*firstLine) != 'm'   ) )
			return NO_FORMAT;

		// Check whitespace for rest of line
		++firstLine;
		while( *firstLine && isspace(*firstLine) )
			++firstLine;
		if( *firstLine ) return BAD_FORMAT;

		this->knowM = this->knowN = true;

		return GOOD;
	}

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
		this->ms->readWhiteSpace();
	        *(this->sin) >> m;
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;
       
		this->ms->readWhiteSpace();
	        *(this->sin) >> n;
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;
       
		this->ms->readWhiteSpace();
	        this->ms->getField().read(*(this->sin),v);
	        if( this->sin->eof() ) return END_OF_FILE;
	        if( !this->sin->good() ) return BAD_FORMAT;

		if( m == 0 && n == 0 ) return END_OF_MATRIX;

		m -= _base;
		n -= _base;

		if( m < 0 || m >= this->_m ||
		    n < 0 || n >= this->_n ) return BAD_FORMAT;

		return GOOD;
	}

    public:
    	SMSReader( int base = 1 ) {
		_base = base;
	}

	const char* getName() const {return LinBox__SMS_H::name;}
	const char* shortName() const
	{ return LinBox__SMS_H::shortname; }

	bool isSparse() const { return true; }
};

}

#endif // __FORMAT_SPARSE_ROW_H
