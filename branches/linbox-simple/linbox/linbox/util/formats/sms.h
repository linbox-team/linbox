/* MatrixStreamReader for sms matrix format
 * Dan Roche, 1-12-05
 * This file should not be included directly.  Rather, include it in
 * linbox/blackbox/matrix-stream-readers.h
 */

#ifndef __SMS_H
#define __SMS_H

#include <sstream>

namespace LinBox__SMS_H 
	{ const char* name = "SMS Sparse Integer Matrix Format"; }

namespace LinBox {

template <class Field>
class SMSReader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
	int _base;
	
    protected:

	MatrixStreamError initImpl() {
		bool retGood;
		
		try {
		    	if( !this->readSomeWhiteSpace() ||
			    !this->readObject( this->_m ) ||
			    !this->readWhiteSpace() ||
		    	    !this->readObject( this->_n ) ||
			    !this->readWhiteSpace() ) return NO_FORMAT;
			this->knowM = this->knowN = true;
			if( this->sin->get() != 'M' ) return NO_FORMAT;
			retGood = this->readBreaks();
		} catch( MatrixStreamError e ) {
			return e;
		}

		if( this->_m < 1 || this->_n < 1 ) return BAD_FORMAT;

		if( retGood) return GOOD;
		else return NO_FORMAT;
	}

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
		bool retGood;
		
		try {
			if( !this->readSomeWhiteSpace() ||
			    !this->readObject( m ) ||
			    !this->readWhiteSpace() ||
			    !this->readObject( n ) ||
			    !this->readWhiteSpace() ||
			    !readElement( v ) ) return BAD_FORMAT;
			if( m == 0 && n == 0 && this->ms->getField().isZero(v) )
				return END_OF_MATRIX;
			retGood = this->readBreaks();
		} catch( MatrixStreamError e ) { return e; }

		m -= _base;
		n -= _base;
		if( m < 0 || m >= this->_m || n < 0 || n >= this->_n )
			return BAD_FORMAT;

		if( retGood ) return GOOD;
		else return BAD_FORMAT;
	}

    public:
    	SMSReader( int base = 1 ) {
		_base = base;
	}

	const char* getName() const { return LinBox__SMS_H::name; }

	bool isSparse() const { return true; }

};

}

#endif
