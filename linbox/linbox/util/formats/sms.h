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
	int _base;
	
    protected:

	MatrixStreamError initImpl() {
		bool retGood;
		
		try {
		    	if( !readSomeWhiteSpace() ||
			    !readObject( _m ) ||
			    !readWhiteSpace() ||
		    	    !readObject( _n ) ||
			    !readWhiteSpace() ) return NO_FORMAT;
			if( sin->get() != 'M' ) return NO_FORMAT;
			retGood = readBreaks();
		} catch( MatrixStreamError e ) {
			return e;
		}

		if( _m < 1 || _n < 1 ) return BAD_FORMAT;

		if( retGood) return GOOD;
		else return NO_FORMAT;
	}

	MatrixStreamError nextTripleImpl( int& m, int& n, Element& v ) {
		bool retGood;
		
		try {
			if( !readSomeWhiteSpace() ||
			    !readObject( m ) ||
			    !readWhiteSpace() ||
			    !readObject( n ) ||
			    !readWhiteSpace() ||
			    !readElement( v ) ) return BAD_FORMAT;
			if( m == 0 && n == 0 && ms->getField().isZero(v) )
				return END_OF_MATRIX;
			retGood = readBreaks();
		} catch( MatrixStreamError e ) { return e; }

		m -= _base;
		n -= _base;
		if( m < 0 || m >= _m || n < 0 || n >= _n )
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
