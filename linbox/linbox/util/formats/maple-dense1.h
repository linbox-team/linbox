/* Format file for dense maple matrices defined in the LinearAlgebra style:
 * Matrix(rows,columns,[[entry,entry,entry...],[entry,entry...],...],...)
 * written by Dan Roche, 1-18-05
 */

#ifndef __MAPLE_DENSE_1_H
#define __MAPLE_DENSE_1_H

#include <string>
#include <queue>

namespace __LinBox_MAPLE_DENSE_1 
	{const char* name = "Dense Maple LinearAlgebra package matrix format";}

namespace LinBox {

template<class Field>
class MapleDense1Reader :public MatrixStreamReader<Field> {
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
	std::vector<char*> tokens;
	int currentM, currentN;

	MatrixStreamError nextTripleImpl( int& m, int& n, Element& v ) {
	    char t;
	    try {
		std::stringstream tempS;
		t = (*readUntil(tokens,&tempS))[1];
		if( t != ',' && t != ']' ) return BAD_FORMAT;
		ms->getField().read(tempS,v);
		//if( tempS.fail() ) return BAD_FORMAT;
		if( !readSomeWhiteSpace(true) ) return BAD_FORMAT;
		m = currentM;
		n = currentN++;
		if( t == ']' ) {
			if( !readSomeWhiteSpace(true) ) return BAD_FORMAT;
			++currentM;
			if( _n == -1 || _n < currentN ) _n = currentN;
			currentN = 0;
			sin->get(t);
			if( t == ']' ) {
				atEnd = true;
				if( _m == -1 || _m < currentM ) _m = currentM;
				tokens.pop_back();
				tokens.pop_back();
				readUntil(tokens);
				return GOOD;
			}
			else if( t != ',' ) return BAD_FORMAT;
			if( !readSomeWhiteSpace(true) ) return BAD_FORMAT;
			if( sin->get() != '[' ) return BAD_FORMAT;
		}
		if( !readSomeWhiteSpace(true) ) return BAD_FORMAT;
	    }
	    catch( MatrixStreamError e ) { return e; }
	    return GOOD;
	}

	MatrixStreamError initImpl() {
		std::string temp;
		char t;
		std::vector<char*> line1tokens;
		line1tokens.push_back("\0(");
		line1tokens.push_back("\0\n");
		try {
		    if( !readSomeWhiteSpace(true) ) return NO_FORMAT;
		    std::stringstream tempS;
		    if( (*readUntil(line1tokens,&tempS, 160 ))[1] != '(' )
		    	return NO_FORMAT;
		    while( tempS >> temp );
		    if( temp != "Matrix" ||
			!readSomeWhiteSpace(true) ) return NO_FORMAT;
		    sin->get(t);
		    if( t != '[' ) {
			sin->putback(t);
		    	if( !readObject(_m) ||
			    !readSomeWhiteSpace(true) ||
			    sin->get() != ',' ||
			    !readSomeWhiteSpace(true) ) return NO_FORMAT;
			sin->get(t);
		        if( t != '[' ) {
			    sin->putback(t);
			    if( !readObject(_n) ||
			        !readSomeWhiteSpace(true) ||
			    	sin->get() != ',' ||
				!readSomeWhiteSpace(true) ) return NO_FORMAT;
			    sin->get(t);
			    if( t != '[' ) return NO_FORMAT;
			}
			else _n = _m;
		    }
		    if( !readSomeWhiteSpace(true) ||
		    	sin->get() != '[' ||
			!readSomeWhiteSpace(true) ) return NO_FORMAT;
		}
		catch( MatrixStreamError e ) { return e; }

		return GOOD;
	}
    	
    public:
    	MapleDense1Reader() :tokens(8)
	{
		currentM = currentN = 0;
		tokens[0] = "\"\"";
		tokens[1] = "''";
		tokens[2] = "{}";
		tokens[3] = "[]";
		tokens[4] = "()";
		tokens[5] = "<>";
		tokens[6] = "\0\n";
		tokens[7] = "\0,";
	}

    	const char* getName() const { return __LinBox_MAPLE_DENSE_1::name; }

	bool isSparse() const { return false; }
}; // end of class MapleDense1Reader

} // end of namespace LinBox

#endif
