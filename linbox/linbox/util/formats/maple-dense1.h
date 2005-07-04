/* Format file for dense maple matrices defined in the LinearAlgebra style:
 * Matrix(rows,columns,[[entry,entry,entry...],[entry,entry...],...],...)
 * written by Dan Roche, 1-18-05
 */

#ifndef __MAPLE_DENSE_1_H
#define __MAPLE_DENSE_1_H

#include <string>
#include <linbox/util/matrix-stream.h>

namespace __LinBox_MAPLE_DENSE_1 
	{const char* name = "Dense Maple LinearAlgebra package matrix format";}

namespace LinBox {

template<class Field>
class MapleDense1Reader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
	std::vector<char*> tokens;
	size_t currentM, currentN;

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
	    char t;
	    try {
		std::stringstream tempS;
		t = (*this->readUntil(tokens,&tempS))[1];
		if( t != ',' && t != ']' ) return BAD_FORMAT;
		this->ms->getField().read(tempS,v);
		//if( tempS.fail() ) return BAD_FORMAT;
		if( !this->readSomeWhiteSpace(true) ) return BAD_FORMAT;
		m = currentM;
		n = currentN++;
		if( t == ']' ) {
			if( !this->readSomeWhiteSpace(true) ) return BAD_FORMAT;
			++currentM;
			if( !this->knowN || this->_n < currentN ) {
				this->_n = currentN;
				this->knowN = true;
			}
			currentN = 0;
			this->sin->get(t);
			if( t == ']' ) {
				this->atEnd = true;
				if( !this->knowM || this->_m < currentM ) {
					this->_m = currentM;
					this->knowM = true;
				}
				tokens.pop_back();
				tokens.pop_back();
				this->readUntil(tokens);
				return GOOD;
			}
			else if( t != ',' ) return BAD_FORMAT;
			if( !this->readSomeWhiteSpace(true) ) return BAD_FORMAT;
			if( this->sin->get() != '[' ) return BAD_FORMAT;
		}
		if( !this->readSomeWhiteSpace(true) ) return BAD_FORMAT;
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
		    if( !this->readSomeWhiteSpace(true) ) return NO_FORMAT;
		    std::stringstream tempS;
		    if( (*this->readUntil(line1tokens,&tempS, 160 ))[1] != '(' )
		    	return NO_FORMAT;
		    while( tempS >> temp );
		    if( temp != "Matrix" ||
			!this->readSomeWhiteSpace(true) ) return NO_FORMAT;
		    this->sin->get(t);
		    if( t != '[' ) {
			this->sin->putback(t);
		    	if( !this->readObject(this->_m) ||
			    !this->readSomeWhiteSpace(true) ||
			    this->sin->get() != ',' ||
			    !this->readSomeWhiteSpace(true) ) return NO_FORMAT;
			this->sin->get(t);
		        if( t != '[' ) {
			    this->sin->putback(t);
			    if( !this->readObject(this->_n) ||
			        !this->readSomeWhiteSpace(true) ||
			    	this->sin->get() != ',' ||
				!this->readSomeWhiteSpace(true) ) return NO_FORMAT;
			    this->sin->get(t);
			    if( t != '[' ) return NO_FORMAT;
			}
			else this->_n = this->_m;
			this->knowM = this->knowN = true;
		    }
		    if( !this->readSomeWhiteSpace(true) ||
		    	this->sin->get() != '[' ||
			!this->readSomeWhiteSpace(true) ) return NO_FORMAT;
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
