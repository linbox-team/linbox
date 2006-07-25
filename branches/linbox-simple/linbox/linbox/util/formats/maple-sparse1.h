/* Format file for sparse maple matrices defined in the LinearAlgebra style:
 * Matrix(rows,columns,{(i,j)=v,(i,j)=v,...},...)
 * written by Dan Roche, 1-18-05
 */

#ifndef __MAPLE_SPARSE_1_H
#define __MAPLE_SPARSE_1_H

#include <string>
#include <linbox/util/matrix-stream.h>

namespace __LinBox_MAPLE_SPARSE_1 
	{const char* name = "Sparse Maple LinearAlgebra package matrix format";}

namespace LinBox {

template<class Field>
class MapleSparse1Reader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
	std::vector<char*> tokens;
	size_t currentM, currentN;

	MatrixStreamError nextTripleImpl( size_t& m, size_t& n, Element& v ) {
	    char t;
	    try {
	    	if( this->sin->get() != '(' ||
		    !this->readSomeWhiteSpace(true) ||
		    !this->readObject(m) ||
		    !this->readSomeWhiteSpace(true) ||
		    this->sin->get() != ',' ||
		    !this->readSomeWhiteSpace(true) ||
		    !this->readObject(n) ||
		    !this->readSomeWhiteSpace(true) ||
		    this->sin->get() != ')' ||
		    !this->readSomeWhiteSpace(true) ||
		    this->sin->get() != '=' ||
		    !this->readSomeWhiteSpace(true) ) return BAD_FORMAT;
		if( m > currentM ) currentM = m;
		if( n > currentN ) currentN = n;
		std::stringstream tempS;
		t = (*this->readUntil(tokens,&tempS))[1];
		this->ms->getField().read(tempS,v);
		//if( tempS.fail() ) return BAD_FORMAT;
		if( !this->readSomeWhiteSpace(true) ) return BAD_FORMAT;
		if( t == '}' ) {
			this->atEnd = true;
			tokens.pop_back();
			tokens.pop_back();
			this->readUntil(tokens);
			if( !this->knowM ) {
				this->_m = currentM;
				this->knowM = true;
			}
			if( !this->knowN ) {
				this->_n = currentN;
				this->knowN = true;
			}
		}
	    }
	    catch( MatrixStreamError e ) { return e; }

	    --m;
	    --n;

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
		    if( t != '{' ) {
			this->sin->putback(t);
		    	if( !this->readObject(this->_m) ||
			    !this->readSomeWhiteSpace(true) ||
			    this->sin->get() != ',' ||
			    !this->readSomeWhiteSpace(true) ) return NO_FORMAT;
			this->knowM = true;
			this->sin->get(t);
		        if( t != '{' ) {
			    this->sin->putback(t);
			    if( !this->readObject(this->_n) ||
			        !this->readSomeWhiteSpace(true) ||
			    	this->sin->get() != ',' ||
				!this->readSomeWhiteSpace(true) ) return NO_FORMAT;
			    this->sin->get(t);
			    if( t != '{' ) return NO_FORMAT;
			}
			else this->_n = this->_m;
			this->knowN = true;
		    }
		    if( !this->readSomeWhiteSpace(true) ) return NO_FORMAT;
		}
		catch( MatrixStreamError e ) { return e; }

		return GOOD;
	}

    public:
    	MapleSparse1Reader() :tokens(8)
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

	bool isSparse() const { return true; }

    	const char* getName() const { return __LinBox_MAPLE_SPARSE_1::name; }
}; // end of class MapleSparse1Reader

} // end of namespace LinBox

#endif
