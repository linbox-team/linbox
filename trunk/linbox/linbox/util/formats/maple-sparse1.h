/* Format file for sparse maple matrices defined in the LinearAlgebra style:
 * Matrix(rows,columns,{(i,j)=v,(i,j)=v,...},...)
 * written by Dan Roche, 1-18-05
 */

#ifndef __MAPLE_SPARSE_1_H
#define __MAPLE_SPARSE_1_H

#include <string>
#include <queue>
using std::string;
using std::queue;
using std::pair;

namespace __LinBox_MAPLE_SPARSE_1 
	{const char* name = "Sparse Maple LinearAlgebra package matrix format";}

namespace LinBox {

template<class Field>
class MapleSparse1Reader :public MatrixStreamReader<Field> {
    public:
    	typedef typename MatrixStreamReader<Field>::Element Element;
    private:
	vector<char*> tokens;
	int currentM, currentN;

	MatrixStreamError nextTripleImpl( int& m, int& n, Element& v ) {
	    char t;
	    try {
	    	if( sin->get() != '(' ||
		    !readSomeWhiteSpace(true) ||
		    !readObject(m) ||
		    !readSomeWhiteSpace(true) ||
		    sin->get() != ',' ||
		    !readSomeWhiteSpace(true) ||
		    !readObject(n) ||
		    !readSomeWhiteSpace(true) ||
		    sin->get() != ')' ||
		    !readSomeWhiteSpace(true) ||
		    sin->get() != '=' ||
		    !readSomeWhiteSpace(true) ) return BAD_FORMAT;
		if( m > currentM ) currentM = m;
		if( n > currentN ) currentN = n;
		stringstream tempS;
		t = (*readUntil(tokens,&tempS))[1];
		ms->getField().read(tempS,v);
		//if( tempS.fail() ) return BAD_FORMAT;
		if( !readSomeWhiteSpace(true) ) return BAD_FORMAT;
		if( t == '}' ) {
			atEnd = true;
			tokens.pop_back();
			tokens.pop_back();
			readUntil(tokens);
			if( _m == -1 ) _m = currentM;
			if( _n == -1 ) _n = currentN;
		}
	    }
	    catch( MatrixStreamError e ) { return e; }

	    --m;
	    --n;

	    return GOOD;
	}
    	
	MatrixStreamError initImpl() {
		string temp;
		char t;
		vector<char*> line1tokens;
		line1tokens.push_back("\0(");
		line1tokens.push_back("\0\n");
		try {
		    if( !readSomeWhiteSpace(true) ) return NO_FORMAT;
		    stringstream tempS;
		    if( (*readUntil(line1tokens,&tempS, 160 ))[1] != '(' )
		    	return NO_FORMAT;
		    while( tempS >> temp );
		    if( temp != "Matrix" ||
			!readSomeWhiteSpace(true) ) return NO_FORMAT;
		    sin->get(t);
		    if( t != '{' ) {
			sin->putback(t);
		    	if( !readObject(_m) ||
			    !readSomeWhiteSpace(true) ||
			    sin->get() != ',' ||
			    !readSomeWhiteSpace(true) ) return NO_FORMAT;
			sin->get(t);
		        if( t != '{' ) {
			    sin->putback(t);
			    if( !readObject(_n) ||
			        !readSomeWhiteSpace(true) ||
			    	sin->get() != ',' ||
				!readSomeWhiteSpace(true) ) return NO_FORMAT;
			    sin->get(t);
			    if( t != '{' ) return NO_FORMAT;
			}
			else _n = _m;
		    }
		    if( !readSomeWhiteSpace(true) ) return NO_FORMAT;
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
