#include "linbox/util/formats/matrix-stream-readers.h"

namespace LinBox {

template<class Field>
bool MatrixStream<Field>::readWhiteSpace()
{
	char c;
	in.get(c);

	while(true) {
		if( in.eof() ) {
			return false;
		}

		switch(c) {
		    case '\n':
			++lineNumber;
		    case '\r':
		    	if( in.peek() != '\n' ) ++lineNumber;
		    case ' ':
		    case '\t':
		    case '\v':
		    case '\f':
			break;
		    default:
		    	in.putback(c);
			return true;
		}

		in.get(c);
	}
}

template<class Field>
void MatrixStreamReader<Field>::saveTriple(size_t m, size_t n, const Element& v ) {
	static std::pair<std::pair<size_t,size_t>,Element> temp;
	temp.first.first = m;
	temp.first.second = n;
	temp.second = v;
	savedTriples.push(temp);
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::init
	(char* firstLine, std::istream* i, MatrixStream<Field>* m )
{
	if( !i || !m || !firstLine ) throw "Bad istream or MatrixStream";
	sin = i;
	ms = m;
	return initImpl(firstLine);
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::nextTriple
	(size_t& m, size_t& n, Element& v) 
{
	if( savedTriples.size() == 0 ) {
		if( atEnd ) {
			if( lastError <= GOOD ) lastError = END_OF_MATRIX;
			return lastError;
		}
		if( lastError > GOOD ) return lastError;
		lastError =  nextTripleImpl(m,n,v);
		return lastError;
	}
	m = savedTriples.front().first.first;
	n = savedTriples.front().first.second;
	v = savedTriples.front().second;
	savedTriples.pop();
	return GOOD;
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::saveNext() {
	if( lastError > GOOD ) return lastError;
	if( atEnd ) {
		lastError = END_OF_MATRIX;
		return lastError;
	}
	size_t m, n;
	Element v;
	lastError = nextTripleImpl(m,n,v);
	if( lastError <= GOOD ) saveTriple(m,n,v);
	return lastError;
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::getRows(size_t& m) {
	MatrixStreamError toRet = GOOD;
	while( !knowM ) {
		if( atEnd ) return END_OF_MATRIX;
		toRet = saveNext();
		if( toRet > GOOD ) return toRet;
	}
	m = _m;
	return toRet;
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::getColumns(size_t& n) {
	MatrixStreamError toRet = GOOD;
	while( !knowN ) {
		if( atEnd ) return END_OF_MATRIX;
		toRet = saveNext();
		if( toRet > GOOD ) return toRet;
	}
	n = _n;
	return toRet;
}

template<class Field>
void MatrixStream<Field>::init() {
	lineNumber = 1;

	//Skip comments
	readWhiteSpace();
	while( !in.eof() && in.peek() == '#' ) {
		char c;
		while( in.get(c) ) {
			if( c == '\n' ) break;
			if( c == '\r' ) {
				if( in.peek() == '\n' )
					in.get();
				break;
			}
		}
		++lineNumber;
		readWhiteSpace();
	}

	//Get first line
	firstLine = new char[FIRST_LINE_LIMIT];
	in.getline(firstLine,FIRST_LINE_LIMIT);

	//Initialize readers
	currentError = NO_FORMAT;
	__MATRIX_STREAM_READERDEFS
	delete [] firstLine;

	if( !reader ) return;
	else if( currentError > GOOD )
		errorLineNumber = lineNumber;
}

template<class Field>
void MatrixStream<Field>::addReader( MatrixStreamReader<Field>* r ) {
	if( currentError == GOOD ) {
		delete r;
		return;
	}

	MatrixStreamError mse = r->init( firstLine, &in, this );
	if( mse < currentError ) {
		if( reader ) delete reader;
		reader = r;
		currentError = mse;
	}
}

template<class Field>
MatrixStream<Field>::MatrixStream(const Field& fld, std::istream& i )
	:reader(NULL),in(i),f(fld)
{
	init();
	if( currentError > GOOD ) throw currentError;
}

template<class Field>
bool MatrixStream<Field>::nextTriple(size_t& m, size_t& n, Element& v) {
	if( currentError > GOOD ) return false;

	do {
		currentError = reader->nextTriple(m,n,v);
	} while( f.isZero(v) && currentError == GOOD );

	if( currentError != GOOD )
		errorLineNumber = lineNumber;

	return( currentError == GOOD );
}

template<class Field>
bool MatrixStream<Field>::getRows(size_t& m) {
	MatrixStreamError mse = reader->getRows(m);

	if( currentError > GOOD ) 
		return (mse == GOOD);
	else if( mse > GOOD ) {
		currentError = mse;
		errorLineNumber = lineNumber;
		return false;
	}
	else return true;
}

template<class Field>
bool MatrixStream<Field>::getColumns(size_t& n) {
	MatrixStreamError mse = reader->getColumns(n);

	if( currentError > GOOD ) 
		return (mse == GOOD);
	else if( mse > GOOD ) {
		currentError = mse;
		errorLineNumber = lineNumber;
		return false;
	}
	else return true;
}

template<class Field>
bool MatrixStream<Field>::getDimensions( size_t& m, size_t& n ) {
	return( getRows(m) && getColumns(n) );
}

template<class Field>
MatrixStreamError MatrixStream<Field>::reportError
		( const char* func, int line ) const
	{
        	std::cerr << std::endl
		         << "ERROR (" << func << ":" << line << "): "
			 << "Problem reading matrix:" << std::endl;
		switch( getError() ) {
		    case END_OF_MATRIX:
		    	std::cerr << "There is no more data in the matrix file.";
			break;
		    case END_OF_FILE:
		    	std::cerr << "An EOF was encountered unexpectedly in reading the data.";
			break;
		    case BAD_FORMAT:
		    	std::cerr << "There is a formatting error in the matrix.";
			break;
		    case NO_FORMAT:
		    	std::cerr << "The matrix format is not recognized or supported.";
			break;
		    case GOOD: break;
		    default: break;
		}
		std::cerr << std::endl << "At line number: " << lineNumber << std::endl
		          << "Matrix format is " << getFormat() << std::endl;
		return currentError;
	}

template<class Field>
int MatrixStream<Field>::getLineNumber() const {
	if( currentError > GOOD ) return errorLineNumber;
	else return lineNumber;
}

} // end of namespace LinBox
