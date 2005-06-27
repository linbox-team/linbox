#include "formats/matrix-stream-readers.h"
#include <stack>

namespace LinBox {

template<class Field>
bool MatrixStreamReader<Field>::genericWSReader
	(bool allowBreaks, bool some, int* breaks, bool stopAfterBreaks ) {
	char c;
	sin->get(c);
	bool any = false;
	if( stopAfterBreaks && !breaks ) return false;

	while(true) {
		if( sin->eof() ) {
			if( !moreData() ) throw END_OF_FILE;
			sin->get(c);
		}

		switch(c) {
		    case ' ':
		    case '\t':
		    	any = true;
			break;
		    case '\n':
		    	if( stopAfterBreaks ) {
			    if( --(*breaks) == -1 ) {
			    	sin->putback(c);
				return true;
			    }
			}
			else if(breaks) --(*breaks);
			++lineNumber;
		    case '\v':
		    case '\f':
		    case '\r':
		    	if( !allowBreaks ) return false;
			any = true;
			break;
		    default:
		    	sin->putback(c);
			if( stopAfterBreaks ) return( *breaks == 0 );
			else return ( some || any );
		}

		sin->get(c);
	}
}

template<class Field>
bool MatrixStreamReader<Field>::readUntil(char c, stringstream* ss, int limit) {
	char x;
	sin->get(x);
	
	while(true) {
		if( sin->eof() ) {
			if( !moreData() ) throw END_OF_FILE;
			sin->get(x);
		}
		if( x == '\n' ) ++lineNumber;
		if( x == c ) return true;
		else if( ss ) ss->put(c);
		if( --limit == 0 ) return false;
		sin->get(x);
	}
}

template<class Field>
vector<char*>::const_iterator MatrixStreamReader<Field>::readUntil
	(const vector<char*>& cm, stringstream* ss, int limit)
{
	char x;
	sin->get(x);
	std::stack<char> matches;

	while(true) {
		if( sin->eof() ) {
			if( !moreData() ) throw END_OF_FILE;
			sin->get(x);
		}
		if( x == '\n' ) ++lineNumber;
		if( !matches.empty() && x == matches.top() ) matches.pop();
		else for( vector<char*>::const_iterator iter = cm.begin();
		         iter != cm.end(); ++iter ) {
		    	if( x == (*iter)[0] && x != '\0' ) {
				matches.push((*iter)[1]);
				break;
			}
			else if( matches.empty() && x == (*iter)[1] )
				return iter;
		}
		if( ss ) ss->put(x);
		if( --limit == 0 ) return cm.end();
		sin->get(x);
	}
}

template<class Field>
bool MatrixStreamReader<Field>::readElement( Element& x ) {
	std::streampos pos = sin->tellg();
	ms->getField().read(*sin,x);
	while( sin->rdstate() ) {
		if( !sin->eof() ) return false;
		if( !moreData() ) {
			if( !sin->fail() ) return true;
			else throw END_OF_FILE;
		}
		sin->clear();
		sin->seekg(pos);
		ms->getField().read(*sin,x);
	}
	return true;
}

template<class Field>
template<class Object>
bool MatrixStreamReader<Field>::readObject( Object& o ) {
	std::streampos pos = sin->tellg();
	(*sin) >> o;
	while( sin->rdstate() ) {
		if( !sin->eof() ) return false;
		if( !moreData() ) {
			if( !sin->fail() ) return true;
			else throw END_OF_FILE;
		}
		sin->clear();
		sin->seekg(pos);
		(*sin) >> o;
	}
	return true;
}

template<class Field>
bool MatrixStreamReader<Field>::moreData() {
	return ms->addChars(&sin);
}

template<class Field>
void MatrixStreamReader<Field>::saveTriple(int m, int n, const Element& v ) {
	static pair<pair<int,int>,Element> temp;
	temp.first.first = m;
	temp.first.second = n;
	temp.second = v;
	savedTriples.push(temp);
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::init
	(istream* i, MatrixStream<Field>* m )
{
	if( !i || !m ) throw "Bad istream or MatrixStream";
	lineNumber = 0;
	sin = i;
	ms = m;
	return initImpl();
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::nextTriple
	(int& m, int& n, Element& v) 
{
	if( savedTriples.empty() ) {
		if( atEnd ) {
			if( lastError <= GOOD ) lastError = END_OF_MATRIX;
			return lastError;
		}
		if( lastError > GOOD ) return lastError;
		lastError =  nextTripleImpl(m,n,v);
		if( lastError > GOOD ) atEnd = true;
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
	if( atEnd ) {
		if( lastError <= GOOD ) lastError = END_OF_MATRIX;
		return lastError;
	}
	int m, n;
	Element v;
	lastError = nextTripleImpl(m,n,v);
	if( lastError <= GOOD ) saveTriple(m,n,v);
	else atEnd = true;
	return lastError;
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::getRows(int& m) {
	MatrixStreamError toRet = GOOD;
	while( _m < 0 ) {
		if( atEnd ) return END_OF_MATRIX;
		toRet = saveNext();
		if( toRet > GOOD ) return toRet;
	}
	m = _m;
	return toRet;
}

template<class Field>
MatrixStreamError MatrixStreamReader<Field>::getColumns(int& n) {
	MatrixStreamError toRet = GOOD;
	while( _n < 0 ) {
		if( atEnd ) return END_OF_MATRIX;
		toRet = saveNext();
		if( toRet > GOOD ) return toRet;
	}
	n = _n;
	return toRet;
}

template<class Field>
void MatrixStream<Field>::init() {
	//Initialize readers
	__MATRIX_STREAM_READERDEFS

	directStream = false;
	if( !addChars() ) {
		currentError = END_OF_FILE;
		return;
	}
	currentError = NO_FORMAT;
	errorLineNumber = 0;
	MatrixStreamError tError;

	for( typename RVector::iterator iter = readers.begin();
	     iter != readers.end(); ) {
		tError = iter->first->init(iter->second,this);
		if( tError <= currentError ) {
			if(currentError <= GOOD )
				currentError = AMBIGUOUS_FORMAT;
			else currentError = tError;
		}
		if( tError > GOOD ) {
			int ln = iter->first->getLineNumber();
			if( ln > errorLineNumber ) errorLineNumber = ln;
			delete iter->first;
			delete iter->second;
			iter = readers.erase(iter);
		}
		else ++iter;
	}
}

template<class Field>
void MatrixStream<Field>::addReader( MatrixStreamReader<Field>* r ) {
	RPair p;
	p.first = r;
	p.second = new stringstream;
	readers.push_back(p);
}

template<class Field>
bool MatrixStream<Field>::automaticResolve() {
    MatrixStreamError tError;
    typename RVector::iterator fBad = readers.end();
    bool anyGood;
    while( readers.size() > 1 ) {
	anyGood = false;
	for( typename RVector::iterator iter = readers.begin();
	     iter != readers.end(); ) {
		tError = iter->first->saveNext();
		if( tError == END_OF_MATRIX ) {
			for(typename RVector::iterator iter2 = iter + 1;
			    iter2 != readers.end(); ) {
				delete iter2->first;
				delete iter2->second;
				iter2 = readers.erase(iter2);
			}
			for(typename RVector::iterator iter3 = readers.begin();
			    iter3 != iter; ++iter3 ) {
			    	delete iter3->first;
				delete iter3->second;
			}
			readers.erase(readers.begin(),iter);
			break;
		}
		if( tError > GOOD ) {
			int ln = iter->first->getLineNumber();
			if( ln > errorLineNumber ) errorLineNumber = ln;
			if(anyGood) {
			    delete iter->first;
			    delete iter->second;
			    iter = readers.erase(iter);
			}
			else {
			    anyGood = true;
			    fBad = iter;
			    ++iter;
			}
		}
		else {
			++iter;
			anyGood = true;
		}
	}
	if( readers.size() > 1 && fBad < readers.end() ) {
		delete fBad->first;
		delete fBad->second;
		readers.erase(fBad);
		fBad = readers.end();
	}
    }

    return !readers.empty();
}

template<class Field>
MatrixStream<Field>::MatrixStream(const Field& fld, istream& i, char delim )
	:in(i),readers(0),f(fld)
{
	delimiter = delim;
	init();
}

template<class Field>
MatrixStream<Field>::~MatrixStream() {
	for( typename RVector::iterator iter = readers.begin();
	     iter != readers.end(); ++iter ) {
		delete iter->first;
		delete iter->second;
	}
	readers.clear();
}

template<class Field>
bool MatrixStream<Field>::nextTriple(int& m, int& n, Element& v) {
	currentError = NO_FORMAT;

	if( readers.empty() ) return false;

	else if( readers.size() > 1 ) {
	    MatrixStreamError tError;
	    
	    for( typename RVector::iterator iter = readers.begin();
	         iter != readers.end(); )
	    {
		tError = iter->first->saveNext();
		if( tError <= currentError ) {
			if(currentError <= GOOD )
				currentError = AMBIGUOUS_FORMAT;
			else currentError = tError;
		}
		if( tError > GOOD ) {
			int ln = iter->first->getLineNumber();
			if( ln > errorLineNumber ) errorLineNumber = ln;
			delete iter->first;
			delete iter->second;
			iter = readers.erase(iter);
		}
		else ++iter;
	    }
	    
	    if( readers.size() == 0 ) return false;
	    else if( readers.size() > 1 && !automaticResolve() )
	    	return false;
	}
	
	currentError = readers.begin()->first->nextTriple(m,n,v);
	return( currentError <= GOOD );
}

template<class Field>
bool MatrixStream<Field>::getRows(int& m) {
	currentError = NO_FORMAT;

	if( readers.empty() ) return false;

	m = -1;
	int tempM;
	MatrixStreamError tError;
	for( typename RVector::iterator iter = readers.begin();
	     iter != readers.end(); ) {
		tError = iter->first->getRows(tempM);
		if( tError <= currentError ) currentError = tError;
		if( tError <= GOOD ) {
		    if( m == -1 || tempM == m ) {
		    	m = tempM;
			++iter;
		    }
		    else {
		    	if( !automaticResolve() ) return false;
			currentError = readers.begin()->first->getRows(m);
			return ( currentError <= GOOD );
		    }
		}
		else {
			int ln = iter->first->getLineNumber();
			if( ln > errorLineNumber ) errorLineNumber = ln;
			delete iter->first;
			delete iter->second;
			iter = readers.erase(iter);
		}
	}
	if( readers.size() > 1 ) currentError = AMBIGUOUS_FORMAT;
	return( currentError <= GOOD );
}

template<class Field>
bool MatrixStream<Field>::getColumns(int& n) {
	currentError = NO_FORMAT;

	if( readers.empty() ) return false;

	n = -1;
	int tempN;
	MatrixStreamError tError;
	for( typename RVector::iterator iter = readers.begin();
	     iter != readers.end(); ) {
		tError = iter->first->getColumns(tempN);
		if( tError <= currentError ) currentError = tError;
		if( tError <= GOOD ) {
		    if( n == -1 || tempN == n ) {
		    	n = tempN;
			++iter;
		    }
		    else {
		    	if( !automaticResolve() ) return false;
			currentError = readers.begin()->first->getColumns(n);
			return ( currentError <= GOOD );
		    }
		}
		else {
			int ln = iter->first->getLineNumber();
			if( ln > errorLineNumber ) errorLineNumber = ln;
			delete iter->first;
			delete iter->second;
			iter = readers.erase(iter);
		}
	}
	if( readers.size() > 1 ) currentError = AMBIGUOUS_FORMAT;
	return( currentError <= GOOD );
}

template<class Field>
bool MatrixStream<Field>::getDimensions( int& m, int& n ) {
	return( getRows(m) && getColumns(n) );
}

template<class Field>
int MatrixStream<Field>::getLineNumber() const {
	if( readers.empty() ) return errorLineNumber;
	else return readers.begin()->first->getLineNumber();
}

template<class Field>
bool MatrixStream<Field>::addChars(istream** eofReached) {
	if( directStream ) return false;
	
	//else

	char x;
	in.get(x);

	if( readers.size() == 1 && eofReached ) {
		*eofReached = &in;
		directStream = true;
		if( in.eof() ) return false;
		else {
			in.putback(x);
			return true;
		}
	}

	if( in.eof() ) return false;

	int nread;

	if( x == delimiter ) {
		if( eofReached ) (*eofReached)->clear();
		for( typename RVector::iterator iter = readers.begin();
		     iter != readers.end(); ++iter )
			iter->second->put(delimiter);
		return true;
	}

	in.putback(x);
	in.get(buffer,BUFLIMIT-1,delimiter);
	nread = in.gcount();

	in.get(x);
	if( in.eof() && nread == 0 ) return false;

	if( eofReached ) (*eofReached)->clear();
	if( !in.eof() ) {
		buffer[nread] = x;
		buffer[++nread] = '\0';
	}

	for( typename RVector::iterator iter = readers.begin();
	     iter != readers.end(); ++iter )
		iter->second->write(buffer,nread);
	return true;
}

template<class Field>
const char* MatrixStream<Field>::getFormat() const {
	if( readers.empty() ) return "No valid format";
	else if( readers.size() > 1 ) return "Ambiguous Format";
	else return readers.begin()->first->getName();
}

template<class Field>
bool MatrixStream<Field>::isSparse() const {
	if( readers.size() == 0 ) return false;
	for( typename RVector::const_iterator iter = readers.begin();
	     iter != readers.end(); ++iter )
		if( !iter->first->isSparse() ) return false;
	return true;
}

template<class Field>
bool MatrixStream<Field>::isDense() const {
	if( readers.size() == 0 ) return false;
	for( typename RVector::const_iterator iter = readers.begin();
	     iter != readers.end(); ++iter )
		if( iter->first->isSparse() ) return false;
	return true;
}

} // end of namespace LinBox
