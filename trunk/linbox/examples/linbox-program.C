/* This is just a simple program to implement a number of LinBox functions
 * without the need to download the entire LinBox library.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <ctype.h>

using std::vector;
using std::pair;
using std::cout;
using std::endl;

#include <linbox/integer.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/algorithms/matrix-rank.h>
#include <linbox/algorithms/matrix-mod.h>
#include <linbox/solutions/rank.h>
#include <linbox/field/unparametric.h>
#include <linbox/field/modular-int32.h>
#include <linbox/field/modular.h>
#include <linbox/field/modular-short.h>
#include <linbox/field/modular-byte.h>

using namespace LinBox;

// For Debugging
template< class T >
std::ostream& operator<<( std::ostream& out, const vector<T>& v ) {
	typename vector<T>::const_iterator iter = v.begin();
	while( iter != v.end() ) {
		out << *iter << ' ';
		++iter;
	}
	out << endl;
}

class VectorBase {
    public:
	virtual ~VectorBase() {}
        virtual void doClear() = 0;
        virtual void reset() const = 0;
};

class RingVectorBase :public VectorBase{
    public:
        virtual bool getNext( integer& ) const = 0;
};

class VectorVectorBase :public VectorBase{
    public:
        virtual bool getNext( vector<size_t>*&, RingVectorBase*& ) = 0;
};

template< class T >
class GenericVector :public VectorBase, public vector<T> {
    protected:
        mutable typename vector<T>::const_iterator iter;
    public:
        GenericVector() :vector<T>::vector() { iter = begin(); }
        GenericVector(size_t n) :vector<T>::vector(n) { iter = begin(); }
        GenericVector(const GenericVector<T>& v) :vector<T>::vector(v)
                { iter = begin() + (v.iter - v.begin()); }

        void doClear() { clear(); }
        void reset() const { iter = begin(); }

        const T* next() const {
                if( iter == end() ) return NULL;
                return &(*(iter++));
        }
};

template< class Ring >
class RingVector :public GenericVector<Ring>, public RingVectorBase {
    public:
	RingVector() :GenericVector<Ring>::GenericVector() {}
        RingVector( size_t n ) :GenericVector<Ring>::GenericVector(n) {}

	void doClear() { GenericVector<Ring>::doClear(); }
	void reset() const { GenericVector<Ring>::reset(); }

        bool getNext( integer& a ) const {
                const Ring* ptr = next();
                if( ptr ) a = *ptr;
                return (bool)ptr;
        }
};

template< class Ring >
class VectorVector
        :public GenericVector<pair<vector<size_t>,RingVector<Ring> > >,
         public VectorVectorBase
{
    public:
	typedef pair<vector<size_t>,RingVector<Ring> > Element;

	VectorVector() :GenericVector<Element>::GenericVector() {}
        VectorVector(size_t n) :GenericVector<Element>::GenericVector(n) {}

	void doClear() { GenericVector<Element>::doClear(); }
	void reset() const { GenericVector<Element>::reset(); }

        bool getNext( vector<size_t>* &v1, RingVectorBase* &v2 ) {
                Element* ptr = const_cast<Element*>( next() );
                if( !ptr ) return false;
                v1 = &(ptr->first);
                v2 = &(ptr->second);
                return true;
        }
};

//For Debugging
template< class Ring >
std::ostream& operator<<( std::ostream& out, const VectorVector<Ring>& v ) {
	out << "VectorVector:" << endl;
	size_t rowInd;
	vector<size_t>::const_iterator colIter;
	typename RingVector<Ring>::const_iterator valIter;
	for( rowInd = 0; rowInd < v.size(); ++rowInd )
		for( colIter = v[rowInd].first.begin(),
		     valIter = v[rowInd].second.begin();
		     colIter != v[rowInd].first.end(),
		     valIter != v[rowInd].second.end();
		     ++colIter, ++valIter )
			out << rowInd << ' '
			    << *colIter << ' '
			    << *valIter << endl;
	out << endl;
	return out;
}

template< class Ring >
std::ostream& operator<<( std::ostream& out, const RingVector<Ring>& v ) {
	return out << (vector<Ring>)v;
}


std::ostream* outPtr = &cout;

/** @memo An abstract base class for integer matrices.
  * @doc This class has little functionality and is used just so I can treat
  * any matrix input the same way.
  */
class IntegerMatrix
{
    public:
	virtual ~IntegerMatrix() {}

	virtual unsigned long rank() = 0;

	int printRank() {
		*outPtr << rank() << std:: endl;
		return 0;
	}
};

template< class Ring >
class SparseIntegerMatrix : public IntegerMatrix,
			    public SparseMatrix< UnparametricField<Ring>,
			                         pair< vector<size_t>,
                                                       vector<Ring>
                                                     >
                                               >
{
    public:
	typedef pair< vector<size_t>, vector<Ring> > Row;
	typedef UnparametricField<Ring> Field;

	SparseIntegerMatrix( const VectorVector<Ring>& vec, size_t m, size_t n,
		             const Field &f = Field() )
		:SparseMatrix<Field,Row>::SparseMatrix(f,m,n)
	{
		typename SparseMatrixBase<Field,Row>::RowIterator i;
		typename VectorVector<Ring>::const_iterator j;

		for(i = rowBegin(), j = vec.begin();
		    i != rowEnd() && j != vec.end();
		    ++i, ++j ) {
			i->first = j->first;
			i->second = j->second;
		}
	}

	unsigned long rank() {
		UnparametricField<Ring> F;
		MatrixRank<UnparametricField<Ring>,Modular<Ring> > mr(F);
		Modular<Ring> ModField( mr.rp.randomPrime() );
		SparseMatrix<Modular<Ring> >* Ap;
		MatrixMod::mod< UnparametricField<Ring>, Modular<Ring> >
			(Ap,*this,ModField);
		unsigned long toRet;
		LinBox::rank(toRet,*Ap,Ap->field());
		return toRet;
	}
};

/** @memo A simple extension of DenseMatrix
  * @doc This is just a simple class that inherits from DenseMatrix (found in
  * linbox/blackbox/dense.h), and from the abstract IntegerMatrix class defined
  * above.  It is templatized over a ring (which will be used to make an
  * UnparametricField).  Also, it allows initialization from an element vector,
  * which is the way the data is stored internally anyway.
  */
template< class Ring >
class DenseIntegerMatrix : public IntegerMatrix,
			   public DenseMatrix< UnparametricField<Ring> >
{
    public:
/** @memo Constructor
  * @param r A vector containing all the elements of the matrix, in row order.
  * @param m Number of rows
  * @param n Number of columns
  * @param f (optional) Unparametric Field over the ring.
  */
	DenseIntegerMatrix( const vector<Ring> &r, size_t m, size_t n,
		  const UnparametricField<Ring> &f = UnparametricField<Ring>() )
		    :DenseMatrix< UnparametricField<Ring> >::DenseMatrix(f,m,n)
	{
		_rep = r;
		_ptr = &_rep[0];
	}

	unsigned long rank() {
		MatrixRank<UnparametricField<Ring>,Modular<Ring> > mr(_F);
		return mr.rank( *this );
	}
};

//Most the following code about arguments is adapted from tests/test-common.C

enum ArgumentType {
        TYPE_ISNONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE, TYPE_STRING
};

struct Argument
{
        char             c;
        char            *example;
        char            *helpString;
        ArgumentType     type;
        void            *data;
};

enum CommandType { HELP, RANK, SMITH_FORM, DETERMINANT, MINPOLY };

struct Command
{
	char*		name;
	char*		alias;
	char*		helpString;
	CommandType	ct;
};

enum RingType
	{ INT8, UINT8, INT16, UINT16, INT32, UINT32, INT64, UINT64, INTEGER };

integer ringMaxSizes[] = {"127","255","32767","65535","2147483647","4294967295",
			  "9223372036854775807","18446744073709551615","0"};
integer ringMinSizes[] = {"-128","0","-32768","0","-2147483648","0",
			  "-9223372036854775808","0","0"};
const int nRings = 9;

enum { AMAX, AMIN, AROWS, ACOLS, ASIZE, ACOUNTER };

// Some global variables
IntegerMatrix* matrixIn = NULL;
void* vectorIn = NULL;
bool vectorIsSparse;
RingType currentRing = INT32;
bool inputSwitching = true;

class RingBase {
    public:
	virtual void switchAndContinueDense( std::istream&, integer [],
				             RingVectorBase&, integer& ) = 0;
	virtual void switchAndContinueSparse( std::istream&, integer [],
                                              VectorVectorBase&,
					      unsigned long&, integer& ) = 0;
	virtual void readMatrix( std::istream&, bool ) = 0;
};

template<class T> RingType continueReadingDense
	(std::istream&,integer [],RingVector<T>&,integer&);
template<class T> RingType continueReadingSparse
	(std::istream&,integer [],VectorVector<T>&,unsigned long&,integer&);
template<class T> inline std::istream& readInt( std::istream&, T& );
RingBase* ringArray[nRings];

template< class Ring >
class RingSpecific :public RingBase {
    public:
	void copyAndClear( RingVector<Ring>& newVec, RingVectorBase& oldVec ){
		integer buf;
		oldVec.reset();
		while( oldVec.getNext( buf ) )
			newVec.push_back( (Ring) buf );
		oldVec.doClear();
	}

	void switchAndContinueDense( std::istream& in, integer array[],
			             RingVectorBase& vec, integer& val)
	{
	        RingVector<Ring> newVec;
		newVec.reserve( (size_t) array[ASIZE] );
		copyAndClear( newVec, vec );
		newVec.push_back( (Ring) val );
        	RingType temp = continueReadingDense(in,array,newVec,val);
        	if( temp != currentRing ) {
                	currentRing = temp;
                	(ringArray[currentRing])->
				switchAndContinueDense(in,array,newVec,val);
        	}
	}

	void switchAndContinueSparse( std::istream& in, integer array[],
				      VectorVectorBase& vec,
				      unsigned long& row, integer& val) {
		VectorVector<Ring> newVec( (size_t) array[AROWS] );
		typename VectorVector<Ring>::iterator i = newVec.begin();
		vector<size_t>* v1Ptr;
		RingVectorBase* v2Ptr;
		vec.reset();
		while( vec.getNext( v1Ptr, v2Ptr ) ) {
			i->first = *v1Ptr;
			copyAndClear( i->second, *v2Ptr );
			++i;
		}
		vec.doClear();
		newVec[row].second.push_back( (Ring) val );
		RingType temp = continueReadingSparse(in,array,newVec,row,val);
		if( temp != currentRing ) {
			currentRing = temp;
			(ringArray[currentRing])->switchAndContinueSparse
				(in,array,newVec,row,val);
		}
	}

	void readMatrix( std::istream& in, bool isDense ) {
        	integer array[6];
        	array[AMAX] = 0;
        	array[AMIN] = 0;
        	integer val;
        	readInt(in,val);
        	array[AROWS] = val;
        	readInt(in,val);
        	array[ACOLS] = val;
        	array[ASIZE] = array[AROWS] * array[ACOLS];
        	array[ACOUNTER] = 0;
		RingType temp;
		if( isDense ) {
			RingVector<Ring> vec;
			vec.reserve( (size_t) array[ASIZE] );
			if( (temp = continueReadingDense(in,array,vec,val))
			    != currentRing ) {
				currentRing = temp;
				ringArray[currentRing]->switchAndContinueDense
					(in,array,vec,val);
			}
		}
		else {
			VectorVector<Ring> vec((size_t) array[AROWS] );
			unsigned long row;
			if( (temp = continueReadingSparse(in,array,vec,row,val))
			    != currentRing ) {
				currentRing = temp;
				ringArray[currentRing]->switchAndContinueSparse
					(in,array,vec,row,val);
			}
		}
	}
};

template< class Ring >
struct MatrixTypes {
        typedef DenseIntegerMatrix<Ring> Dense;
};

template <class T>
inline std::istream& readInt( std::istream& in, T& i ) {
	char temp;
	while( !in.eof() && in.get(temp) && !isdigit(temp) );
	if( in.eof() ) in.setstate( std::ios_base::failbit );
	else {
		in.putback(temp);
		in >> i;
	}
	return in;
}

template< class Ring >
RingType continueReadingSparse( std::istream& in, integer array[],
				VectorVector<Ring>& vec,
				unsigned long& rowInd, integer& val ) {
    if( currentRing == (RingType)(nRings-1) ) inputSwitching = false;
    unsigned long colInd;
    if( inputSwitching ) {
	int r = (int)currentRing;
	while( readInt(in,rowInd) && readInt(in,colInd) && readInt(in,val) ) {
		++array[ACOUNTER];
		if( rowInd >= array[AROWS] || colInd >= array[ACOLS] ) {
			std::cerr <<"ERROR: Dimension mismatch in matrix input."
				  << endl;
			exit(2);
		}
		vec[rowInd].first.push_back( colInd );
		if( val > array[AMAX] ) {
			array[AMAX] = val;
			if( val > ringMaxSizes[r] ) break;
		}
		else if( val < array[AMIN] ) {
			array[AMIN] = val;
			if( val < ringMinSizes[r] ) break;
		}
		if( array[ACOUNTER] == array[ASIZE] / 1000 ) {
			int i = 0;
			for( ;
			     i < nRings-1 &&
			     array[AMAX] > ringMaxSizes[i] &&
			     array[AMIN] < ringMinSizes[i];
			     ++i );
			if( i != r ) return (RingType)i;
		}
		vec[rowInd].second.push_back( (Ring) val );
	}
	if( val > ringMaxSizes[r] ) {
		while( ++r < nRings-1 && val > ringMaxSizes[r] );
		return (RingType)r;
	}
	else if( val < ringMinSizes[r] ) {
		while( ++r < nRings-1 && val < ringMinSizes[r] );
		return (RingType)r;
	}
    }
    else {
	Ring value;
	while( readInt(in,rowInd) && readInt(in,colInd) && readInt(in,value) ) {
                if( rowInd >= array[AROWS] || colInd >= array[ACOLS] ) {
                        std::cerr <<"ERROR: Dimension mismatch in matrix input."
                                  << endl;
                        exit(2);
                }
		vec[rowInd].first.push_back( colInd );
		vec[rowInd].second.push_back( value );
	}
    }
    matrixIn = new SparseIntegerMatrix<Ring>
                           (vec, (size_t) array[AROWS], (size_t) array[ACOLS] );
    return currentRing;
}

template< class Ring >
RingType continueReadingDense( std::istream& in, integer array[],
			       RingVector<Ring>& vec, integer& val ) {
    if( currentRing == (RingType)(nRings-1) ) inputSwitching = false;
    if( inputSwitching ) {
	int r = (int)currentRing;
        while( readInt(in,val) ) {
                ++array[ACOUNTER];
                if( val > array[AMAX] ) {
			array[AMAX] = val;
       		        if( val > ringMaxSizes[r] ) break;
               	}
               	else if( val < array[AMIN] ) {
			array[AMIN] = val;
                       	if( val < ringMinSizes[r] ) break;
               	}
                if( array[ACOUNTER] == array[ASIZE] / 10 ) {
			int i = 0;
                        for( ;
                             i < nRings-1 &&
                             array[AMAX] > ringMaxSizes[i] &&
                             array[AMIN] < ringMinSizes[i];
                             i++ );
                        if( i != r ) return (RingType)i;
                }
		vec.push_back((Ring)val);
        }
        if( val > ringMaxSizes[r] ) {
                while( ++r < nRings-1 && val > ringMaxSizes[r] );
		return (RingType)r;
	}
        else if( val < ringMinSizes[r] ) {
                while( ++r < nRings-1 && val < ringMinSizes[r] );
		return (RingType)r;
	}
    }
    else {
	Ring value;
	while( readInt(in,value) ) {
		++array[ACOUNTER];
		vec.push_back(value);
	}
    }
    if( array[ACOUNTER] == array[ASIZE] )
            matrixIn = new DenseIntegerMatrix<Ring>
			   (vec,(size_t)array[AROWS],(size_t)array[ACOLS]);
    else {
            std::cerr << "ERROR: Dimension mismatch in matrix input."
                      << endl;
            exit(2);
    }

    return currentRing;
}

void readFile( std::istream& in ) {
	char buf[15];
	bool failed = true;
	while( failed && in >> std::setw(15) >> buf ) {
		if( !strcasecmp( buf, "DenseMatrix" ) ) {
			(ringArray[currentRing])->readMatrix(in,true);
			failed = false;
		}
		else if( !strcasecmp( buf, "SparseMatrix" ) ) {
			(ringArray[currentRing])->readMatrix(in,false);
			failed = false;
		}
	}
	if( failed ) {
		std::cerr << "ERROR: File type invalid/not specified."
			  << endl;
		exit(2);
	}
}

/* Display a help message on command usage */

void printHelpMessage (const char *program, Argument *args, Command *coms)
{
        int i, l;

        cout << "Usage: " << program << " command [options] [inputFile...]"
	     << endl;
        cout << endl;
	cout << "Where command is one of the following:" << endl;
	for (i = 0; coms[i].name != NULL; i++) {
		cout << "  " << coms[i].name << "(" << coms[i].alias << ")";
		l = 20 - strlen(coms[i].name) - strlen(coms[i].alias);
		do cout << ' '; while (--l > 0);
		cout << coms[i].helpString << endl;
	}

        cout << endl << "[options] are the following:" << endl;

        for (i = 0; args[i].c != '\0'; i++) {
                cout << "  " << args[i].example;
                l = 10 - strlen (args[i].example);
                do cout << ' '; while (--l > 0);
                cout << args[i].helpString << endl;
        }


        cout   << endl
	       << "inputFiles can be in any one of the following formats:"
	       << endl << endl;

	cout   << "  Dense Matrix:"
	       << endl << "     "
	       << "The word " << '"'
	       << "DenseMatrix" << '"'
	       << " (case insensitive) must appear before"
	       << endl << "     "
	       << "any  matrix  data.  After this, all  non-digits are ignored."
	       << endl << "     "
	       << "The first two numbers must be the row and column dimensions,"
	       << endl << "     "
	       << "respectively.  After this,  each entry of the  matrix should"
	       << endl << "     "
	       << "given, in row order  (all of first row,  then all of  second"
	       << endl << "     "
	       << "row, etc.)."
	       << endl << endl;

        cout   << "  Sparse Matrix:"
               << endl << "     "
               << "The  word  " << '"'
               << "SparseMatrix" << '"'
               << "  (case  insensitive)  must  appear"
               << endl << "     "
               << "before  any  matrix  data.  After  this, all  non-digits are"
               << endl << "     "
               << "ignored.  The first two  numbers must be the  row and column"
               << endl << "     "
               << "dimensions,  respectively.  After this, numbers will be read"
	       << endl << "     "
	       << "in  triples,  indicating the  row index,  column index,  and"
	       << endl << "     "
	       << "value  of an entry,  respectively.  The triples may be given"
	       << endl << "     "
	       << "in any order."
               << endl << endl;

        cout   << "  Dense Vector:"
               << endl << "     "
               << "The word " << '"'
               << "DenseVector" << '"'
               << " (case insensitive) must appear before"
               << endl << "     "
               << "any  vector  data.  After this, all  non-digits are ignored."
               << endl << "     "
               << "The first  number  must be  the size  of the  vector.  After"
               << endl << "     "
               << "this, each entry of the vector should be given, in order."
               << endl << endl;

        cout   << "  Sparse Vector:"
               << endl << "     "
               << "The  word  " << '"'
               << "SparseVector" << '"'
               << "  (case  insensitive)  must  appear"
               << endl << "     "
               << "before  any  matrix  data.  After  this, all  non-digits are"
               << endl << "     "
               << "ignored.  The first  number must be the  size of the vector."
               << endl << "     "
	       << "After this, each numbers will be read in doubles, indicating"
	       << endl << "     "
	       << "the index and value of an entry, respectively.  The doubles"
	       << endl << "     "
	       << "may be given in any order."
               << endl;
}

Argument *findArgument (Argument *args, char c)
{
        int i;

        for (i = 0; args[i].c != '\0' && args[i].c != c; i++);

        if (args[i].c != '\0')
                return &(args[i]);
        else
                return (Argument *) 0;
}

/* Parse command line arguments */

void parseArguments (int argc, char **argv, Argument *args, Command* coms,
		     CommandType& comm)
{
        int i=-1;
        Argument *current;
	comm = HELP;

	if( argc >= 2 &&
	    strcmp( argv[1], "h" ) &&
	    strcmp( argv[1], "?" ) &&
	    strcmp( argv[1], "help" ))
		while( coms[++i].name != NULL )
			if( !(strcmp( argv[1], coms[i].name ) &&
			      strcmp( argv[1], coms[i].alias ))) {
				comm = coms[i].ct;
				break;
			}

	if( comm == HELP ) {
		printHelpMessage(argv[0], args, coms);
		return;
	}

        for (i = 2; i < argc && argv[i][0] == '-'; i++) {
                if ( (current = findArgument (args, argv[i][1])) 
		    != (Argument *) 0 ) {
                        switch (current->type) {
                        case TYPE_ISNONE:
                                *(bool *) current->data = true;
                                break;

                        case TYPE_INT:
                                *(int *) current->data = atoi (argv[i+1]);
                                i++;
                                break;

                        case TYPE_INTEGER:
                                *(integer *) current->data = integer(argv[i+1]);
                                i++;
                                break;

                        case TYPE_DOUBLE:
                                *(double *) current->data = atof (argv[i+1]);
                                i++;
                                break;
			case TYPE_STRING:
				*(char**) current->data = 
					new char[strlen(argv[i+1])+1];
				strcpy( *(char**) current->data, argv[i+1]);
				i++;
				break;
                        }
                } else {
                        std::cerr << "ERROR: Bad argument " << argv[i] << endl;
                        exit(1);
                }
	}
	if( i == argc ) readFile( std::cin );
	else for( std::ifstream fin; i < argc; ++i ) {
		fin.open( argv[i] );
		if( !fin ) {
			std::cerr << "ERROR: Could not open file " << argv[i]
				  << endl;
			exit(1);
		}
		readFile( fin );
		fin.close();
	}
	if( comm != HELP && !matrixIn ) {
		std::cerr << "ERROR: No matrix input." << endl;
		exit(2);
	}
}

int main(int argc, char** argv) {
	RingSpecific<int8> r0;
	RingSpecific<uint8> r1;
	RingSpecific<int16> r2;
	RingSpecific<uint16> r3;
	RingSpecific<int32> r4;
	RingSpecific<uint32> r5;
	RingSpecific<int64> r6;
	RingSpecific<uint64> r7;
	RingSpecific<integer> r8;

	ringArray[INT8] = &r0;
	ringArray[UINT8] = &r1;
	ringArray[INT16] = &r2;
	ringArray[UINT16] = &r3;
	ringArray[INT32] = &r4;
	ringArray[UINT32] = &r5;
	ringArray[INT64] = &r6;
	ringArray[UINT64] = &r7;
	ringArray[INTEGER] = &r8;

	char* outFileName = NULL;
	integer min = 0, max = 0;

	static Argument args[] = {
		{'m',"-m minValue","Explicitly state minimum value in input",
		 TYPE_INTEGER,&min},
		{'M',"-M maxValue","Explicitly state maximum value in input",
		 TYPE_INTEGER,&max},
		{'o',"-o fileName","Output to a file (default is standard out)",
		 TYPE_STRING,&outFileName},
	};
	static Command coms[] = {
		{"rank","r","Compute the rank of a given matrix.",RANK},
		{"smith-form","sf",
		 "Compute the Smith normal form of a given matrix.",SMITH_FORM},
	};

	CommandType comm;
	parseArguments( argc, argv, args, coms, comm );

	std::ofstream* fout = NULL;
	if(outFileName) {
		fout = new std::ofstream(outFileName);
		if( fout ) outPtr = fout;
		else std::cerr << "Error opening output file " << outFileName
			       << ".  Standard out will be used." << endl;
	}

	if( min < max ) {
		inputSwitching = false;
		int i = 0;
		for( ;
		     i < nRings &&
		     min < ringMinSizes[i] &&
		     max > ringMaxSizes[i];
		     ++i );
		currentRing = (RingType)i;
	}

	if( comm == HELP ) exit(0);
	int retVal = 0;

	switch( comm ) {
		case RANK: retVal = matrixIn->printRank(); break;
		default: cout << "I don't know how to do that yet." << endl;
	}
	if( matrixIn ) delete matrixIn;
	if( vectorIn ) delete vectorIn;
	if( fout ) {
		fout->close();
		delete fout;
	}

	return retVal;
}
