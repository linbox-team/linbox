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
#include <linbox/blackbox/random-matrix.h>
#include <linbox/blackbox/scompose.h>
#include <linbox/algorithms/matrix-rank.h>
#include <linbox/algorithms/last-invariant-factor.h>
// #include <linbox/algorithms/my-ith-invariant-factor.h>
// #include <linbox/algorithms/my-smith-form.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/algorithms/matrix-mod.h>
#include <linbox/solutions/rank.h>
#include <linbox/solutions/det.h>
#include <linbox/field/unparametric.h>
#include <linbox/field/modular-int32.h>
#include <linbox/field/modular.h>
#include <linbox/field/modular-short.h>
#include <linbox/field/modular-byte.h>
#include <linbox/field/ntl-ZZ.h>
#include <linbox/util/commentator.h>

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
	return out;
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

std::ostream& operator<<( std::ostream& out, const RingVectorBase& v ) {
	integer buf;
	v.reset();
	out << "DenseVector" << endl;
	while( v.getNext( buf ) )
		out << buf << endl;
	return out;
}

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

class IntegerVector
{
    protected:
    	bool dense;
	integer size;

    public:
    	virtual ~IntegerVector() {}
	bool isDense() { return dense; }
	const integer& getSize() { return size; }
};

template< class Ring >
class DenseIntegerVector :public IntegerVector, public vector<Ring>
{
    public: DenseIntegerVector( const vector<Ring> &v )
    	:vector<Ring>::vector(v) {
		dense = true;
		IntegerVector::size = (integer)v.size();
	}
};

template< class Ring >
class SparseIntegerVector :public IntegerVector,
                           public pair< vector<size_t>, vector<Ring> >
{
    public: SparseIntegerVector( const VectorVector<Ring>& v )
	    :pair<vector<size_t>,vector<Ring> >::pair(v[0].first,v[0].second)
	{
	    	dense = false;
		size = (integer) first.size();
	}
};

std::ostream* outPtr = &cout;

/*- @memo An abstract base class for integer matrices.
  * @doc This class has little functionality and is used just so I can treat
  * any matrix input the same way.
  */
class IntegerMatrix
{
    protected:
    	size_t rows, cols;

    public:
	virtual ~IntegerMatrix() {}

	integer getSize() const { return static_cast<integer>(rows) * cols; }
	size_t getRows() const { return rows; }
	size_t getCols() const { return cols; }

	virtual unsigned long rank() const = 0;
	virtual void smithForm(vector<NTL_ZZ::Element>&) const = 0;
	virtual integer& determinant(integer&) const = 0;

	int printRank() const {
		*outPtr << "Rank: " << rank() << endl;
		return 0;
	}
	int printSmithForm() const {
		vector<NTL_ZZ::Element> v;
		smithForm(v);
		*outPtr << "DenseVector" << endl;
		*outPtr << v;
		return 0;
	}
	int printDeterminant() const {
		integer res;
		*outPtr << "determinant: " << determinant(res) << endl;
		return 0;
	}
};

class MyRank
{
    public: unsigned long rank( const IntegerMatrix& m ) const
	{ return m.rank(); }
};

class MyRandom
{
    private:
    	integer max;

    public:
    	MyRandom( integer m ) :max(m) {}
	integer randomPrime() const {
		integer test;
		do integer::nonzerorandom(test,max);
		while( !mpz_probab_prime_p( test.get_rep(), 10 ) );
		return test;
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
		rows = m;
		cols = n;

		typename SparseMatrixBase<Field,Row>::RowIterator i;
		typename VectorVector<Ring>::const_iterator j;

		for(i = rowBegin(), j = vec.begin();
		    i != rowEnd() && j != vec.end();
		    ++i, ++j ) {
			i->first = j->first;
			i->second = j->second;
		}
	}

	integer& determinant(integer& i) const {
		SparseMatrixFactory<Modular<double>,
		                    Ring,
				    LinBox::Vector<Modular<double> >::Sparse,
				    Row>
			factory( *this );
		return LinBox::det( i, factory );
	}
	
	unsigned long rank() const {
	    UnparametricField<Ring> F;
	    unsigned long toRet;
	    if( currentRing >= INT16 ) {
		typedef Modular<Ring> ModField;
		MyRandom rand( ringMaxValues[ currentRing ] );
		ModField modField( rand.randomPrime() );
		SparseMatrix<ModField>* Ap;
		MatrixMod::mod< UnparametricField<Ring>, ModField >
			(Ap,*this,modField);
		LinBox::rank(toRet,*Ap,Ap->field());
	    }
	    else {
	    	typedef Modular<LinBox::int16> ModField;
		MyRandom rand( ringMaxValues[ INT16 ] );
		ModField modField( rand.randomPrime() );
		SparseMatrix<ModField>* Ap;
		MatrixMod::mod< UnparametricField<Ring>, ModField >
			(Ap,*this,modField);
	    	LinBox::rank(toRet,*Ap,Ap->field());
	    }
	    return toRet;
	}
	
	void smithForm(vector<NTL_ZZ::Element>& v) const {
		std::cerr << "No smith form of sparse integer matrix." << endl;
		exit(3);
	}
};

/*- @memo A simple extension of DenseMatrix
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
/*- @memo Constructor
  * @param r A vector containing all the elements of the matrix, in row order.
  * @param m Number of rows
  * @param n Number of columns
  * @param f (optional) Unparametric Field over the ring.
  */
	DenseIntegerMatrix( const vector<Ring> &r, size_t m, size_t n,
		  const UnparametricField<Ring> &f = UnparametricField<Ring>() )
		    :DenseMatrix< UnparametricField<Ring> >::DenseMatrix(f,m,n)
	{
		rows = m;
		cols = n;
		_rep = r;
		_ptr = &_rep[0];
	}

	integer& determinant(integer& i) const {
		DenseMatrixFactory< Modular<double>, Ring > factory( *this );
		return LinBox::det( i, factory );
	}

	unsigned long rank() const {
	    if( currentRing >= INT16 ) {
		MatrixRank<UnparametricField<Ring>,
		           Modular<Ring>,
			   MyRandom> 
			mr(_F,MyRandom(ringMaxValues[currentRing]));
		return mr.rank( *this );
	    }
	    else {
	    	MatrixRank<UnparametricField<Ring>,
		           Modular<LinBox::int16>,
			   MyRandom>
			mr(_F,MyRandom(ringMaxValues[INT16]));
		return mr.rank( *this );
	    }
	}

        void smithForm(vector<NTL_ZZ::Element>& v) const {
/*
                SmithForm<
                        NTL_ZZ,
                        IthInvariantFactor<
                                NTL_ZZ,
                                LastInvariantFactor<
                                        NTL_ZZ,
                                        RationalSolver<NTL_ZZ,
                                                       Modular<Ring>,
                                                       RandomPrime>
                                                   >,
                                SCompose,
                                RandomMatrix
                                          >,
                        MyRank
                         >
                        sf;
                sf.smithForm( v, *this );
*/
        }

};

enum RingType
	{ INT8, UINT8, INT16, UINT16, INT32, UINT32, INT64, UINT64, INTEGER };

integer ringMaxValues[] ={"127","255","32767","65535","2147483647","4294967295",
                          "9223372036854775807","18446744073709551615","0"};
integer ringMinValues[] = {"-128","0","-32768","0","-2147483648","0",
                          "-9223372036854775808","0","0"};
const int nRings = 9;


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

enum CommandType { RANK=0, SMITH_FORM, DETERMINANT, HELP };

struct Command
{
	char*		name;
	char*		alias;
	char*		helpString;
	CommandType	ct;
	bool		denseOnly;
};

struct Algorithm
{
	char*	name;
	char*	helpString;
	int	value;
};

enum { AMAX, AMIN, ASIZE, ACOUNTER, AROWS, ACOLS };

struct InputStorage {
	VectorVectorBase* sparseMatrix;
	VectorVectorBase* sparseVector;
	RingVectorBase* denseMatrix;
	RingVectorBase* denseVector;
};

// Some global variables
IntegerMatrix* matrixIn = NULL;
IntegerVector* vectorIn = NULL;
InputStorage ins = {NULL,NULL,NULL};
RingType currentRing = INT32;
int algorithm = 0;
bool denseOnly = false;
bool inputSwitching = true;
bool computerReadable = false;
bool showAlgorithms = false;

class RingBase {
    public:
	virtual void switchAndContinueDense( std::istream&, integer [],
				             RingVectorBase&, integer&,
					     bool ) = 0;
	virtual void switchAndContinueSparse( std::istream&, integer [],
                                              VectorVectorBase&,
					      unsigned long&, integer&,
					      bool ) = 0;
	virtual void switchAndContinueSparseAsDense( std::istream&,integer [],
	                                             RingVectorBase&,
						     unsigned long&,
						     unsigned long&,
						     integer&, bool ) = 0;
	virtual void readMatrix( std::istream&, bool ) = 0;
	virtual void readVector( std::istream&, bool ) = 0;
	virtual void switchAndSave( bool ) = 0;
};

template<class T> RingType continueReadingDense
	(std::istream&,integer [],RingVector<T>&,integer&,bool);
template<class T> RingType continueReadingSparse
	(std::istream&,integer[],VectorVector<T>&,unsigned long&,integer&,bool);
template<class T> RingType continueReadingSparseAsDense
	(std::istream&,integer[],RingVector<T>&,unsigned long&,unsigned long&,
	 integer&,bool);
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

	RingVector<Ring>* switchDense( const integer& size,
	                               RingVectorBase& oldVec ) {
		RingVector<Ring>* newVec = new RingVector<Ring>;
		newVec->reserve( (size_t) size );
		copyAndClear( *newVec, oldVec );
		delete &oldVec;
		return newVec;
	}

	VectorVector<Ring>* switchSparse( const integer& rows,
	                                  VectorVectorBase& oldVec ) {
		VectorVector<Ring>* newVec =
			new VectorVector<Ring>( (size_t) rows );
		typename VectorVector<Ring>::iterator i = newVec->begin();
		vector<size_t>* v1Ptr;
		RingVectorBase* v2Ptr;
		oldVec.reset();
		while( oldVec.getNext( v1Ptr, v2Ptr ) ) {
			i->first = *v1Ptr;
			copyAndClear( i->second, *v2Ptr );
			++i;
		}
		delete &oldVec;
		return newVec;
	}

	void switchAndSave( bool isMatrix ) {
	    if( isMatrix && ins.denseMatrix ) {
	    	RingVector<Ring>* newVec =
			switchDense(matrixIn->getSize(),*(ins.denseMatrix));
		ins.denseMatrix = newVec;
		IntegerMatrix* temp = new DenseIntegerMatrix<Ring>
			(*newVec,matrixIn->getRows(),matrixIn->getCols());
	    	delete matrixIn;
		matrixIn = temp;
	    }
	    else if( !isMatrix && ins.denseVector ) {
	    	RingVector<Ring>* newVec =
			switchDense(vectorIn->getSize(),*(ins.denseVector));
		ins.denseVector = newVec;
		delete vectorIn;
		vectorIn = new DenseIntegerVector<Ring>(*newVec);
	    }
	    else if( isMatrix && ins.sparseMatrix ) {
	    	VectorVector<Ring>* newVec =
			switchSparse((integer)matrixIn->getRows(),
			             *(ins.sparseMatrix));
	    	ins.sparseMatrix = newVec;
		IntegerMatrix* temp = new SparseIntegerMatrix<Ring>
			(*newVec,matrixIn->getRows(),matrixIn->getCols());
		delete matrixIn;
		matrixIn = temp;
	    }
	    else if( !isMatrix && ins.sparseVector ) {
	    	VectorVector<Ring>* newVec =
			switchSparse((integer)1,*(ins.sparseVector));
	    	ins.sparseVector = newVec;
		delete vectorIn;
		vectorIn = new SparseIntegerVector<Ring>(*newVec);
	    }
	}
	
	void switchAndContinueDense( std::istream& in, integer array[],
			             RingVectorBase& oldVec, integer& val,
				     bool isMatrix )
	{
	        RingVector<Ring>* newVec = switchDense(array[ASIZE],oldVec);
		newVec->push_back( (Ring) val );
        	RingType temp =
			continueReadingDense(in,array,*newVec,val,isMatrix);
        	if( temp != currentRing ) {
			currentRing = temp;
			if((!isMatrix && matrixIn) || (isMatrix && vectorIn) )
				(ringArray[currentRing])->
					switchAndSave(!isMatrix);
                	(ringArray[currentRing])->switchAndContinueDense
				(in,array,*newVec,val,isMatrix);
        	}
	}

	void switchAndContinueSparse( std::istream& in, integer array[],
				      VectorVectorBase& vec,
				      unsigned long& row, integer& val,
				      bool isMatrix ) {
		VectorVector<Ring>* newVec = switchSparse(array[AROWS],vec);
		(*newVec)[row].second.push_back( (Ring) val );
		RingType temp = continueReadingSparse
			(in,array,*newVec,row,val,isMatrix);
		if( temp != currentRing ) {
			currentRing = temp;
			if((!isMatrix && matrixIn) || (isMatrix && vectorIn) )
				(ringArray[currentRing])->
					switchAndSave(!isMatrix);
			(ringArray[currentRing])->switchAndContinueSparse
				(in,array,*newVec,row,val,isMatrix);
		}
	}

        void switchAndContinueSparseAsDense( std::istream& in, integer array[],
                                             RingVectorBase& oldVec,
                                             unsigned long& row,
					     unsigned long& col,
					     integer& val, bool isMatrix )
        {
                RingVector<Ring>* newVec = switchDense(array[ASIZE],oldVec);
                (*newVec)[(size_t)(row*array[ACOLS]+col)] = (Ring) val;
                RingType temp = continueReadingSparseAsDense
			(in,array,*newVec,row,col,val,isMatrix);
                if( temp != currentRing ) {
                        currentRing = temp;
                        if((!isMatrix && matrixIn) || (isMatrix && vectorIn) )
                                (ringArray[currentRing])->
                                        switchAndSave(!isMatrix);
                        (ringArray[currentRing])->switchAndContinueSparseAsDense
                                (in,array,*newVec,row,col,val,isMatrix);
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
			RingVector<Ring>* vecPtr = new RingVector<Ring>;
			vecPtr->reserve( (size_t) array[ASIZE] );
			if( (temp =
			       continueReadingDense(in,array,*vecPtr,val,true))
			    != currentRing ) {
				currentRing = temp;
				if( vectorIn )
				    (ringArray[currentRing])->
				    	switchAndSave(false);
				ringArray[currentRing]->switchAndContinueDense
					(in,array,*vecPtr,val,true);
			}
		}
		else {
		    unsigned long row;
		    if( denseOnly ) {
			unsigned long col;
		    	RingVector<Ring>* vecPtr =
			    new RingVector<Ring>( (size_t) array[ASIZE] );
			for( size_t i = 0; i < (size_t) array[ASIZE]; ++i )
				(*vecPtr)[i] = 0;
			if( (temp = continueReadingSparseAsDense
				(in,array,*vecPtr,row,col,val,true))
		    	    != currentRing ) {
				currentRing = temp;
				if( vectorIn )
				    (ringArray[currentRing])->
				    	switchAndSave(false);
				ringArray[currentRing]->
				    switchAndContinueSparseAsDense
				    	(in,array,*vecPtr,row,col,val,true);
			}
		    }
		    else {
			VectorVector<Ring>* vecPtr =
				new VectorVector<Ring>((size_t) array[AROWS] );
			if( (temp = continueReadingSparse
			 		(in,array,*vecPtr,row,val,true))
			    != currentRing ) {
				currentRing = temp;
				if( vectorIn )
				    (ringArray[currentRing])->
				    	switchAndSave(false);
				ringArray[currentRing]->switchAndContinueSparse
					(in,array,*vecPtr,row,val,true);
			}
		    }
		}
	}

	void readVector( std::istream& in, bool isDense ) {
		integer array[6];
		array[AMAX] = array[AMIN] = array[ACOUNTER] = 0;
		integer val;
		readInt( in, val );
		array[ASIZE] = array[ACOLS] = val;
		array[AROWS] = 1;
		RingType temp;
		if( isDense ) {
			RingVector<Ring>* vecPtr = new RingVector<Ring>;
			vecPtr->reserve( (size_t) array[ASIZE] );
			if( (temp =
			       continueReadingDense(in,array,*vecPtr,val,false))
			    != currentRing ) {
				currentRing = temp;
				if(matrixIn)
				    ringArray[currentRing]->
					switchAndSave(true);
				ringArray[currentRing]->switchAndContinueDense
					(in,array,*vecPtr,val,false);
			}
		}
		else {
		    unsigned long row;
		    if( denseOnly ) {
		    	unsigned long col;
                        RingVector<Ring>* vecPtr =
                        	new RingVector<Ring>( (size_t) array[ASIZE] );
                        for( size_t i = 0; i < (size_t) array[ASIZE]; ++i )
                                (*vecPtr)[i] = 0;
                        if( (temp = continueReadingSparseAsDense
				(in,array,*vecPtr,row,col,val,false))
			    != currentRing ) {
				currentRing = temp;
				if( matrixIn )
                                    (ringArray[currentRing])->
                                        switchAndSave(true);
                                ringArray[currentRing]->
				    switchAndContinueSparseAsDense
					(in,array,*vecPtr,row,col,val,false);
			}
		    }
		    else {
			VectorVector<Ring>* vecPtr = new VectorVector<Ring>(1);
			if( (temp = continueReadingSparse
					(in,array,*vecPtr,row,val,false))
			    != currentRing ) {
				currentRing = temp;
				if( matrixIn )
				    (ringArray[currentRing])->
				    	switchAndSave(true);
				ringArray[currentRing]->switchAndContinueSparse
					(in,array,*vecPtr,row,val,false);
			}
		    }
		}
	}
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
				unsigned long& rowInd, integer& val,
				bool isMatrix ) {
    if( currentRing == (RingType)(nRings-1) ) inputSwitching = false;
    unsigned long colInd;
    rowInd = 0;
    if( inputSwitching ) {
	int r = (int)currentRing;
	while( (isMatrix ? readInt(in,rowInd) : in)
	       && readInt(in,colInd) && readInt(in,val) ) {
		++array[ACOUNTER];
		if( rowInd >= array[AROWS] || colInd >= array[ACOLS] ) {
			std::cerr <<"ERROR: Dimension mismatch in matrix input."
				  << endl;
			exit(2);
		}
		vec[rowInd].first.push_back( colInd );
		if( val > array[AMAX] ) {
			array[AMAX] = val;
			if( val > ringMaxValues[r] ) break;
		}
		else if( val < array[AMIN] ) {
			array[AMIN] = val;
			if( val < ringMinValues[r] ) break;
		}
		if( array[ACOUNTER] == array[ASIZE] / 1000 ) {
			int i = 0;
			for( ;
			     i < nRings-1 &&
			     array[AMAX] > ringMaxValues[i] &&
			     array[AMIN] < ringMinValues[i];
			     ++i );
			if( i != r ) return (RingType)i;
		}
		vec[rowInd].second.push_back( (Ring) val );
	}
	if( val > ringMaxValues[r] ) {
		while( ++r < nRings-1 && val > ringMaxValues[r] );
		return (RingType)r;
	}
	else if( val < ringMinValues[r] ) {
		while( ++r < nRings-1 && val < ringMinValues[r] );
		return (RingType)r;
	}
    }
    else {
	Ring value;
	while( (isMatrix ? readInt(in,rowInd) : in )
	       && readInt(in,colInd) && readInt(in,value) ) {
                if( rowInd >= array[AROWS] || colInd >= array[ACOLS] ) {
                        std::cerr <<"ERROR: Dimension mismatch in matrix input."
                                  << endl;
                        exit(2);
                }
		vec[rowInd].first.push_back( colInd );
		vec[rowInd].second.push_back( value );
	}
    }
    if( isMatrix ) {
    	ins.sparseMatrix = &vec;
    	matrixIn = new SparseIntegerMatrix<Ring>
                           (vec, (size_t) array[AROWS], (size_t) array[ACOLS] );
    }
    else {
    	ins.sparseVector = &vec;
	vectorIn = new SparseIntegerVector<Ring>(vec);
    }
    return currentRing;
}

template< class Ring >
RingType continueReadingDense( std::istream& in, integer array[],
			       RingVector<Ring>& vec, integer& val,
			       bool isMatrix ) {
    if( currentRing == (RingType)(nRings-1) ) inputSwitching = false;
    if( inputSwitching ) {
	int r = (int)currentRing;
        while( readInt(in,val) ) {
                ++array[ACOUNTER];
                if( val > array[AMAX] ) {
			array[AMAX] = val;
       		        if( val > ringMaxValues[r] ) break;
               	}
               	else if( val < array[AMIN] ) {
			array[AMIN] = val;
                       	if( val < ringMinValues[r] ) break;
               	}
                if( array[ACOUNTER] == array[ASIZE] / 10 ) {
			int i = 0;
                        for( ;
                             i < nRings-1 &&
                             array[AMAX] > ringMaxValues[i] &&
                             array[AMIN] < ringMinValues[i];
                             i++ );
                        if( i != r ) return (RingType)i;
                }
		vec.push_back((Ring)val);
        }
        if( val > ringMaxValues[r] ) {
                while( ++r < nRings-1 && val > ringMaxValues[r] );
		return (RingType)r;
	}
        else if( val < ringMinValues[r] ) {
                while( ++r < nRings-1 && val < ringMinValues[r] );
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
    if( array[ACOUNTER] == array[ASIZE] ) {
	if( isMatrix ) {
    	    ins.denseMatrix = &vec;
            matrixIn = new DenseIntegerMatrix<Ring>
			   (vec,(size_t)array[AROWS],(size_t)array[ACOLS]);
    	}
	else {
	    ins.denseVector = &vec;
	    vectorIn = new DenseIntegerVector<Ring>( vec );
	}
    }
    else {
            std::cerr << "ERROR: Dimension mismatch in matrix input."
                      << endl;
            exit(2);
    }

    return currentRing;
}

template< class Ring >
RingType continueReadingSparseAsDense( std::istream& in, integer array[],
                                       RingVector<Ring>& vec,
                                       unsigned long& rowInd,
				       unsigned long& colInd,
				       integer& val, bool isMatrix ) {
    if( currentRing == (RingType)(nRings-1) ) inputSwitching = false;
    rowInd = 0;
    if( inputSwitching ) {
        int r = (int)currentRing;
        while( (isMatrix ? readInt(in,rowInd) : in)
               && readInt(in,colInd) && readInt(in,val) ) {
                ++array[ACOUNTER];
                if( rowInd >= array[AROWS] || colInd >= array[ACOLS] ) {
                        std::cerr <<"ERROR: Dimension mismatch in matrix input."
                                  << endl;
                        exit(2);
                }
                if( val > array[AMAX] ) {
                        array[AMAX] = val;
                        if( val > ringMaxValues[r] ) break;
                }
                else if( val < array[AMIN] ) {
                        array[AMIN] = val;
                        if( val < ringMinValues[r] ) break;
                }
                if( array[ACOUNTER] == array[ASIZE] / 1000 ) {
                        int i = 0;
                        for( ;
                             i < nRings-1 &&
                             array[AMAX] > ringMaxValues[i] &&
                             array[AMIN] < ringMinValues[i];
                             ++i );
                        if( i != r ) return (RingType)i;
                }
		vec[(size_t)(rowInd*array[ACOLS]+colInd)] = (Ring) val;
        }
        if( val > ringMaxValues[r] ) {
                while( ++r < nRings-1 && val > ringMaxValues[r] );
                return (RingType)r;
        }
        else if( val < ringMinValues[r] ) {
                while( ++r < nRings-1 && val < ringMinValues[r] );
                return (RingType)r;
        }
    }
    else {
        Ring value;
        while( (isMatrix ? readInt(in,rowInd) : in )
               && readInt(in,colInd) && readInt(in,value) ) {
                if( rowInd >= array[AROWS] || colInd >= array[ACOLS] ) {
                        std::cerr <<"ERROR: Dimension mismatch in matrix input."
                                  << endl;
                        exit(2);
                }
		vec[(size_t)(rowInd*array[ACOLS]+colInd)] = (Ring) val;
        }
    }
    if( isMatrix ) {
        ins.denseMatrix = &vec;
        matrixIn = new DenseIntegerMatrix<Ring>
                           (vec, (size_t) array[AROWS], (size_t) array[ACOLS] );
    }
    else {
        ins.denseVector = &vec;
        vectorIn = new DenseIntegerVector<Ring>(vec);
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
		else if( !strcasecmp( buf, "DenseVector" ) ) {
			(ringArray[currentRing])->readVector(in,true);
			failed = false;
		}
		else if( !strcasecmp( buf, "SparseVector" ) ) {
			(ringArray[currentRing])->readVector(in,false);
			failed = false;
		}
	}
	if( failed ) {
		std::cerr << "ERROR: File type invalid/not specified."
			  << endl;
		exit(2);
	}
}

void printAlgorithms( CommandType comm, Algorithm algs[][3] ) {
	int i, l;

	Algorithm* a = algs[(int)comm];

	if( computerReadable ) {
	    for(i = 0; a[i].name != NULL; ++i) {
	    	cout << a[i].name << endl;
		cout << a[i].helpString << endl;
	    }
	    return;
	}

	for (i = 0; a[i].name != NULL; i++) {
		cout << "  " << a[i].name;
		l = 20 - strlen(a[i].name);
		do cout << ' '; while (--l > 0);
		cout << a[i].helpString << endl;
	}
}

/* Display a help message on command usage */

void printHelpMessage (const char *program, Argument *args, Command *coms)
{
        int i, l;

	if( computerReadable ) {
	    for(i = 0; coms[i].name != NULL; ++i) {
	    	cout << coms[i].name << endl;
		cout << coms[i].helpString << endl;
	    }
	    return;
	}

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
                l = 15 - strlen (args[i].example);
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
		     Algorithm algs[][3], CommandType& comm)
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
				if(coms[i].denseOnly) denseOnly = true;
				break;
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

	if( showAlgorithms ) {
		printAlgorithms( comm, algs );
		comm = HELP;
		return;
	}
	if( comm == HELP ) {
		printHelpMessage(argv[0], args, coms);
		return;
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
	if( !matrixIn ) {
		std::cerr << "ERROR: No matrix input." << endl;
		exit(2);
	}
}

int main(int argc, char** argv) {
	RingSpecific<LinBox::int8> r0;
	RingSpecific<LinBox::uint8> r1;
	RingSpecific<LinBox::int16> r2;
	RingSpecific<LinBox::uint16> r3;
	RingSpecific<LinBox::int32> r4;
	RingSpecific<LinBox::uint32> r5;
	RingSpecific<LinBox::int64> r6;
	RingSpecific<LinBox::uint64> r7;
	RingSpecific<LinBox::integer> r8;

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
	char* algorithmName = NULL;
	integer min = 0, max = 0;

	static Argument args[] = {
		{'a',"-a name","Specify the algorithm to use.",TYPE_STRING,
		 &algorithmName},
		{'A',"-A",
		 "Show a list of available algorithms for the given command",
		 TYPE_ISNONE,&showAlgorithms},
		{'c',"-c",
		 "Output in computer-readable (not human-readable) format",
		 TYPE_ISNONE,&computerReadable},
		{'m',"-m minValue","Explicitly state minimum value in input",
		 TYPE_INTEGER,&min},
		{'M',"-M maxValue","Explicitly state maximum value in input",
		 TYPE_INTEGER,&max},
		{'o',"-o fileName","Output to a file (default is standard out)",
		 TYPE_STRING,&outFileName},
		{'\0',NULL,NULL,(ArgumentType)0,NULL}
	};

	static Command coms[] = {
		{"rank","r","Compute the rank of a given matrix.",RANK,false},
		{"smith-form","sf",
		 "Compute the Smith normal form of a given matrix.",SMITH_FORM,
		 true},
		{"determinant","d","Compute the determinant of a given matrix.",
		 DETERMINANT,false},
		{NULL,NULL,NULL,(CommandType)0,false}
	};

	static Algorithm algs[][3] = {
	    	{
		    {"elimination","Use Gaussian elimination and a randomly-chosen prime, Monte Carlo.",0},
		    {"wiedemann","Use Wiedemann method, again with a field constructed from a randomly-chosen prime.",1},
		    {NULL,NULL,0}
		},
		{
		    {"not working","Smith form is not working right now.",0},
		    {NULL,NULL,0}
		},
		{
		    {"chinese","Compute modulo a number of primes and reconstruct solution with Chinese remaindering.",0},
		    {NULL,NULL,0}
		},
	};

	CommandType comm;
	parseArguments( argc, argv, args, coms, algs, comm );

	if( computerReadable )
	    commentator.setBriefReportParameters
	    	(Commentator::OUTPUT_PIPE,true,true,true);

	if( algorithmName != NULL ) {
		int i = (int)comm;
		for( int j = 0; algs[i][j].name != NULL; ++j ) {
			if( !strcmp( algorithmName, algs[i][j].name ) ) {
				algorithm = algs[i][j].value;
				break;
			}
		}
	}

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
		     min < ringMinValues[i] &&
		     max > ringMaxValues[i];
		     ++i );
		currentRing = (RingType)i;
	}

	if( comm == HELP ) exit(0);
	int retVal = 0;

	switch( comm ) {
		case RANK: retVal = matrixIn->printRank(); break;
		case SMITH_FORM: retVal = matrixIn->printSmithForm(); break;
		case DETERMINANT: retVal = matrixIn->printDeterminant(); break;
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
