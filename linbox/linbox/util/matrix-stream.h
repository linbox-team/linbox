/* matrix-stream.h
 * Take a file or stream of a matrix, return (row, column, value) triples from
 * that matrix.  Automatically determine matrix format.
 * 
 * NOTE to developers:  To add a new format, write a subclass of
 *                      MatrixStreamReader and put it in the formats
 *                      subdirectory.  Then you must add this new format to
 *                      the file formats/matrix-stream-readers.h in two
 *                      different places.  It should be clear how to do this if
 *                      you look at that file.  Once you do this, the new format
 *                      will automatically be used.
 *
 * by Dan Roche, 1-12-05
 */

#ifndef __MATRIX_STREAM_H
#define __MATRIX_STREAM_H

#include <iostream>
#include <vector>
#include <sstream>
#include <queue>

namespace LinBox { // namespace in which all LinBox code resides

//  Error codes enumeration
enum MatrixStreamError {
	AMBIGUOUS_FORMAT=-1, // It is not yet clear which format this file is in
	GOOD, // There is one format that has not yet encountered any errors
	END_OF_MATRIX, // The end of tha matrix has been reached; no more data
	END_OF_FILE, // An EOF was encountered unexpectedly in reading the data
	BAD_FORMAT, // This indicates that the format is recognized, but the
	            // file itself is not properly formatted.  For instance,
		    // indices out of bounds.
	NO_FORMAT  // Indicates that the format of the given file is not
	           // recognized by any of the readers.
};

template <class Field> class MatrixStream;

/** An abstract base class to represent readers for specific formats.  For each
 * format that is to be supported, make an extension of this class that
 * implements protected methods nextTripleImpl and initImpl.
 \ingroup util
 */
template <class Field>
class MatrixStreamReader {
    public:
    	typedef typename Field::Element Element;
    private:
    /** This is a function that performs the work for the protected
     * readWhiteSpace et al methods below.
     * @param allowBreaks If false, function will return true only if no
     *                    newline characters are encountered.
     * @param some If false, function will return true only if at least one
     *             whitespace character is read.
     * @param breaks Pointer to an initialized int.  If stopAfterBreaks is
     *               false, then this int will be returned with the number of
     *               breaks encountered.  If stopAfterBreaks is true, then this
     *               should be initialized to the number of breaks after which
     *               to stop.
     * @param stopAfterBreaks If true, indicates that reading should stop after
     *                        *breaks breaks are encountered.
     */
    	bool genericWSReader( bool allowBreaks, bool some, int* breaks = NULL,
	                      bool stopAfterBreaks = false );
    /** The std::queue used to save triples that have been read but have not yet been
     * requested with a call to nextTriple
     */
	std::queue<std::pair<std::pair<size_t,size_t>,Element> > savedTriples;

    /** The last error returned by saveNext(), to be returned by nextTriple
     * after savedTriples has been emptied.
     */
     	MatrixStreamError lastError;

    /* The protected members and methods are provided to save on implementation
     * time for making the readers. */

    protected:

    /** The stream that provides the data to the reader.
     * NOTE: subclasses should NOT use this stream directly except for one-byte
     * reading as in sin->get().  This stream only contains a portion of the
     * matrix data; this data must be replentished with calls to moreData().
     * If at all possible use sin->get() only and use the various read...
     * methods below to read data.
     */
    	std::istream* sin;

    /** A pointer to the MatrixStream that is using this reader.  Useful to
      * get an instance of the field via ms->getField().
      */
	MatrixStream<Field>* ms;
	
    /** The lineNumber is recorded in case the user wants to know at which line
      * an error occurs.  This will be updated automatically by any of the read
      * methods below if they encounter breaks; it is up to the subclasses to 
      * increment lineNumber if they read any newline characters independently.
      */
	int lineNumber;

    /** The number of rows in the matrix.  This will be set by default to 0.
      */
	size_t _m;

    /** Indicates whether the number above is accurate */
    	bool knowM;
	
    /** Number of columns in the matrix.  Similar requirements as _m above. */
	size_t _n;

    /** Indicates whether the number above is accurate */
    	bool knowN;

    /** Indicates that the end of the matrix has been reached; no more calls to
      * nextTripleImpl will be made once this value is true.  This will
      * automatically be set to true if nextTripleImple returns END_OF_MATRIX.
      */
	bool atEnd;

    /** Read white space.  Function returns true if and only if at least one
     * character of white space is read.  After a successful call, there will
     * be at least one character available on the stream.
     * @param allowBreaks If false, the function will return false if it
     *                    encounters a newline.
     * @throws END_OF_FILE if an EOF is encountered before any non-whitespace
     *                     characters.
     */
	bool readWhiteSpace( bool allowBreaks = false )
		{ return genericWSReader( allowBreaks, false ); }
	
    /** Read white space.  Does not require that any white space characters be
     * read.  After a successful call, there will be at least one character
     * available on the stream.
     * @param allowBreaks If false, the function will return false if it
     *                    encounters a newline.
     * @throws END_OF_FILE if an EOF is encountered before any non-whitespace
     *                     characters.
     */
	bool readSomeWhiteSpace( bool allowBreaks = false )
		{ return genericWSReader( allowBreaks, true ); }

    /** Read white space.  Function returns true if and only if at least one
     * character of white space is read.  After a successful call, there will
     * be at least one character available on the stream.
     * @param breaks Upon successful return, this int will contain the number
     *               of breaks (newlines) that were read.
     * @param allowBreaks If false, the function will return false if it
     *                    encounters a newline.
     * @throws END_OF_FILE if an EOF is encountered before any non-whitespace
     *                     characters.
     */
	bool readWhiteSpace( int& breaks, bool allowBreaks = true )
		{ return genericWSReader( allowBreaks, false, &breaks ); }

    /** Read white space.  Does not require that any white space characters be
     * read.  After a successful call, there will be at least one character
     * available on the stream.
     * @param breaks Upon successful return, this int will contain the number
     *               of breaks (newlines) that were read.
     * @param allowBreaks If false, the function will return false if it
     *                    encounters a newline.
     * @throws END_OF_FILE if an EOF is encountered before any non-whitespace
     *                     characters.
     */
	bool readSomeWhiteSpace( int& breaks, bool allowBreaks = true )
		{ return genericWSReader( allowBreaks, true, &breaks ); }

    /** Read up to breaks breaks.  Reading will stop on the first non-whitespace
     * character or first newline after breaks newlines.  After a successful
     * call, there will be at least one character available on the stream.
     * @param breaks The number of breaks to read.
     * @throws END_OF_FILE if an EOF is encountered before any non-whitespace
     *                     characters.
     */
	bool readBreaks(int breaks = 1)
		{ return genericWSReader( true, false, &breaks, true ); }

    /** Read up to a given character.
     * @param c The character to read until.  This character is extracted and
     *          not put back.
     * @param ss A pointer to a std::stringstream to which to write all characters
     *           read (except, again, for c)
     * @param limit The maximum number of characters to read.
     * @return True if c is read before limit characters are seen, false o.w.
     */
	bool readUntil(char c, std::stringstream* ss = NULL, int limit = -1);

    /** Read until an unmatched character.
     * @param cm A std::vector of char* strings, each containing (at least) two
     *           characters: the left-match and the right-match, e.g. "()".
     *           If the character has no match (e.g. ','), the first character
     *           in the array should be null, e.g. "\0,".
     * @param ss A pointer to a std::stringstream to which to write all characters
     *           read (except for the unmatched one)
     * @param limit The maximum number of characters to read.
     * @return The iterator that points to the std::pair in cm that was unmatched.
     */
	std::vector<char*>::const_iterator
		readUntil(const std::vector<char*>& cm, std::stringstream* ss = NULL,
		          int limit = -1 );

    /** Read a field element.  Uses the read method of the field for the parent
     * MatrixStream object.
     * @param x Where to store the read element.
     * @return True iff the read was successful
     */
	bool readElement( Element& x );

    /** Read any object.  Uses the overloaded stream extraction operator >>,
     * which must be defined for this type.
     * @param o Where to store the read object
     * @return True iff the read was successful.
     */
	template<class Object> bool readObject( Object& o );

    /** Try and get more data from the underlying input stream.
      * Should be called when an EOF is reached on input.
      * @return true iff more data was available and has been added to sin.
      */
	bool moreData();
	
    /** Save the triple (m,n,v) onto the savedTriples std::queue. */
	void saveTriple(size_t m, size_t n, const Element& v);
	
    /** Read the next triple of row index, column index, value and store it in
     * the given references.
     * @return A MatrixStreamError indicating the success or failure of the
     *         operation
     */
	virtual MatrixStreamError nextTripleImpl(size_t&,size_t&,Element&) = 0;
	
    /** Read the beginning (header) of the matrix from the stream and attempt to
     * determine if it is of this reader's type.
     * @return A MatrixStreamError indicating the success or failure of the
     *         operation
     */
	virtual MatrixStreamError initImpl() = 0;

    /** A protected constructor that is called automatically when subclasses
     * are instantiated. */
	MatrixStreamReader() {
		sin = NULL;
		ms = NULL;
		_m = _n = 0;
		knowM = knowN = false;
		lineNumber = -1;
		atEnd = false;
		lastError = GOOD;
	}

    public:

    /** Get a unique string describing this format. */
    	virtual const char* getName() const = 0;

    /** Determine if this format is sparse or dense.
     * @return true if it is a sparse format, false if it is a dense one
     */
	virtual bool isSparse() const = 0;

    /** Initialize this MatrixStreamReader.  Calls the initImpl method of the
      * subclass. */
	MatrixStreamError init( std::istream*, MatrixStream<Field>* );

    /** Get the next triple of row index, column index, value and store it into
     * the three referenced variables.  Uses the nextTripleImpl method of the
     * subclass. */
	MatrixStreamError nextTriple( size_t&, size_t&, Element& );

    /** Reads the next triple from the subclass nextTripleImpl method and saves
     * it to the savedTriples std::queue rather than returning it.  The error
     * returned is that given from the subclass method. */
	MatrixStreamError saveNext();
	
    /** Get the number of rows in this matrix, store it in the given int. */
	MatrixStreamError getRows(size_t&);

    /** Get the number of columns in this matrix, store it in the given int. */
	MatrixStreamError getColumns(size_t&);

    /** Get the line number that this reader is currently on. */
	virtual int getLineNumber() const { return lineNumber; }

    /** Virtual destructor. */
	virtual ~MatrixStreamReader() {
		while( !savedTriples.empty() ) savedTriples.pop();
	}
};

template <class Field>
class MatrixStream {
    public:
    	typedef typename Field::Element Element;

    private:
    	typedef std::pair<MatrixStreamReader<Field>*,std::stringstream*> RPair;
	typedef std::vector< RPair > RVector;

    /** The maximum number of characters read for each chunk by the addChars
      * method */
	static const int BUFLIMIT = 160;
	
    /** The maximum number of triples to read and save in an attempt to resolve
     * ambiguity by the automaticResolve() method. */
	static const int RESOLVE_LIMIT = 50;

    /** A buffer used to read and write data. */
	char buffer[BUFLIMIT];

    /** The underlying input stream from which data is being read. */
    	std::istream& in;

    /** A std::vector holding all the MatrixStreamReaders that have not yet failed*/
	RVector readers;

    /** The current state of matrix streaming. */
	MatrixStreamError currentError;

    /** The line number on which the last error occurred */
	int errorLineNumber;

    /** The Field to use in reading elements. */
	const Field& f;

    /** True if there is exactly one MatrixStreamReader in readers and it is
     * directly connected to in, rather than using a buffered std::stringstream.
     * False otherwise. */
	bool directStream;

    /** The chunk delimiter used by addChars */
	char delimiter;

    /** To ensure no one makes a copy of an instance of this class */
    	MatrixStream( const MatrixStream<Field>& ) {}

    /** Called by the constructors to get things going. */
	void init();

    /** Adds the given MatrixStreamReader to the readers std::vector. */
	void addReader( MatrixStreamReader<Field>* );

    /** Read ahead in the matrix up to RESOLVE_LIMIT triples in an attempt to
     * resolve ambiguity in formats.
     */
	bool automaticResolve();
	
    public:
    
    /** Constructor from an input stream.
     * @param fld The Field used to read Elements from the matrix.
     * @param in The input stream from which to read
     * @param delim The chunk delimited used by addChars.  This should be the
     *              last character in the matrix data.  If only one matrix is
     *              stored in a given file, this can be anything.  Otherwise, if
     *              delim is not the last character in the matrix, there could
     *              be extra data read from in and not put back.  Default is
     *              newline.
     */
    	MatrixStream( const Field& fld, std::istream& i, char delim = '\n' );
	
    /** Destructor */
	~MatrixStream();
	
    /** Read the next triple of row index, column index, value and store it in
     * the three referenced elements.
     * @return true iff the read succeeded.
     */
	bool nextTriple( size_t&, size_t&, Element& );

    /** Get the number of rows in the matrix and store it in the given size_t.
     * @return true iff the operation succeeded.
     */
	bool getRows(size_t&);

    /** Get the number of columns in the matrix and store it in the given size_t.
     * @return true iff the operation succeeded.
     */
	bool getColumns(size_t&);
	
    /** Get the number of rows and columns in the matrix and store them in the
     * given ints.
     * @return true iff the operation succeeded.
     */
	bool getDimensions( size_t&, size_t& );

    /** Get the current state of the stream.  Especially useful if called after
     * one of the four above operations on failure to get some information on
     * what caused the failure.
     */
	MatrixStreamError getError() const { return currentError; }

    /** Report the error to the error stream and return it.  Designed for throw
     * operations.
     */
     	MatrixStreamError reportError(const char*, int) const;

    /** If there is a reader in the GOOD state, get the line number it is on.
     * Otherwise, get the line number on which the last error occurred.
     */
	int getLineNumber() const;

    /** Get the Field that was passed to the constructor. */
	const Field& getField() const { return f; }

    /** Get characters from in and give them to the MatrixStreamReaders.
     * Called by the moreData() function of MatrixStreamReader.
     * @param eofReached A pointer to a pointer of the std::istream that reached an
     *                   EOF and initiated this function call.
     * @return true iff more characters were actually added to the stream
     */
	bool addChars( std::istream** eofReached = NULL );

    /** Get a brief description of the format of the matrix being read. */
	const char* getFormat() const;

    /** Tell if the matrix being read is sparse.
     * @return true iff the matrix is sparse.
     * NOTE: a return of false does NOT mean the matrix is dense.  Use isDense()
     */
	bool isSparse() const;

    /** Tell if the matrix being read is dense.
     * @return true iff the matrix is dense.
     * NOTE: a return of false does NOT mean the matrix is sparse.  Use
     *       isSparse()
     */
	bool isDense() const;
}; // end of class MatrixStream

} // end of namespace LinBox

#include "matrix-stream.inl"

#endif // MATRIX_STREAM_H
