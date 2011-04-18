/* Copyright (C) 2005 LinBox
 * Written by Dan Roche
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_matrix_stream_H
#define __LINBOX_matrix_stream_H

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
 */

#include <iostream>
#include <queue>
#include <vector>

namespace LinBox 
{ 
// namespace in which all LinBox code resides
//  Error codes enumeration
enum MatrixStreamError {
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
     */
    	std::istream* sin;

    /** A pointer to the MatrixStream that is using this reader.  Useful to
      * get an instance of the field via ms->getField().
      */
	MatrixStream<Field>* ms;
	
    /** The number of rows in the matrix.  This will be set by default to 0.
      */
	size_t _m;

    /** Indicates whether the row dimension is accurate */
    	bool knowM;
	
    /** Number of columns in the matrix.  Similar requirements as _m above. */
	size_t _n;

    /** Indicates whether the column dimension is accurate */
    	bool knowN;

    /** Indicates that the end of the matrix has been reached; no more calls to
      * nextTripleImpl will be made once this value is true.  This will
      * automatically be set to true if nextTripleImpl returns END_OF_MATRIX.
      */
	bool atEnd;

    /** Save the triple (m,n,v) onto the savedTriples std::queue. */
	void saveTriple(size_t m, size_t n, const Element& v);
	
    /** Read the next triple of row index, column index, value and store it in
     * the given references.
     * @return A MatrixStreamError indicating the success or failure of the
     *         operation
     */
	virtual MatrixStreamError nextTripleImpl(size_t&,size_t&,Element&) = 0;
	
    /** Read the first line of the matrix from the stream and attempt to
     * determine if it is of this reader's type.
     * @return A MatrixStreamError indicating the success or failure of the
     *         operation
     */
	virtual MatrixStreamError initImpl(const char* firstLine) = 0;

    /** A protected constructor that is called automatically when subclasses
     * are instantiated. */
	MatrixStreamReader() {
		sin = NULL;
		ms = NULL;
		_m = _n = 0;
		knowM = knowN = false;
		atEnd = false;
		lastError = GOOD;
	}

    public:

    /** Get a unique string describing this format. */
    	virtual const char* getName() const = 0;

    /** Get a (possibly) shortened version of the format name. */
    	virtual const char* shortName() const
	{ return getName(); }

    /** Determine if this format is sparse or dense.
     * @return true if it is a sparse format, false if it is a dense one
     */
	virtual bool isSparse() const = 0;

    /** Initialize this MatrixStreamReader.  Calls the initImpl method of the
      * subclass. */
	MatrixStreamError init(const char*,std::istream*,MatrixStream<Field>*);

    /** Get the next triple of row index, column index, value and store it into
     * the three referenced variables.  Uses the nextTripleImpl method of the
     * subclass. */
	MatrixStreamError nextTriple( size_t&, size_t&, Element& );

    /** Get the whole matrix as a dense (row-major) array of elements.
     * By default, this implementation just repeatedly calls nextTriple to
     * read in the whole matrix. Subclasses of dense formats should override
     * this behavior.
     * @param array The array to fill with entries. May be resized as needed.
     */
    	virtual MatrixStreamError getArray( std::vector<Element> &array );

    /** Reads the next triple from the subclass nextTripleImpl method and saves
     * it to the savedTriples std::queue rather than returning it.  The error
     * returned is that given from the subclass method. */
	MatrixStreamError saveNext();
	
    /** Get the number of rows in this matrix, store it in the given int. */
	MatrixStreamError getRows(size_t&);

    /** Get the number of columns in this matrix, store it in the given int. */
	MatrixStreamError getColumns(size_t&);

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
    /** The underlying reader for this matrix stream. */
	MatrixStreamReader<Field>* reader;

    /** The first line of text used to init. */
        char* firstLine;

    /** The maximum number of characters to read for the first line (to send to
     * the init functions).
     */
	static const int FIRST_LINE_LIMIT = 160;
	
    /** The underlying input stream from which data is being read. */
    	std::istream& in;

    /** The lineNumber is recorded in case the user wants to know at which line
      * an error occurs.  This will be updated automatically by any of the read
      * methods below if they encounter breaks; it is up to the subclasses to 
      * increment lineNumber if they read any newline characters independently.
      */
	int lineNumber;

    /** The current state of matrix streaming. */
	MatrixStreamError currentError;

    /** The line number on which the last error occurred */
	int errorLineNumber;

    /** True once any read functions have been called. */
    	bool readAnythingYet;

    /** The Field to use in reading elements. */
	const Field& f;

    /** To ensure no one makes a copy of an instance of this class */
	MatrixStream( const MatrixStream<Field>& )  ;// BB si {} pbm d'initialisation

    /** Called by the constructors to get things going. */
	void init();

    /** Adds the given MatrixStreamReader to the readers std::vector. */
	void addReader( MatrixStreamReader<Field>* );

    public:

    /** Constructor from an input stream.
     * @param fld The Field used to read Elements from the matrix.
     * @param in The input stream from which to read
     * @throws MatrixStreamError if an error occurs in reading the
     *         first line (i.e. on initialization).
     */
    	MatrixStream( const Field& fld, std::istream& i );
	
    /** Destructor */
	~MatrixStream() { delete reader; }
	
    /** Re initiliaze after one matrix has been read. */
	void newmatrix();

    /** Read some white space (if there is any). Using this method is preferable
     * to letting the input stream handle whitespace skipping because this 
     * method will update the line number when breaks are encountered.
     * @return true iff there is more data after the whitespace.
     */
   	bool readWhiteSpace();

    /** Read the next triple of row index, column index, value and store it in
     * the three referenced elements.
     * @return true iff the read succeeded.
     */
	bool nextTriple( size_t&, size_t&, Element& );

    /** Get the whole matrix as a dense (row-major) array of elements.
     * @param array The array to fill with entries. May be resized as needed.
     */
    	bool getArray( std::vector<Element> &array );

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
     * nextTriple or one of the get operations on failure to get some
     * information on what caused the failure.
     */
	MatrixStreamError getError() const { return currentError; }

    /** Report the error to the error stream and return it.  Designed for throw
     * operations.
     */
     	MatrixStreamError reportError(const char*, int) const;

    /** If the reader is in the GOOD state, get the line number it is on.
     * Otherwise, get the line number on which the last error occurred.
     */
	int getLineNumber() const;

    /** Get the Field that was passed to the constructor. */
	const Field& getField() const { return f; }

    /** Get a brief description of the format of the matrix being read. */
	const char* getFormat() const { return reader->getName(); }

    /** Get a very brief description of the matrix format. */
    	const char* getShortFormat() const { return reader->shortName(); }

    /** Tell if the matrix being read is sparse.
     * @return true if the matrix is sparse, false if it is dense.
     */
	bool isSparse() const { return reader->isSparse(); }
}; // end of class MatrixStream

}  // end of namespace LinBox

#include "matrix-stream.inl"

#endif // __LINBOX_matrix_stream_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
