#ifndef __FIELD_READER_ANALYZER
#define __FIELD_READER_ANALYZER

#include <string>

using std::string;

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/integer.h"


namespace LinBox {

	// empty class definition for error throw below
	class ReaderNotAField {};
	class NoFiniteExtendedInLinBox {};
	class NoTableFiniteInLinBox {};

	class FieldReaderAnalyzer {

	public:
		FieldReaderAnalyzer();
		FieldReaderAnalyzer(Reader &);
		~FieldReaderAnalyzer() {}
		
		// This method resets the analyzer with a new reader
		bool reset(Reader &);

		// These are the "generic field types" that are supported in
		// LinBox.  Each one of them has a default (sort of, though
		// we can't really do FiniteExtended without NTL, Givaro or
		// Lidia, w/o NTL the best we can do for Real is
		// Unparametric<double>, though w/ integer we could do
		// Unparametric<integer>.  Note that of course NotField
		// doesn't represent a field but the fact that the
		// user screwed up
		//
		static const int Finite = 0;
		static const int FiniteExtended = 1;
		static const int Rational = 2;
		static const int Real = 3;
		static const int Integer = 4;
		static const int Unknown = 5;
		
		static const int NotField = -1;
		
		// This method returns true if the Reader last used to
		// reset the analyzer was actually a field
		bool isField() const;

		// This method returns one of the values above, depending
		// on which field was last used to set the analyzer
		int whatType() const;

		// This method returns the string of the implDetail attribute
		// of the field, if the object analyzed was a field.  If it
		// wasn't a field, return the empty string (ie, the string
		// created by "string s;"
		//
		string implDetail() const;
		
		// returns the field's cardinality if the field is actually
		// a field, returns -1 otherwise
		integer cardinality() const;


		// Okay, this is the complicated method.
		// This method tries to create an instance of the field 
		// given in the Reader on the heap.  As arguments it takes
		// a functor templatized on the field type and a data 
		// structure type, and an instance of this data structure
		// type to be passed as the second argument of the functor.
		// This DSType is meant to allow the user to pass some kind
		// of data structure in to their functor call
		//
		// Note, if the Reader given to be analyzed does not describe
		// a field (or the reader attempts to initalize the field
		// and fails for whatever reason), an exception is thrown
		// Another special note for the author of functors to
		// work in this thing:  it is expected that the operator()
		// of this functor (or the function itself, if it's a function)
		// be templatized on the field type and take a Field* as it's
		// first argument, and expect to take a reference to a 
		// DSType as it's second argument, and to return a void*
		// It is also assumed that the responsibility for cleaning
		// up the field created is passed on to the functor (which
		// can in turn pass this responsibility on to whoever it likes)
		//
		// Sorry, I didn't mean to make this so complicated, but
		// I can think of 3 very different uses for this thing, 
		// and short of implementing this in Scheme, this is the best
		// I can think of to make the FieldReaderAnalyzer work in
		// all of them
		//
		template<class FunctorType, class DSType>
		void* makeField(const FunctorType &, DSType &) const;

		// This one is the same as the one above, but is designed
		// for functors that don't need another data argument to
		// be passed in.
		//
		template<class FunctorType>
		void* makeField(const FunctorType &) const;	


	private:
		mutable Reader _R;
		bool _isField;
		int _fieldType;
		string _implDetail;
		integer _card;
	};


			
	
