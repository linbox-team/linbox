/* Copyright (C) 2010 LinBox
 * Written by <?>
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

#ifndef __LINBOX_reader_H
#define __LINBOX_reader_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <stack>

// for strchar
#include <cstring>
// for isdigit
#include <cctype>

#include "linbox/integer.h"
#include "linbox/util/xml/linbox-writer.h"

#include "xml-tree.h"

/*
  LinBox basic_reader Class -
  A helpful tool that lets a LinBox user who wants to convert an XML structure
  into a LinBox object to be able to do so.  The Reader object provides the 
  following tools:

  - Parse an std::istream containing an XML structure into a Node structure
    bool parse(std::istream &, const char* encoding = "US-ASCII")

  - Check whether a given std::string is an integer
    bool isInt(const std::string&)

  - Check whether a given std::string is a std::list of space seperated ints
    bool isIntList(const std::string&)

  - Convert a given std::string to an int if possible
    bool toInt(int&, const std::string&)

  - Convert a given std::string to a GMP int if possible
    bool toGMP(integer&, const std::string &);

  - Convert a given std::string to a std::vector of ints if possible
    bool toIntList(std::vector<int>&, const std::string &)

  - Convert a given std::string to a std::vector of LinBox integers if possible
    bool toGMPList(std::vector<integer>&, const std::string &);

  - Read a Node's Tag
    bool tagName(std::string &);

  - given an attribute name, get the value of that attribute (return false
    if it isn't there)
    bool getAttribValue(const std::string &name, std::string &value);

  - Move through a Node's children

      - (we always start at the first child)
      - Move to the next child, if possible
        bool getNextChild(); // returns false if at last child

      - Move to the previous child, if possible
        bool getPrevChild(); // False if at first child

      - Check whether at the beginning or the end
        bool firstChild();
	bool lastChild();

      - Check whether there are children at all?
        bool haveChildren();

      - Check whether the current child is a DataNode, a TagNode
        bool isChildDataNode();
	bool isChildTagNode();

      - If a DataNode, return the std::string comprising that Data
        bool getChildText(std::string &dataHolder);

      - If a TagNode, go into that child (the previous place is saved)
        bool traverseChild();

      - Go back up to a Parent Node if possible
        bool upToParent();

    The class will inherit from the XMLTree class (as the reader is nothing
more than an XMLTree with some API that makes it easy for us to read the
data and convert to LinBox objects), so we inherit to get the parse functions
and some internal data strucures to work with.  We also have a currentNode*
for the current node we are on, a child std::list iterator for the current child
we are on, and a std::stack of a struct containing a TagNode* and a child std::list
iterator which we use to save our place

*/

namespace LinBox {

	struct TagTreeFrame {
		node* ptr;
		std::list<node*>::iterator it;
	};


	class BasicReader : private XMLTree {
		//		friend class Writer;


	  public:
		// Constructors
		BasicReader();
		BasicReader(std::istream &, const char* encoding = "US-ASCII");
		BasicReader(const BasicReader &);
		BasicReader(const Writer &);
		~BasicReader();

		// Note - By design this operator _DOES NOT_ create a true
		// copy, it is a hack designed to have multiple Reader objects
		// all tied together at the hip, the idea being that the user
		// creates one Reader for the whole of the XML, and any
		// cloned Readers that operate on a portion all modify
		// their error state
		// (oh, and only the original is allowed to call parse)
		//
		const BasicReader &operator=(const BasicReader &);

		// the parser
		bool parse(std::istream &, const char* encoding = "US-ASCII");

		// methods for finding out what happened when you tried reading the XML
		bool initalized() const;
		bool haveError() const;


		const std::string &getErrorString() const;
		int getErrorCode() const;
		int getParseErrorLine() const;
		void setErrorString(const std::string &);
		void setErrorCode(int);
		void setParseErrorLine(int);

		static const int NO_ERROR = 0;
		static const int NO_INIT = 1;
		static const int FAILED_INIT = 2;
		static const int IMPROPER_COPY_OVERWRITE = 12;

		// static helper functions
		static bool isNum(const std::string &);
		static bool isUnsignedNum(const std::string &);

		// template version which the four below overload
		template<class Num>
		static bool toNum(Num &, const std::string &);

		template<class Field>
		bool toNum(const Field &, typename Field::Element &e);

		//static bool toInt(int&, const std::string &);
		//static bool toLong(long&, const std::string &);
		//static bool toSizet(size_t&, const std::string &);
		//static bool toGMP(integer &, const std::string &);

		static bool isNumVector(const std::string &);
		static bool isUnsignedNumVector(const std::string &);

		// template version of the four below
		template<class Num>
		static bool toNumVector(std::vector<Num>&, const std::string &, bool = false);

		template<class Field>
		bool toNumVector(const Field &, std::vector<typename Field::Element> &);

		//static bool toIntVector(std::vector<int>&, const std::string &);
		//static bool toLongVector(std::vector<long>&, const std::string &);
		//static bool toSizetVector(std::vector<size_t>&, const std::string &);
		//static bool toGMPVector(std::vector<integer>&, const std::string &);

		// get the tag name & attributes
		bool isTag() const;
		bool isText() const;

		bool getTagName(std::string &name) const;
		bool getAttribValue(const std::string &name, std::string &value) const;

		// API for moving through the tree & accessing parts of the
		// subtree

		bool getNextChild();
		bool getPrevChild();

		bool isFirstChild() const;
		bool isLastChild() const;
		bool haveChildren() const;

		size_t numChildren() const;
		
		bool isChildText() const;
		bool isChildTag() const;

		bool getChildText(std::string &dataHolder) const;
		bool traverseChild();
		bool upToParent();


		// These members correspond to size_t N calls to
		// upToParent(), traverseChild(), getNextChild() and
		// getPrevChild().  Note, the given # of calls will be made
		// regardless of whether or not errors begin to be reported
		BasicReader &Up(size_t);
		BasicReader &Down(size_t);
		BasicReader &Right(size_t);
		BasicReader &Left(size_t);
		
		bool getText(std::string &text) const;



	protected:
		node* currentNode;
		TagNode* TNode;
		DataNode* DNode;
		std::list<node*>::iterator currentChild;
		std::stack<TagTreeFrame> parentStore;


	  private:
		std::string* errorString;
		int* ErrorCode, *parseErrorLine, *ref_count;
		bool original;


		// takes a char and returns the integer value of the char
		// (subtracts off the ascii value of '0')
		static int char2digit(char);



	};


	/* LinBox Reader - is a BasicReader w/ an expanded API for easier XML validation.
	 * Uses the concept of an "expect" method - you expect a certain aspect of the current peice of the 
	 * parsed document to have a certain form, and if so, translate that form into a useable value.  Automates
	 * the process of validating a strucutre, I just tell it what I expect and where.  If what I'm expecting isn't there,
	 * the expect method writes an error std::string describing the problem and returns false (and I assume the read function
	 * or method using this Reader returns false as well.
	 * Provides the following API:
	 * - expectTagName : check that the Tag Name of the currentNode is equal to the provided std::string
	 * - expectAttribute"type" : take a name and a reference to a variable of "type" and check that the 1)
	 *     there is an attribute of the name given and 2) it's value is of the type given.  If both of these hold
	 *     true, convert the value to the type given and initalize the reference to that type.  The valid types in this case
	 *     are: std::string, Int, GMP, Long and Sizet
	 * - expectNChildren : Expect a certain # of children
	 * - expectChildData"type": Take in a reference to a variable of "type" and checks that 1) The currentChild is a DataNode,
	 *    2) the data of the currentChild is of the type specified.  If so, convert the Data of the child to the given type
	 *    and initalize the reference to that value.  Valid types are: std::string, Int, GMP, Long, Sizet, IntVector, GMPVector,
	 *    LongVector, SizetVector
	 * - expectChildTag: Expect that the currentChild is a TagNode
	 * 
	 * Using public inheritance, it is ensured that this Reader has the basic contructors, parse(), traverseChild(), upToParent()
	 * etc are still there for use.
	 */

	// forward declaration for ReaderIterator
	class ReaderIterator;


	class Reader : public BasicReader {


	public:
		typedef ReaderIterator iterator;
		typedef const ReaderIterator const_iterator;

		Reader();
		Reader(std::istream &, const char* encoding = "US-ASCII");
		Reader(const Reader&);
		~Reader();

		// Note - This operator= _DOES NOT_ create a true copy,
		// but instead creates a new "hydra-head" on the underlying
		// XML Tag object.  This operator is intended to be used 
		// mainly to allow multiple instances of a Reader to work on
		// differing peices of an XML document and share Error state.
		// It's used A) In conjunction with the ReaderIterators
		// to spawn little Readers to be used when the user is calling
		// *iter on iter of type ReaderIterator, B) In conjunction with
		// FieldReaderAnalyzer and BlackboxReaderFactory classes, which
		// both have internal Readers.  You hand each of these a Reader
		// and it mainly makes a copy of the tag structure in memory
		//
		const Reader &operator=(const Reader &);


		// New Error codes for this type of Reader
		static const int TAGNAME_ERROR = 3; 
		static const int MISSING_ATT = 4;
		static const int WRONG_ATT_TYPE = 5;
		static const int WRONG_N_CHILD = 6;
		static const int CHILD_NOT_DATA = 7;
		static const int WRONG_DATA_TYPE = 8;
		static const int CHILD_NOT_TAG = 9;
		static const int NOT_TEXT = 10;
		static const int NOT_TAG = 11;
		static const int OVFLOW = 13;
		static const int OTHER = 14;

		// Expect a Tag, that's it
		bool expectTag();
		bool checkTag() const;
		bool expectText();
		bool checkText() const;

		
		// Expect a Tag with a certain name
		bool expectTagName(const std::string &);
		bool expectTagName(const char*);


		// Check for a Tag with a certain name
		bool checkTagName(const std::string &) const;
		bool checkTagName(const char*) const;

		template<class Num>
		bool expectTagNum(Num &);

		template<class Num>
		bool checkTagNum(Num &);

		template<class Num>
		bool expectTagNumVector(std::vector<Num> &);

		template<class Num>
		bool checkTagNumVector(std::vector<Num> &);

		// template version that uses a field
		template<class Field>
		bool expectFieldNum(const Field &, typename Field::Element &);

		template<class Field>
		bool checkFieldNum(const Field &, typename Field::Element &);

		template<class Field>
		bool expectFieldNumVector(const Field &, std::vector<typename Field::Element>&);

		template<class Field>
		bool checkFieldNumVector(const Field &, std::vector<typename Field::Element>&);


		// Check for a Text Node
		bool expectTextString(std::string &);

		// template version of all of these
		template<class Num>
		bool expectTextNum(Num &);

		//bool expectTextInt(int &);
		//bool expectTextGMP(integer &);
		//bool expectTextLong(long &);
		//bool expectTextSizet(size_t &);


		template<class Num>
		bool expectTextNumVector(std::vector<Num>&, bool = false);


		//bool expectTextIntVector(std::vector<int>&);
		//bool expectTextGMPVector(std::vector<integer>&);
		//bool expectTextLongVector(std::vector<long>&);
		//bool expectTextSizetVector(std::vector<size_t>&);
		

		bool checkTextString(std::string &) const; 

		template<class Num>
		bool checkTextNum(Num&) const;

		//bool checkTextInt(int &) const;
		//bool checkTextGMP(integer &) const;
		//bool checkTextLong(long &) const;
		//bool checkTextSizet(size_t &) const;

		template<class Num>
		bool checkTextNumVector(std::vector<Num>&, bool = false) const;

		//bool checkTextIntVector(std::vector<int>&) const;
		//bool checkTextGMPVector(std::vector<integer>&) const;
		//bool checkTextLongVector(std::vector<long>&) const;
		//bool checkTextSizetVector(std::vector<size_t>&) const;
		

		// Expect an attribute with a certain value type, and get the value in that type
		bool expectAttributeString(const std::string &name, std::string &value);

		template<class Num>
		bool expectAttributeNum(const std::string &name, Num &value);

		//bool expectAttributeInt(const std::string &name, int &value);
		//bool expectAttributeGMP(const std::string &name, integer &value);
		//bool expectAttributeLong(const std::string &name, long &value);
		//bool expectAttributeSizet(const std::string &name, size_t &value);


		bool checkAttributeString(const std::string &name, std::string &value) const;
		template<class Num>
		bool checkAttributeNum(const std::string &name, Num &value) const;

		//bool checkAttributeInt(const std::string &name, int &value) const;
		//bool checkAttributeGMP(const std::string &name, integer &value) const;
		//bool checkAttributeLong(const std::string &name, long &value) const;
		//bool checkAttributeSizet(const std::string &name, size_t &value) const;

		/* same as above, just use a char* instead for the name
		   At the moment I'm not sure why I'm removing these, but I have this feeling like it's
		   a good idea, so I'm going to go with it
	        bool expectAttributeString(const char*name, std::string &value);
	        bool expectAttributeInt(const char* name, int &value);
	        bool expectAttributeGMP(const char* name, integer &value);
		bool expectAttributeLong(const char* name, long &value);
		bool expectAttributeSizet(const char* name, size_t &value); 
		*/

		// Expect a certain number of children
		bool expectNChildren(const size_t &N);

		// Expect at least a minimum number of children
		bool expectAtLeastNChildren(const size_t &N);


		// Expect the current child to be a DataChild, and for that data to be of a certain type
		// initalize for that type
		bool expectChildTextString(std::string &);
		
		template<class Num>
		bool expectChildTextNum(Num &);

		//bool expectChildTextInt(int &);
		//bool expectChildTextGMP(integer &);
		//bool expectChildTextLong(long &);
		//bool expectChildTextSizet(size_t&);

		template<class Num>
		bool expectChildTextNumVector(std::vector<Num>&, bool = false);

		//bool expectChildTextIntVector(std::vector<int>&);
		//bool expectChildTextGMPVector(std::vector<integer>&);
		//bool expectChildTextLongVector(std::vector<long>&);
		//bool expectChildTextSizetVector(std::vector<size_t>&);

		// Check for each of these functions, set no error if they are false
		bool checkChildTextString(std::string &) const;

		template<class Num>
		bool checkChildTextNum(Num &) const;

		//bool checkChildTextInt(int &) const;
		//bool checkChildTextGMP(integer &) const;
		//bool checkChildTextLong(long &) const;
		//bool checkChildTextSizet(size_t&) const;

		template<class Num>
		bool checkChildTextNumVector(std::vector<Num> &, bool = false) const;
		
		//bool checkChildTextIntVector(std::vector<int> &) const;
		//bool checkChildTextGMPVector(std::vector<integer> &) const;
		//bool checkChildTextLongVector(std::vector<long> &) const;
		//bool checkChildTextSizetVector(std::vector<size_t> &) const;
	       
		
		// Expect the next child to be a Tag Node
		bool expectChildTag();
		
		// Check for it
		bool checkChildTag() const;

		// ReaderIterator access methods.  Creates a ReaderIterator for the first child
		// If this is a text node & not a tag node, this ReaderIterator is empty, and as such will
		// have no state
		//
		iterator begin();
		const_iterator begin() const;
		iterator end();
		const_iterator end() const;
		
		Reader & Up(size_t N = 1);
		Reader &Down(size_t N = 1);
		Reader &Right(size_t N = 1);
		Reader &Left(size_t N = 1);


	protected:



		
	};




	// Default constructor.  Just call the XMLTree default constructor
	//
	BasicReader::BasicReader() : XMLTree() {
		errorString = new std::string;
		parseErrorLine = new int;
		ErrorCode = new int;
		ref_count = new int(0);

		setErrorString("Reader not initalized.");
		*parseErrorLine = 0;
		setErrorCode(NO_INIT);
		original = true;
		currentNode = NULL;
		TNode = NULL;
		DNode = NULL;
		currentChild = NULL;

	}

	BasicReader::BasicReader(const Writer &W) : XMLTree() {
		errorString = new std::string;
		parseErrorLine = new int;
		ErrorCode = new int;
		ref_count = new int(0);

		setErrorString("No Error");
		*parseErrorLine = 0;
		setErrorCode(NO_ERROR);
		original = true;
		TNode = W.currentNode->clone();
		currentNode = TNode;
		DNode = NULL;
		currentChild = TNode->children.begin();
	}


	// std::istream constructor.  Call the XMLTree std::istream constructor and
	// if the parsing construction was successful, initalize the
	// currenNode ptr and currentChild
	//
	BasicReader::BasicReader(std::istream &dataStream, const char* encoding) : XMLTree() {
		errorString = new std::string();
		parseErrorLine = new int;
		ErrorCode = new int;
		ref_count = new int(0);
		original = true;

		try {
			XMLTree::parse(dataStream, encoding);
			setErrorString("No Error.");
			setErrorCode(NO_ERROR);
			*parseErrorLine = 0;
			currentNode = &Tree;
			TNode = dynamic_cast<TagNode*>(currentNode);
			DNode = NULL;
			currentChild = TNode->children.begin();
		}
		catch(Error e) {
			setErrorString(e.error_descrip);
			setErrorCode(FAILED_INIT);
			*parseErrorLine = e.param;
			
		}
			
	}

	// copy constructor - Makes two readers that point to the same tag structure, and also have the
	// same error conditions & stuff.  Useful primarily for the construction of ReaderIterator s
	// NOTE - it's okay that we're setting the initFlag to the
	// value of the copied BasicReader
	// we only allow after the fact parsing on original objects, not
	// on copied Readers.  Trying to parse a copied Reader will
	// set error flags in the BasicReader
	//
	BasicReader::BasicReader(const BasicReader &Bin) : XMLTree() {
		errorString = Bin.errorString;
		ErrorCode = Bin.ErrorCode;
		parseErrorLine = Bin.parseErrorLine;
		ref_count = Bin.ref_count;
		*ref_count = *ref_count + 1;

		currentNode = Bin.currentNode;
		TNode = Bin.TNode;
		DNode = Bin.DNode;
		currentChild = Bin.currentChild;
		parentStore = Bin.parentStore;
		initFlag = true;
		original = false;
	}
	
	const BasicReader &BasicReader::operator=(const BasicReader &Bin) {

		// prevent a disasterous bug
		if(this != &Bin) {

			// first delete all the old dynamically allocated memory
			delete errorString;
			delete ErrorCode;
			delete parseErrorLine;
			delete ref_count;

			
			errorString = Bin.errorString;
			ErrorCode = Bin.ErrorCode;
			parseErrorLine = Bin.parseErrorLine;
			ref_count = Bin.ref_count;
			*ref_count = *ref_count + 1;

			currentNode = Bin.currentNode;
			TNode = Bin.TNode;
			DNode = Bin.DNode;
			currentChild = Bin.currentChild;
			parentStore = Bin.parentStore;
			initFlag = true;
			original = false;
		}

		return *this;
	}
	

	// In order to have a shared Error reporting between Reader iterators
	// and the the spawning Reader, I have to the make the Error
	// reporting shared memory.  For this reason, the destructor is
	// no longer empty.  What we do here is check the shared ref_count
	// If the ref_count is 0, only 1 object is pointing to it, so it's
	// okay to free the memory.  Otherwise, just decrement the reference
	// count
	//
	BasicReader::~BasicReader() {
		if(*ref_count == 0) {
			delete errorString;
			delete parseErrorLine;
			delete ErrorCode;
			delete ref_count;
		}
		else *ref_count = *ref_count - 1;
	}

	// parse - Takes in an std::istream which (we hope) has XML text
	// which we will turn into a LinBox object.  Performs
	// the base calss parse, then initalize the derived class
	// members if the parse was successful
	//
	bool BasicReader::parse(std::istream &dataStream, const char* encoding) {
		if(!original) {
			setErrorString("Tried to re-parse a copy of the Reader, not okay.");
			setErrorCode(IMPROPER_COPY_OVERWRITE);
			return false;
		}



		try { 
			XMLTree::parse(dataStream, encoding); 
			setErrorString("No Error."); 
			setErrorCode(NO_ERROR); 
			*parseErrorLine = 0; 
			currentNode = &Tree; 
			TNode = dynamic_cast<TagNode*>(currentNode); 
			DNode = NULL; 
			currentChild = TNode->children.begin();
			return true; 
		} 
		catch(Error e) {  
			setErrorString(e.error_descrip); 
			*parseErrorLine = e.param;
			setErrorCode(FAILED_INIT);
			return false; 
		}
	}

	// initalized - checks if the parser is properly initalized, and if not returns false
	bool BasicReader::initalized() const {
		return initFlag;
	}


	// haveError - Checks whether the BasicReader has an error code that isn't NO_ERROR (0).  If it does,
	// the BasicReader does haveError
	bool BasicReader::haveError() const {
		return *ErrorCode != NO_ERROR;
	}


	// getErrorString - returns the Error std::string associated with the current state of the
	// BasicReader.  If there is no error, it will return a std::string containing "No Error".
	// If the reader is not yet initalized, it will return a std::string containing "Reader not initalized
	// If there was a problem parsing the XML, the text description of that problem is returned along w/
	// the line number that the error occured on
	//
	const std::string &BasicReader::getErrorString() const {
		return *errorString;
	}

	// getParseErrorLine - returns the line of XML associated w/ the XML parsing error, if such an
	// error occured.  If no error has occured (or no attempt has yet been made to initalize the reader)
	// a 0 will be returned
	//
	int BasicReader::getErrorCode() const {
		return *ErrorCode;
	}

	int BasicReader::getParseErrorLine() const {
		return *parseErrorLine;
	}


	// setErrorString - sets the errorCode of the class.  Only available to derived classes
	void BasicReader::setErrorString(const std::string &inString) {
		*errorString = inString;
		return;
	}

	// setErrorCode - sets the ErrorCode of the object
	void BasicReader::setErrorCode(int code) {
		*ErrorCode = code;
		return;
	}


	// isInt - The first public utility function.  Takes in a 
	// std::string and checks whether the contens of that std::string are a number
	// a linear call, it works as follows:  First check if the first
	// character is either a digit or a -.  Then, check whether all
	// subsequent characters are numbers.  If not, return false.
	// if we get to the end, return true (also, a std::string that has
	// a - at the front MUST have another character behind it
	//
	bool BasicReader::isNum(const std::string &source) {
		
		size_t i;

		if( source.length() == 0) {
			return false;
		}
		else if(source[0] == '-') {
			if(source.length() == 1)
				return false;
			else {
				for(i = 1; i < source.length(); ++i) {
					if(!isdigit(source[i]))
						return false;
				}
				return true;
			}
		}
		else {
			for(i = 0; i < source.length(); ++i) {
				if(!isdigit(source[i]))
					return false;
			}
			return true;
		}
	}

	// isUnsignedNum - a simpler version of the isNum predicate that
	// checks solely for strings representing non-negative numbers.
	// used by toSizet
	bool BasicReader::isUnsignedNum(const std::string &source) {
		size_t i;

		for(i = 0; i < source.length(); ++i) 
			if(!isdigit(source[i])) return false;
		
		return true;
	}

			
	// toNum - Takes a std::string and converts that std::string to an int
	// essentially uses an istreamstream to attempt the conversion.
	// If the type in question isn't built-in or doesn't have an
	// overloaded operator>>, it's your own fault
	// Note that in this general case, there is no overflow protection,
	// it just returns true.  Good Luck! :-)
	//
	template<class Num>
	bool BasicReader::toNum(Num &dest, const std::string &source) {
	
		std::istringstream iss(source);
		iss >> dest;

		return true;
	}


	// the old function
	//		size_t i;
	//		Num oldValue, newValue;
	//		int temp;
	//
	//		if(source[0] == '-') {
	//			oldValue = newValue = 0;
	//			for(i = 1; newValue >= oldValue && i < source.length(); ++i ) {
	//				temp = BasicReader::char2digit(source[i]);
	//				oldValue = newValue;
	//				newValue = 10 * oldValue + temp;
	//			}
	//			if(newValue < oldValue) 
	//				return false;

	//			dest = -1 * newValue;
	//			return true;
	//		}
	//		else {
	//			oldValue = newValue = 0;
	//			for(i = 0; newValue >= oldValue && i < source.length(); ++i) {
	//				temp = BasicReader::char2digit(source[i]);
	//				oldValue = newValue;
	//				newValue = 10 * oldValue + temp;
	//			}
	//			if(newValue < oldValue)
	//				return false;
	//
	//			dest = newValue;
	//			return true;
	//		}
	//	}


	// toLong - Just like toLong, but convert to a long instead of
	// to an int
	//
	// Thanks to template version above, this function
	// becomes unecessary

	//	bool BasicReader::toLong(long &dest, const std::string &source) {
	
	//	size_t i;
	//	long oldValue, newValue, temp;

	//		if(source[0] == '-') {
	//			oldValue = newValue = 0;
	//			for(i = 1; newValue >= oldValue && i < source.length(); ++i ) {
	//		temp = BasicReader::char2digit(source[i]);
	//		oldValue = newValue;
	//		newValue = 10 * oldValue + temp;
	//	}
	//	if(newValue < oldValue) 
	//		return false;
	//
	//	dest = -1 * newValue;
	//	return true;
	//}
	//else {
	//	oldValue = newValue = 0;
	//	for(i = 0; newValue >= oldValue && i < source.length(); ++i) {
	//		temp = BasicReader::char2digit(source[i]);
	//		oldValue = newValue;
	//		newValue = 10 * oldValue + temp;
	//	}
	//	if(newValue < oldValue)
	//		return false;
	//
	//	dest = newValue;
	//return true;
	//		}
	//}


	// toNum - template specialization for size_t
	// Fails if an arithmetic overflow occurs
	//
	template<>
	bool BasicReader::toNum(size_t & dest, const std::string &source) {

		size_t oldValue, newValue, temp, i;

		oldValue = newValue = 0;
		for(i = 0; newValue >= oldValue && i < source.length(); ++i) {
			temp = BasicReader::char2digit(source[i]);
			oldValue = newValue;
			newValue = 10 * oldValue + temp;
		}

		if(newValue < oldValue) return false;
		
		dest = newValue;
		return true;
	}
		

	// toNum - template specialization that 
	// takes a std::string and attempts to convert that std::string to
	// a GMP int.  This one is alot easier, as there is no
	// question of the std::string being too long, as well as
	// having a majority of the infrastructure for the conversion
	// already available
	//
	template<>
	bool BasicReader::toNum(integer &dest, const std::string &source) {

		dest = Integer(source.c_str());
		return true;
	}

	// this function presumes that the current tag represents a field
	// element, and using the field, try to initalize this num
	template<class Field>
	bool BasicReader::toNum(const Field &F, typename Field::Element &e) {

		return F.fromTag(*this, e);
	}

	// This function presumes that every child of the current point in 
	// the tree is a field element, and initalizes a std::vector of these
	// elements
	//
	template<class Field>
	bool BasicReader::toNumVector(const Field &F, std::vector<typename Field::Element> &v) {

		typedef typename Field::Element Element;
		Element e;

		v.clear();

		if(haveChildren()) {

			traverseChild();
			if(F.fromTag(*this, e))
				v.push_back(e);
			else {
				upToParent();
				return false;
			}
			upToParent();

			while(!isLastChild()) {
				getNextChild();
				traverseChild();
				if( F.fromTag(*this, e)) 
					v.push_back(e);
				else {
					upToParent();
					Left(v.size());
					return false;
				} 
				upToParent();
			}
		}

		Left(v.size() - 1);
		return true;
	}
		




	// An number std::vector is a std::string w/ the following property: The std::string
	// contains one or more integers in std::string form seperated by a 
	// single whitespace character.  This function takes in a std::string
	// and determines whether the std::string has such a property.
	// The function works by making use of the std::string API (namely the find and substr)
	// functions to pick out substrings between single whitespace characters, then checking
	// whether the substring is a valid number using isNum
	//
	//
	bool BasicReader::isNumVector(const std::string &source) {
		
		size_t start = 0, offset;

		if(source.length() == 0) return false;

		while(start < source.length() ) {
			// offset first gets the location of the next space character
			offset = source.find(' ', start);
			// if we are on the last substring
			if(offset == std::string::npos) offset = source.length() - start;
			// otherwise
			else offset -= start;
			
			// check whether this substring is a Num, if not false
			if(!isNum(source.substr(start, offset))) return false;

			// set the start of the next substring to be beyond the beginning of the current substring
			start += offset + 1;
		}
		return true;
	}
		
	// isUnsignedNumVector - A simpler version of isNumVector that
	// checks only for std::string representations of vectors of unsigned
	// (non-negative) integer values.  Included mainly for
	// toSizetVector
      	//
	bool BasicReader::isUnsignedNumVector(const std::string &source) {

		size_t start = 0, offset;

		if(source.length() == 0) return false;
		while(start < source.length() ) {
			offset = source.find(' ', start);
			if(offset == std::string::npos) offset = source.length() - start;
			else offset -= start;

			if(!isUnsignedNum(source.substr(start, offset))) return false;
			start += offset + 1;
		}

		return true;
	}

		

	// toNumVector - Takes a std::string which represents a NumVector
	// and attempts to convert this to an actual Num std::vector
	// This function returns false if one of the numbers in this
	// number std::vector is too big to be held in a Num.  Notice that
	// this function WILL overwrite the std::vector given with all number
	// values before the first number that is too large, if the 
	// std::string is in fact a number std::vector
	//
	template<class Num>
	bool BasicReader::toNumVector(std::vector<Num> &vect , const std::string &source, bool addOne) {
		
		Num holder;
		size_t start = 0, offset;

		vect.clear();
		while(start < source.length()) {
			offset = source.find(' ', start);
			if(offset == std::string::npos) offset = source.length() - start;
			else offset -= start;

			if(!toNum(holder, source.substr(start, offset))) return false;
			else {
				if(addOne) 
					vect.push_back(holder + 1);
				else
					vect.push_back(holder);
				
				start += offset + 1;
			}
		}

		return true;
	}
		

	// due to the introduction of templates, this member has
	// become depricated
	// toLongVector - This function is about the same as toIntVector above except that
	// it takes a std::vector of long instead of a std::vector of int.  Everything else is the same
	//

	//	bool BasicReader::toLongVector(std::vector<long> &vect, const std::string &source) {
		
	//	long holder;
	//	size_t start = 0, offset;

	//		vect.clear();
	//              while(start < source.length()) {
	//	offset = source.find(' ', start);
	//		if(offset == std::string::npos) offset = source.length() - start;
	//		else offset -= start;

	//			if(!toLong(holder, source.substr(start, offset))) return false; // overflow
	//			else {
	//				vect.push_back(holder);
	//				start += offset + 1;
	//			}
	//		}

	//		return true;
	//	}

	
	// this method has become depricated due to the introduction of
	// templates, and has been removed
	// toNumVector - size_t specialization
	//Take in a std::string that represents an unsigned number std::vector and return a std::vector of size_t
	// comprising the numbers of that std::list.  This one is similar (but shorter) than the 2 functions above.
	// It again will overwrite the current std::vector so long as no number overflows
	//

	//	bool BasicReader::toSizetVector(std::vector<Num> &vect, const std::string & source) {

	//		size_t holder, start = 0, offset;

	//		vect.clear();
	//		while(start < source.length()) {
	//			offset = source.find(' ', start);
	//			if(offset == std::string::npos) offset = source.length() - start;
	//			else offset -= start;

	//			if(!toSizet(holder, source.substr(start, offset))) return false; // overflow
	//			else {
	//				vect.push_back(holder);
	//				start += offset + 1;
	//			}
//		}
//		return true;
//	}


	// so depricated, so removed
	// toGMPVector - Takes in a std::string that represents a number std::vector and return a std::vector of LinBox integer
	// objects.  This setup works by making substrings of the original std::string, each substring representing a
	// number in the number std::vector, then using the std::string constructor of the LinBox integer class to create
	// the integer.  This function will always return a fully initalized std::vector (provided the original std::string is in
	// the proper format) as there is no problem with overflow.
	//
	//	bool BasicReader::toGMPVector(std::vector<integer> &vect, const std::string &source) {

	//		integer holder;
	//		size_t start = 0, offset;

	//		vect.clear();
	//		while(start < source.length()) {
	//			offset = source.find(' ', start);
	//			if(offset == std::string::npos) offset = source.length() - start;
	//			else offset -= start;
			
	//			toGMP(holder, source.substr(start, offset));
	//			vect.push_back(holder);
	//			start += offset + 1;
	//		}
		
	//		return true;
	//	}

	// isTag - Checks whether the current node of the reader is a Tag or not
	// works because at any one time, only one of TNode & DNode won't be null
	bool BasicReader::isTag() const {
		return TNode != NULL;
	}

	// isText - Checks whether the current node of the readfer is Text or not
	bool BasicReader::isText() const {
		return DNode != NULL;
	}



	// getTagName - Returns the name of the current Tag
	bool BasicReader::getTagName(std::string &name) const {
		if(isTag()) {
			name = TNode->tag;
			return true;
		}
		else return false;
		
	}


	// getAttribValue - Takes in the name of an attribute and
	// returns the value of that attribute, if it is there
	// If the attribute is not there, return false
	//
	bool BasicReader::getAttribValue(const std::string &name, std::string &value) const {


		if(isTag()) {

			std::list<std::string>::iterator it;
			for(it = TNode->attrib.begin(); it != TNode->attrib.end() && *it != name; ++it) ++it;
			
			if(it == TNode->attrib.end()) return false;
			else {
				++it;
				value = *it;
				return true;
			}
		}
		else return false;
	}

	// getText - Takes in a std::string and if that std::string is text, return the std::string
	bool BasicReader::getText(std::string &writeable) const {
		if(isText()) {
			writeable = DNode->data;
			return true;
		}
		else return false;
	}

	
	// getNextChild - Moves to the next child in the Child std::list of the currentNode.  If we are already at the last child,
	// hold here and return false
	bool BasicReader::getNextChild() {
		if( !haveChildren() || isLastChild() ) return false;
		else ++currentChild;
		return true;
	}

	// getPrevChild - Moves to the previous child in the Child std::list of the currentNode.
	// if we are already at the first node, return false and go no further
	bool BasicReader::getPrevChild() {
		if( !haveChildren() || isFirstChild() ) return false;
		else --currentChild;
		return true;
	}

	// numChildren - returns the number of children of this Tag Node
	size_t BasicReader::numChildren() const {
		return (isTag() ? TNode->children.size() : 0);

	}
	

	// isFirstChild - Checks whether the current child is in fact the first child of the
	// current Tag
	bool BasicReader::isFirstChild() const {
		return (haveChildren() ? currentChild == TNode->children.begin() : false);
	}

	//isLastChild - Checks whether the current child is in fact the last child of the
	// current Tag
	bool BasicReader::isLastChild() const {
		if(haveChildren()) {

			std::list<node*>::iterator it = currentChild;
			++it;
			return it == TNode->children.end();
		}
		else return false;
	}

	// haveChildren - Returns true if the currentNode has children
	bool BasicReader::haveChildren() const {
		if(isTag()) 
			return !TNode->children.empty();
		else
			return false;
	}

	
	// isChildDataNode - Check whether the Child we are currently on is a DataNode 
	bool BasicReader::isChildText() const {
		if(haveChildren()) 
			return (*currentChild)->Ident == Data;
		else
			return false;
	}

	// isChildTagNode - Check whether the Child we are currently on is a TagNode
	bool BasicReader::isChildTag() const {
		if(haveChildren())
			return (*currentChild)->Ident == Tag;
		else
			return false;
	}

	// getChildText - If the currentChild is a DataNode, this function writes the
	// contents of the DataNode to the input std::string and returns true.  Otherwise, this
	// function returns false
	bool BasicReader::getChildText(std::string &dataHolder) const {
		node* nPtr;
		DataNode* dPtr;

		if(!isChildText() ) return false;

		nPtr = *currentChild;
		dPtr = dynamic_cast<DataNode*>(nPtr);
		dataHolder = dPtr->data;
		return true;
	}

	// traverseChild - If the currentChild is a TagNode, set the currentNode to
	// this child (and std::stack the parent and the child std::list iterator)
	//
	bool BasicReader::traverseChild() {

		TagTreeFrame forStack;
		std::list<node*>::iterator i;
		
		if( !haveChildren()) return false;

		forStack.ptr = currentNode;
		forStack.it = currentChild;
		parentStore.push(forStack);

		currentNode = *currentChild;
		if(currentNode->Ident == Tag) {
			TNode = dynamic_cast<TagNode*>(currentNode);
			currentChild = TNode->children.begin();
			DNode = NULL;
		}
		else {
			DNode = dynamic_cast<DataNode*>(currentNode);
			TNode = NULL;
			// don't reset currentChild, it's okay
			currentChild = i;
		}
		return true;
	}


	// upToParent - Goes up one level in the tag tree, taking an only TagTree frame off the std::stack
	// and reseting our node & child members with the contents of the std::stack.  This function returns
	// false if the std::stack is empty
	//
	bool BasicReader::upToParent() {

		TagTreeFrame restoreFrame;
		if( parentStore.empty() ) return false;
		
		restoreFrame = parentStore.top();
		parentStore.pop();
		currentNode = restoreFrame.ptr;
		currentChild = restoreFrame.it;
		if(currentNode->Ident == Tag) {
			TNode = dynamic_cast<TagNode*>(currentNode);
			DNode = NULL;
		}
		else {
			DNode = dynamic_cast<DataNode*>(currentNode);
			TNode = NULL;
		}
		return true;
	}

	// Up - given a number N, perform N upToParent() calls
	// NOTE - Does not return a boolean indicating whether operation
	// was sucessful
	BasicReader &BasicReader::Up(size_t N)
	{
		size_t i;
		for(i = 0; i < N; ++i) upToParent();
		return *this;
	}

	// Down - same as above, only do traverseChild() rather than
	// upToParent
	BasicReader &BasicReader::Down(size_t N)
	{
		size_t i;
		for(i = 0; i < N; ++i) traverseChild();
		return *this;
	}

	// Right - same as above, only do getNextChild() rather than
	// upToParent
	BasicReader &BasicReader::Right(size_t N)
	{
		size_t i;
		for(i = 0; i < N; ++i) getNextChild();
		return *this;
	}

	// Left - blah blah blah getPrevChild
	BasicReader &BasicReader::Left(size_t N)
	 {
		 size_t i;
		 for(i = 0; i < N; ++i) getPrevChild();
		 return *this;
	 }

	  
     
	// BasicReader::char2digit - a simple conversion between a single char and an int value
	// performed by subtracting ASCII value 48, '0', from the character.
	//
	int BasicReader::char2digit(char in) {
		return in - 48;
	}




	// Reader() - default constructor - Just calls the constructor of the BasicReader and that's it
	Reader::Reader() : BasicReader() {}

	// Reader(std::istream &, const char* encoding) - call the BasicReader constructor and do nothing else
	Reader::Reader(std::istream &in, const char* encoding) : BasicReader(in, encoding) {}

	Reader::Reader(const Reader &Rin) : BasicReader(Rin) {} // Call the base class copy constructor

	const Reader &Reader::operator=(const Reader &Rin) {
		BasicReader::operator=(Rin);
		return *this;
	}
	

	// ~Reader - This is just an adapter, do nothing
	// okay, according to my reliable sources, the BasicReader
	// destructor is called automatically, so I need to make no
	// specific reference to it.  If that's not true, I'll alter it
	Reader::~Reader() {}

	// expectTag - checks whether the currentNode is a Tag
	bool Reader::expectTag() {
		std::string error;

		if(!initalized() ) return false;
		if(!isTag()) {
			error = "Expecting the next node to be a Tag, instead is Text.";
			setErrorString(error);
			setErrorCode(Reader::NOT_TAG);
			return false;
		}
		else return true;
	}

	// checkTag - the "check" version of the expect above
	bool Reader::checkTag() const {
		if(!initalized() || !isTag()) return false;
		else return true;
	}

	// expectText - checks whether the currentNode is Text or not
	bool Reader::expectText() {
		std::string error;
		if(!initalized()) return false;
		if(!isText()) {
			error = "Was expecting the next node to be Text, instead is a Tag.";
			setErrorString(error);
			setErrorCode(Reader::NOT_TEXT);
			return false;
		}
		else return true;
	}

	// checkText - the "check" version of the expect above
	bool Reader::checkText() const {
		if(!initalized() || !isText() ) return false;
		else return true;
	}

	
	// expectTagName - Takes is a std::string and checks whether the name of the currentTag is equal to it
	bool Reader::expectTagName(const std::string &name) {

		std::string error, true_name;
		if( !initalized() || !expectTag() ) return false;

		getTagName(true_name);
		if(true_name != name) {
			error = "Expecting Tag Name: \"";
			error += name;
			error += "\", Got: \"";
			error += true_name;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::TAGNAME_ERROR);
			return false;
		}
		return true;

	}

	// checkTagName - Checks for the Tag Name, but doesn't set errors if it isn't the same
	bool Reader::checkTagName(const std::string &name) const {
		std::string true_name;
		if(!initalized() || !checkTag() ) return false;
		getTagName(true_name);
		return true_name == name;
	
	}



	// expectTagName - Takes in a char* and checks whether the name of the currentTag is equal to it
	bool Reader::expectTagName(const char* name) {

		std::string error, true_name;
		if( !initalized() || !expectTag() ) return false;

		getTagName(true_name);
		if( true_name != name) {
			error = "Expecting Tag name: \"";
			error += name;
			error += "\", Got: \"";
			error += true_name;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::TAGNAME_ERROR);
			return false;
		}
		return true;
	}



	// checkTagName - Same thing
	bool Reader::checkTagName(const char* name) const {
		std::string true_name;
		if(!initalized() || !checkTag()) return false;
		getTagName(true_name);

		return true_name == name;
	}

	template<class Num>
	bool Reader::expectTagNum(Num &n) {

		std::string s;

		if(!expectTagName("cn") || !expectChildTextString(s)) return false;

		if(!Reader::isNum(s)) {
			s = "Attempting to convert numerical entry to number, got std::string: " + s;
			setErrorString(s);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		if(!Reader::toNum(n, s)) {
			setErrorString("Attempted to convert number, but overflow occured.");
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}


	template<class Num>
	bool Reader::checkTagNum(Num &N) {
		std::string s;

		if(!checkTagName("cn") || !checkChildTextString(s) || !Reader::isNum(s) || !Reader::toNum(N, s)) return false;

		return true;
	}


	template<class Num>
	bool Reader::expectTagNumVector(std::vector<Num> &v) {
		
		v.clear();
		Num n;

		traverseChild();
		if(!checkTagNum(n)) {
			setErrorString("Tried to convert first entry of numerical std::vector to number, failed");
			setErrorCode(Reader::OTHER);
			upToParent();
			return false;
		}
		v.push_back(n);
		upToParent();
		while(getNextChild()) {
			traverseChild();
			if(!checkTagNum(n)) {
				setErrorString("Tried to convert next entry of numerical std::vector to number, failed");
				setErrorCode(Reader::OTHER);
				
				upToParent();
				Left(v.size());
				return false;
			}
			else
				v.push_back(n);
			upToParent();
		}
		Left(v.size() - 1);

		return true;
	}

	template<class Num>
	bool Reader::checkTagNumVector(std::vector<Num> &v) {

		Num n;
		
		traverseChild();
		if(!checkTagNum(n)) {
			upToParent();
			return false;
		}
		else
			v.push_back(n);
		upToParent();
		while(getNextChild()) {
			traverseChild();
			if(!checkTagNum(n)) {
				upToParent();
				Left(v.size());
				return false;
			}
			else
				v.push_back(n);
			upToParent();
		}

		return true;
	}


	// template version that uses a field
	template<class Field>
	bool Reader::expectFieldNum(const Field &F, typename Field::Element &e) {

		if(!toNum(F, e)) {
			setErrorString("Attempted to convert current Tag to field element, but failed!");
			setErrrorCode(Reader::OTHER);
			return false;
		}

		return true;
	}

       	template<class Field>
	bool Reader::checkFieldNum(const Field &F, typename Field::Element &e) {

		return toNum(f, e);
	}


	template<class Field>
	bool Reader::expectFieldNumVector(const Field &F, std::vector<typename Field::Element>&v) {
		if(!toNumVector(F, v)) {
			setErrorString("Attempted to convert current child and subsequent children to field std::vector, failed");
			setErrorCode(Reader::OTHER);

			return false;
		}
		return true;

	}

	template<class Field>
	bool Reader::checkFieldNumVector(const Field & F, std::vector<typename Field::Element>& v) {

		return toNumVector(F, v);
	}



	// expectAttributeString - Takes in a std::string name and a std::string reference.  If the attribute is there, return it.  Otherwise,
	// set the error conditons.  Note for this function, the type returned from the BasicReader function, a std::string, is already
	// in the proper form, so this function will never set a WRONG_ATT_TYPE error
	bool Reader::expectAttributeString(const std::string &name, std::string &value) {

		std::string error, true_name;
		if(!initalized() || !expectTag() ) return false;

		if( !getAttribValue(name, value)) {
			error = "Could not find Attribute: \"";
			error += name;
			error += "\" in Tag: \"";
			getTagName(true_name);
			error += true_name;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::MISSING_ATT);
			return false;
		}

		return true;
	}


	// checkAttributeString - Just checks for an attribute of the name, doesn't set error codes
	bool Reader::checkAttributeString(const std::string &name, std::string & value) const {

		if(!initalized() || !checkTag() ) return false;
		return getAttribValue(name, value);
	}


	// expectAttributeNum - Takes in a std::string name and an int reference.  If the attribute is there, check to see whether
	// it is a num.  If it is, convert it and return the results.  Otherwise, set the appropriate error
	//
	template<class Num>
	bool Reader::expectAttributeNum(const std::string &name, Num &value) {

		std::string hold, error;
		if( !expectAttributeString(name, hold) ) return false;

		if(!Reader::isNum(hold) ) {
			error = "Found Attribute \"";
			error += name;
			error += "\", but was expecting value to be a numeric, instead got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_ATT_TYPE);
			return false;
		}
		if(!Reader::toNum(value, hold)) {
			error = "Tried to convert std::string \"";
			error += hold;
			error += "\" to numeric, but conversion overflowed.";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}

	template<class Num>
	bool Reader::checkAttributeNum(const std::string &name, Num &value) const {
		std::string hold;
		if( !checkAttributeString(name, hold) || !Reader::isNum(hold)) return false;
	       
		return Reader::toNum(value, hold);
	}


	// depricated due to introduction of templates.  Removed
	// expectAttributeGMP - Takes in a std::string name and an GMP reference.  If the attribute is there, check to see whether
	// it is a num.  If it is, convert it and return the results.  Otherwise, set the appropriate error
	//
	//	bool Reader::expectAttributeGMP(const std::string &name, integer &value) {

	//		std::string hold, error;
	//		if( !expectAttributeString(name, hold) ) return false;

	//		if(!Reader::isNum(hold) ) {
	//			error = "Found Attribute \"";
	//			error += name;
	//			error += "\", but was expecting value to be a GMP, instead got \"";
	//			error += hold;
	//			error += "\".";
	//			setErrorString(error);
	//			setErrorCode(Reader::WRONG_ATT_TYPE);
	//			return false;
	//		}
	//		Reader::toGMP(value, hold);
	//		return true;
	//	}

	//	bool Reader::checkAttributeGMP(const std::string &name, integer &value) const {
	//		std::string hold;
	//		if( !checkAttributeString(name, hold) || !Reader::isNum(hold) ) return false;
		
	//		Reader::toGMP(value, hold);
	//		return true;
	//	}


	// expectAttributeLong - Takes in a std::string name and a long reference.  If the attribute is there, check to see whether
	// it is a num.  If it is, convert it and return the results.  Otherwise, set the appropriate error
	//
	//	bool Reader::expectAttributeLong(const std::string &name, long &value) {
	//
	//		std::string hold, error;
	//		if( !expectAttributeString(name, hold)) return false;

	//		if(!Reader::isNum(hold) ) {
	//			error = "Found Attribute \"";
	//			error += name;
	//			error += "\", but was expecting value to be an long, instead got \"";
	//			error += hold;
	//			error += "\".";
	//			setErrorString(error);
	//			setErrorCode(Reader::WRONG_ATT_TYPE);
	//			return false;
	//		}
	//		Reader::toLong(value, hold);
	//		return true;
	//	}


	//	bool Reader::checkAttributeLong(const std::string &name, long &value) const {
	//		std::string hold;
	//		if(!checkAttributeString(name, hold) || !Reader::isNum(hold)) return false;
	//
	//		Reader::toLong(value, hold);
	//		return true;
	//	}

	// template specialization for size_t
	// expectAttributeNum - Takes in a std::string name and an size_t reference.  If the attribute is there, check to see whether
	// it is a num.  If it is, convert it and return the results.  Otherwise, set the appropriate error
	//
	template<>
	bool Reader::expectAttributeNum(const std::string &name, size_t &value) {

		std::string hold, error;
		if( !expectAttributeString(name, hold)) return false;

		if(!Reader::isUnsignedNum(hold) ) {
			error = "Found Attribute \"";
			error += name;
			error += "\", but was expecting value to be numeric, instead got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_ATT_TYPE);
			return false;
		}
		if(!Reader::toNum(value, hold)) {
			error = "Attempted to convert std::string \"";
			error += hold;
			error += "\" to numeric, but conversion failed.";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}

	template<>
	bool Reader::checkAttributeNum(const std::string &name, size_t &value) const {
		std::string hold;
		if(!checkAttributeString(name, hold) || !Reader::isUnsignedNum(hold) ) return false;

		return Reader::toNum(value, hold);
	}



	// expectNChildren - Takes in a size_t and checks whether the currentNode has that many children.  If so, return true,
	// if not, set error codes & return false
	bool Reader::expectNChildren(const size_t &N) {
	
		std::string error, true_name;
		char buffer[100]; // in case of error message, we need to be able to write N & numChildren to std::string format
		if(!initalized() || !expectTag() ) return false;

		if(N != numChildren() ) {
			error = "Was expecting ";
			sprintf(buffer, "%to", N);
			error += buffer;
			error += " children for Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", instead have ";
			sprintf(buffer, "%to", numChildren());
			error += buffer;
			setErrorString(error);
			setErrorCode(Reader::WRONG_N_CHILD);
			return false;
		}
		return true;
	}



	// expectAtLeastNChildren - Checks that the number of children is greater than or equal to
	// the input N
	bool Reader::expectAtLeastNChildren(const size_t &N) {
		std::string error, true_name;
		char buffer[100]; // in case of error message, we need to be able to write N & numChildren to std::string format
		if(!initalized() || !expectTag() ) return false;

		if(N > numChildren() ) {
			error = "Was expecting at least ";
			sprintf(buffer, "%to", N);
			error += buffer;
			error += " children for Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", instead have ";
			sprintf(buffer, "%to", numChildren());
			error += buffer;
			setErrorString(error);
			setErrorCode(Reader::WRONG_N_CHILD);
			return false;
		}
		return true;
	}


	// expectTextString - Expect the Tag to be Text, and if not return false and write error codes
	bool Reader::expectTextString(std::string &writeable) {
		// we want the error codes from expectText()
		if(!initalized() || !expectText()) return false;

		return getText(writeable); // will be true
	}

	bool Reader::checkTextString(std::string &writeable) const {
		if(!initalized() || !checkText()) return false;

		return getText(writeable); // will be true
	}


	// expectTextNum - Expect The Tag to be Text, and that this Text be convertable to Numeric type
	template<class Num>
	bool Reader::expectTextNum(Num &value) {
		std::string hold, error;
		if(!expectTextString(hold)) return false;

		if(!Reader::isNum(hold)) {
			error = "Was expecting Text to be numeric, instead have: ";
			error += hold;
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}
		if(!Reader::toNum(value, hold)) {
			error = "Tried to convert \"";
			error += hold;
			error += "\" to numeric, but conversion failed.";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}


	// template specialization for size_t
	template<> 
        bool Reader::expectTextNum(size_t& value) {
		std::string hold, error;
		if(!expectTextString(hold)) return false;

		if(!Reader::isUnsignedNum(hold)) {
			error = "Was expecting Text to be unsigned numeric, instead have:";
			error += hold;
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}
		if(!Reader::toNum(value, hold)) {
			error = "Tried to convert std::string \"";
			error += hold;
			error += "\" to numeric, but conversion failed.";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}


	template<class Num>
	bool Reader::checkTextNum(Num &value) const {
		std::string hold;
		if(!checkTextString(hold) || !Reader::isNum(hold)) return false;

		return Reader::toNum(value, hold);
	}
	template<>
	bool Reader::checkTextNum(size_t &value) {
		std::string hold;
		if(!checkTextString(hold) || !Reader::isUnsignedNum(hold)) return false;

		return Reader::toNum(value, hold);
	}


	// expectTextGMP - Expect the Tag to be Text, and that this Text to be convertable to GMP
	//	bool Reader::expectTextGMP(integer &value) {
	//		std::string hold, error;
	//		if(!expectTextString(hold)) return false;
	//
	//		if(!Reader::isNum(hold)) {
	//			error = "Was expecting Text to be a GMP, instead have: ";
	//			error += hold;
	//			setErrorString(error);
	//			setErrorCode(Reader::WRONG_DATA_TYPE);
	//			return false;
	//		}
	//		Reader::toGMP(value, hold);
	//		return true;
	//	}

	// expectTextLong - Expect the Tag to be Text, and that this Text to be convertable to Long
	//	bool Reader::expectTextLong(long &value) {
	//		std::string hold, error;
	//		if(!expectTextString(hold)) return false;
	//
	//		if(!Reader::isNum(hold)) {
	//			error = "Was expecting Text to be a Long, instead have: ";
	//			error += hold;
	//			setErrorString(error);
	//			setErrorCode(Reader::WRONG_DATA_TYPE);
	//			return false;
	//		}
	//		Reader::toLong(value, hold);
	//		return true;
	//	}

	// expectTextSizet - Expect the Tag to be Text, and that this Text convertable to Sizet
	//	bool Reader::expectTextSizet(size_t &value) {
	//		std::string hold, error;
	//		if(!expectTextString(hold)) return false;
	//
	//              if(!Reader::isUnsignedNum(hold)) {
	//			error = "Was expecting Text to be a Size_t, instead have: ";
	//		error += hold;
	//			setErrorString(error);
	//			setErrorCode(Reader::WRONG_DATA_TYPE);
	//			return false;
	//		}

	//		Reader::toSizet(value, hold);
	//		return true;
	//	}

	// expectTextNumVector - Expect the Tag to be Text, and that this Text be convertible to an Num std::vector
	template<class Num>
	bool Reader::expectTextNumVector(std::vector<Num> &value, bool addOne) {
       		std::string hold, error;
		if(!expectTextString(hold)) return false;

		if(!Reader::isNumVector(hold)) {
			error = "Was expecting Text to be an num std::vector, instead have: ";
			error += hold;
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		if(!Reader::toNumVector(value, hold, addOne)) {
			error = "Tried to convert std::string to numeric std::vector, but conversion failed (overflow).";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}

	template<class Num>
	bool Reader::checkTextNumVector(std::vector<Num> &value, bool addOne) const {
		std::string hold;
		if(!checkTextString(hold) || !Reader::isNumVector(hold)) return false;

		return Reader::toNumVector(value, addOne);
	}


	// expectTextGMPVector - Expect the Tag to be Text, and that this Text be convertible to a GMP std::vector
	/*	bool Reader::expectTextGMPVector(std::vector<integer> &value) {
		std::string hold, error;
		if(!expectTextString(hold)) return false;

		if(!Reader::isNumVector(hold)) {
			error = "Was expecting Text to be a GMP std::vector, instead have: ";
			error += hold;
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		Reader::toGMPVector(value, hold);
		return true;
	}

	// expectTextLongVector - Expect the Tag to be Text, and that this Text be convertible to an int std::vector
	bool Reader::expectTextLongVector(std::vector<long> &value) {
		std::string hold, error;
		if(!expectTextString(hold)) return false;

		if(!Reader::isNumVector(hold)) {
			error = "Was expecting Text to be a long std::vector, instead have: ";
			error += hold;
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		Reader::toLongVector(value, hold);
		return true;
	}

	*/
	// template specialization for size_t
	// expectTextNumVector - Expect the Tag to be Text, and that this Text be convertible to a size_t std::vector
	template<>
	bool Reader::expectTextNumVector(std::vector<size_t> &value, bool addOne) {
		std::string hold, error;
		if(!expectTextString(hold)) return false;

		if(!Reader::isUnsignedNumVector(hold)) {
			error = "Was expecting Text to be a numeric std::vector, instead have: ";
			error += hold;
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}
		if(!Reader::toNumVector(value, hold, addOne)) {
			error = "Tried converting std::string to numeric std::vector, but conversion failed (overflow";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}

	template<>
	bool Reader::checkTextNumVector(std::vector<size_t> &value, bool addOne) {
		std::string hold;
		if(!checkTextString(hold) || !Reader::isUnsignedNumVector(hold)) return false;

		return Reader::toNumVector(value, hold, addOne);
	}
	

	// expectChildTextString - Takes in a std::string reference and checks whether the currentChild is a Data Node, and if so
	// initalizes the reference to the returned std::string.  If not, sets the error conditions and returns false;
	// Note, this function won't ever set a WRONG_DATA_TYPE error code, as there is no need to convert the 
	// data given to another format, it's just a std::string
	//
	bool Reader::expectChildTextString(std::string &value) {

		std::string error, true_name;
		if( !initalized() || !expectTag() ) return false;

		if( !isChildText() ) {
			error = "Was expecting current Child of Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\" to be character data.";
			setErrorString(error);
			setErrorCode(Reader::CHILD_NOT_DATA);
			return false;
		}
		getChildText(value);
		return true;
	}

	bool Reader::checkChildTextString(std::string &value) const {
		getChildText(value);
	}


	// expectChildTextNum - Takes in an Num reference and checks whether the currentChild is a Data Node and if so
	// checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
	//
	template<class Num>
	bool Reader::expectChildTextNum(Num &value) {
		
		std::string hold, error, true_name;
		if( !expectChildTextString(hold) ) return false;

		if( !Reader::isNum(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", but was expecting numeric, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		if(!Reader::toNum(value, hold)) {
			error = "Tried to convert std::string \"";
			error += hold;
			error = "\", to numeric, but conversion failed.";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}

	template<class Num>
	bool Reader::checkChildTextNum(Num &value) const {
		std::string hold;
		if(!checkChildTextString(hold) || !Reader::isNum(hold)) return false;

		return Reader::toNum(value, hold);
	}


	// expectChildDataGMP - Takes in a GMP reference and checks whether the currentChild is a Data Node and if so
	// checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
	//
	/*	bool Reader::expectChildTextGMP(integer &value) {
		
		std::string hold, error, true_name;
		if( !expectChildTextString(hold) ) return false;

		if( !Reader::isNum(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", but was expecting GMP, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		Reader::toGMP(value, hold);
		return true;
	}

	bool Reader::checkChildTextGMP(integer &value) const {
		std::string hold;
		if(!checkChildTextString(hold) || !Reader::isNum(hold)) return false;

		Reader::toGMP(value, hold);
		return true;
	}
		
	
	// expectChildDataLong - Takes in an long reference and checks whether the currentChild is a Data Node and if so
	// checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
	//
	bool Reader::expectChildTextLong(long &value) {
		
		std::string hold, error, true_name;
		if( !expectChildTextString(hold)) return false;

		if( !Reader::isNum(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", but was expecting long, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		Reader::toLong(value, hold);
		return true;
	}

	bool Reader::checkChildTextLong(long &value) const {
		std::string hold;
		if( !checkChildTextString(hold) || !Reader::isNum(hold)) return false;

		Reader::toLong(value, hold);
		return true;
	}
	*/

	// template specialization for size_t
	// expectChildTextNum - Takes in an size_t reference and checks whether the currentChild is a Data Node and if so
	// checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
	template<>
	bool Reader::expectChildTextNum(size_t &value) {
		
		std::string hold, error, true_name;
		if( !expectChildTextString(hold)) return false;

		if( !Reader::isUnsignedNum(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", but was expecting numeric, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		
		if(!Reader::toNum(value, hold)) {
			error = "Tried to convert std::string \"";
			error += hold;
			error = "\", but conversion failed.";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}

	template<>
	bool Reader::checkChildTextNum(size_t &value) const {
		std::string hold;
		if(!checkChildTextString(hold) || !Reader::isUnsignedNum(hold)) return false;

		return Reader::toNum(value, hold);
		return true;
	}


	// expectChildDataNumVector - Takes in an Num std::vector reference and checks whether the currentChild is a Data Node and if so
	// checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
	template<class Num>
	bool Reader::expectChildTextNumVector(std::vector<Num> &value, bool addOne) {
		
		std::string hold, error, true_name;
		if( !expectChildTextString(hold) ) return false;

		if( !Reader::isNumVector(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", but was expecting numeric std::vector, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		if(!Reader::toNumVector(value, hold, addOne)) {
			error = "Tried to convert std::string to numeric std::vector, but conversion failed (overflow).";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}


	template<class Num>
	bool Reader::checkChildTextNumVector(std::vector<Num> &value, bool addOne) const {
		std::string hold;
		if(!checkChildTextString(hold) || !Reader::isNumVector(hold)) return false;

		return Reader::toNumVector(value, hold, addOne);
	}

	// expectChildDataGMPVector - Takes in a GMP reference and checks whether the currentChild is a Data Node and if so
	// checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
      	//
	/*	bool Reader::expectChildTextGMPVector(std::vector<integer> &value) {
		
		std::string hold, error, true_name;
		if( !expectChildTextString(hold) ) return false;

		if( !Reader::isNumVector(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(true_name);
			error += true_name;
			error += "\", but was expecting GMP std::vector, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		Reader::toGMPVector(value, hold);
		return true;
	}


	bool Reader::checkChildTextGMPVector(std::vector<integer> &value) const {
		std::string hold;
		if(!checkChildTextString(hold) || !Reader::isNumVector(hold)) return false;

		Reader::toGMPVector(value, hold);
		return true;
	}



	// expectChildDataLongVector - Takes in a long std::vector reference and checks whether the currentChild is a Data Node and if so
	// checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
	//
	bool Reader::expectChildTextLongVector(std::vector<long> &value) {
		
		std::string hold, error, name;
		if( !expectChildTextString(hold) ) return false;

		if( !Reader::isNumVector(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(name);
			error += name;
			error += "\", but was expecting long std::vector, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		Reader::toLongVector(value, hold);
		return true;
	}


	bool Reader::checkChildTextLongVector(std::vector<long> &value) const {
		std::string hold;
		if(!checkChildTextString(hold) || !Reader::isNumVector(hold)) return false;

		Reader::toLongVector(value, hold);
		return true;
	}
	*/

	// template specialization for size_t
	// expectChildDataSizetVector - Takes in an size_t std::vector reference and checks whether the currentChild is a Data Node 
	// and if so checks whether the returns data represents a number.  If it does, initalize the reference to the numerical
	// representation of the data.  Otherwise, set the appropriate error flags
	//
	template<>
	bool Reader::expectChildTextNumVector(std::vector<size_t> &value, bool addOne) {
		
		std::string hold, error, name;
		if( !expectChildTextString(hold)) return false;

		if( !Reader::isUnsignedNumVector(hold) ) {
			error = "Found Child Data of Tag \"";
			getTagName(name);
			error += name;
			error += "\", but was expecting size_t std::vector, got \"";
			error += hold;
			error += "\".";
			setErrorString(error);
			setErrorCode(Reader::WRONG_DATA_TYPE);
			return false;
		}

		if(!Reader::toNumVector(value, hold, addOne)) {
			error = "Tried to convert std::string to numeric std::vector, but conversion failed (overflow).";
			setErrorString(error);
			setErrorCode(Reader::OVFLOW);
			return false;
		}
		return true;
	}

	template<>
	bool Reader::checkChildTextNumVector(std::vector<size_t> &value, bool addOne) const {
		std::string hold;
		if(!checkChildTextString(hold) || !Reader::isUnsignedNumVector(hold)) return false;

		return Reader::toNumVector(value, hold, addOne);
	}

	
	// expectChildTag - Checks whether the current Child is a Tag.  If it isn't, set error conditions & returns false
	// This one is meant to be used in conjunction with traverseChild() and upToParent()
	bool Reader::expectChildTag() {

		std::string error, name;
		if( !initalized() || !expectTag()) return false;

		if( !isChildTag() ) {
			error = "Was expecting current Child of Tag \"";
			getTagName(name);
			error += name;
			error += "\" to be an embeded tag.";
			setErrorString(error);
			setErrorCode(Reader::CHILD_NOT_TAG);
			return false;
		}
		return true;
	}

	bool Reader::checkChildTag() const {
		if(!initalized() || !isChildTag()) return false;
		return true;
	}

	Reader &Reader::Up(size_t N) {
		BasicReader::Up(N);
		return *this;
	}

	Reader &Reader::Down(size_t N) {
		BasicReader::Down(N);
		return *this;
	}

	Reader &Reader::Right(size_t N) {
		BasicReader::Right(N);
		return *this;
	}

	Reader &Reader::Left(size_t N) {
		BasicReader::Left(N);
		return *this;
	}



	/* ReaderIterator - The iterator type of the Reader "container" above.  This object is an Input iterator that
	 * supports operator--() (so that we can go back & forth in the Reader).  The general requirement that I'm making
	 * is that input manipulation of the ReaderIterator will NOT mangle the internal state of the source Reader.
	 * The bad side of this is that if you have exceptional behavior using the Iterator, this does not show up
	 * when using the ReaderIterator, though it'd be nice if it did.  I may tool around w/ Reader to make this possible
	 * As a result, the ReaderIterator, rather than saving a copy of the spawning Reader, actually IS-a Reader (so that
	 * it has it's own internal Reader state variables which are intialized by the source Reader), and it is THESE
	 * state variables that are manipulated.  In addition, the ReaderIterator has a group of flags indicating when it is
	 * "past the end".  Note, if you try to dereference the iterator if it's "past the end", an exception is thrown
	 * (because frankly I can't think of a better way to deal w/ the case where the user has made a ReaderIterator
	 * out of a TextNode, then tries to actually USE that iterator)
	 */
	class ReaderIterator : private Reader {

		typedef Reader    value_type;
		typedef Reader*   pointer;
		typedef Reader&  reference;
		typedef ptrdiff_t difference_type;
		typedef size_t    size_type;

	public:
		ReaderIterator();
		ReaderIterator(const ReaderIterator &);
		ReaderIterator(const Reader &, std::list<node*>::iterator, bool);
		~ReaderIterator();

		bool operator==(const ReaderIterator &) const;
		bool operator!=(const ReaderIterator &) const;

		const ReaderIterator &operator=(const ReaderIterator &);

		reference operator*();
		pointer   operator->();

		ReaderIterator &operator++();
		bool isValidIterator() const;
		
	private:
		bool calledOnText; // in case the user is dumb
		std::list<node*>::iterator eqIter; // used for operator== & operator!=, as well as end()

	};

	// Default constructor.  With no state to point to, this object can iterate nothing, so just use an 
	// uninitalized std::list<node*>::iterator I guess
	ReaderIterator::ReaderIterator() : Reader() {
		calledOnText = true; // not actually true, but true enough to cause the exception to be thrown
		std::list<node*>::iterator it;
		eqIter = it; // get some end iterator
	}

	// Copy constructor.  Makes a copy of the previous iterator
	ReaderIterator::ReaderIterator(const ReaderIterator &x) : Reader(x) { // use polymorphism, call the base copy constr.
		calledOnText = x.calledOnText;
		eqIter = x.eqIter;
	}

	// The "real" constructor used by the Reader class
	ReaderIterator::ReaderIterator(const Reader &x, std::list<node*>::iterator it, bool onText) : Reader(x) { // call the base copy constructor
		eqIter = it;
		calledOnText = onText;
		if(haveChildren()) {
			traverseChild();
		}
			
	}

	// Empty destructor, we have no resources to take care of
	ReaderIterator::~ReaderIterator() {}

	// operator== - part of the input iterator API
	// Note, I am assuming that all "past-end" iterators are equal, which logically makes the
	// "logical-or" work here
	bool ReaderIterator::operator==(const ReaderIterator &rhs) const {
		if(calledOnText && rhs.calledOnText) 
			return true;
		else
			return eqIter == rhs.eqIter;
	}

	// operator!= again part of the input iterator API.  Same as above
	bool ReaderIterator::operator!=(const ReaderIterator &rhs) const {

		return !(*this == rhs);
	}

	const ReaderIterator &ReaderIterator::operator=(const ReaderIterator &rhs) {
		Reader::operator=(rhs);
		calledOnText = rhs.calledOnText;
		eqIter = rhs.eqIter;
		return *this;
	}
		

	// The increment operator.  Increments eqIter
	// don't bother w/ worrying about calledOnText, if the pointer has been mangled
	// than whatever weird results come from these calls won't matter, the flag will
	// still be true and attmepted dereferencing will throw the Error
	ReaderIterator &ReaderIterator::operator++() {
		upToParent();
		getNextChild();
		traverseChild();
		++eqIter;
		return *this;
	}
		
	bool ReaderIterator::isValidIterator() const {
		return !calledOnText;
	}


	// operator* - Get a Reader at the current state level of the Iterator.  Not defined if
	// the iterator in question was initalized using Text Data;
	//
	ReaderIterator::reference ReaderIterator::operator*() {
		if(calledOnText) {
			throw Error("Tried to use a ReaderIterator initalized with XML Text.", Undef_Action);
		}
		else
			return *this;
	}

	// operator-> Return a pointer to a Reader, included for calls like "i->expectTagName("something")", etc
	// Once again, throws an exception if this iterator was initalized using a Reader on Text
	//
	ReaderIterator::pointer ReaderIterator::operator->() {
		if(calledOnText) {
			throw Error("Tried to use a ReaderIterator initalized with XML Text.", Undef_Action);
		}
		else
			return this;
	}


	// Reader "begin" method for Reader iterator.  Creates a new Reader iterator (using of all things a reference
	// to itself.  
	//
	Reader::iterator Reader::begin() {

		std::list<node*>::iterator i;
		if(isTag()) {
			i = TNode->children.begin();
		}
		return ReaderIterator(*this, i, isChildText());
	}

	Reader::const_iterator Reader::begin() const {
		std::list<node*>::iterator i;
		if(isTag()) {
			i = TNode->children.begin();
		}
		return ReaderIterator(*this, i, isChildText());
	}

	// Reader "end" method for Reader iterator.  Creates a ReaderIterator w/ the eqIter value set to the end
	Reader::iterator Reader::end() {
		std::list<node*>::iterator i;
		if(isTag()) {
			i = TNode->children.end();
		}
		return ReaderIterator(*this, i, isChildText());
	}
	
	Reader::const_iterator Reader::end() const {
		std::list<node*>::iterator i;
		if(isTag()) {
			i = TNode->children.end();
		}
		return ReaderIterator(*this, i, isChildText());
	}
			

				
				
} // namespace LinBox

#endif // ifndef __LINBOX_reader_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
