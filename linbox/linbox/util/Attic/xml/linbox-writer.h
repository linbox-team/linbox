#ifndef __LINBOX_WRITER
#define __LINBOX_WRITER

// Important parts of the STL that are needed
#include <iostream>
#include <sstream>
#include <std::stack>
#include <std::list>
#include <std::string>

// For extended Data type
#include "linbox/integer.h"

// for sprintf
#include <cstdio>
// for strlen and strcpy
#include <cstring>

// TagNode definitions
#include "xml-tree.h"

//#include "linbox-reader.h"

// For LinBox integer type
// #include "relative/path/to/linbox/integer.h"

/*

Tag Object Writer - LinBox Writer utility

The Tag Object Writer wraps an XML tag structure TagNode tree or XMLTree,
providing API that simplifies the creation of tag objects


With this utility you can do the following:
- Write the tag object created
write(std::ostream&)

- Specify a tag name:
setTagName(const char*)
setTagName(const std::string&);

- Add attributes
addAttribute(const char*, const char*)
addAttrribute(const std::string&, const std::string&);

- Convert an int or GMP to a char* or std::string
intToString(int)
intToCString(int)
GMPToString(const integer&)
GMPToCString(const integer&)

- Add Data Children
addDataChild(const char*)
addDataChild(const std::string&)
addNumericalList(const vector<int>&) <- For <index>1 2 3 4 1 2 1</index>
addNumericalList(const vector<integer>&)

- Add a sequence promise.  This adds a Promise Child to the 
  TagList's class.  Be sure if you use this method you are using a
  sequence type (vector, std::list, etc) that support iterators

// Note for subsequent calls to these functions, if the previous
// node is a data child, concatenate onto this one.


- Add children:
(Have an internal std::stack which stores parent *ers at previous levels)
addTagChild()
UpToParent()


*/


namespace LinBox 
{


	//	template<class Num>
	//std::string &num2String(std::string &, const Num&);


	struct WriterFrame {
		TagNode* tnode;
		std::list<node*>::iterator iter;
	};

	// forward declaration
	class BasicReader;

  
	class Writer {
		friend class BasicReader;

	public:
		Writer();
		//  Writer(const Writer&);
		~Writer() { rootNode->killself(); }
		
		// API
		bool write(std::ostream &);
		
		void setTagName(const char*);
		void setTagName(const std::string &);
		
		void setAttribute(const char*, const char*);
		void setAttribute(const std::string &, const std::string &);
	
		//		static std::string &numToString<int>(std::string &, const int &);
		
		//		template<>
		//		static std::string &numToString<long>(std::string &, const long &);

		//		template<>
		//		static std::string &numToString<size_t>(std::string &, const size_t &);
		
		//		template<>
		//		static std::string &numToString<integer>(std::string &, const integer &);

		template<class Num>
		static std::string &numToString(std::string &, const Num &);
		

		//		template<>
		//		static char* numToCString<int>(char*, const int &);

		//		template<>
		//		static char* numToCString<long>(char*, const long &);

		//		template<>
		//		static char* numToCString<size_t>(char*, const size_t);

		//		template<>
		//		static char* numToCString(char*, const integer &);

		template<class Num>
		static char* numToCString(char*, const Num &);


		bool addDataChild(const char*);
		bool addDataChild(const std::string &);
		bool insertDataChild(const char*);
		bool insertDataChild(const std::string &);


		//		template<>
		//		bool addNumericalList<num>(const vector<int> &);

		//		template<>
		//		bool addNumericalList<long>(const vector<long> &);

		//		template<>
		//		bool addNumericalList<size_t>(const vector<size_t> &);

		//		template<>
		//		bool addNumericalList<integer>(const vector<integer> &);

		template<class Num>
		bool addNum(const Num &);

		template<class Num>
		bool insertNum(const Num &);

		template<class Num>
		bool addNumericalList(const vector<Num> &, bool = false);
		
		template<class Num>
		bool insertNumericalList(const vector<Num> &, bool = false);

		template<class Field>
		bool addNumericalList(const Field &, const vector<typename Field::Element> &);
		template<class Field>
		bool insertNumericalList(const Field &, const vector<typename Field::Element> &);


		template<class ContainerType>
       		bool addSequencePromise(const ContainerType &);	

		template<class ContainerType>
		bool insertSequencePromise(const ContainerType &);
	
		bool addTagChild();
		bool insertTagChild();

		bool upToParent();

		bool getNextChild();
		bool getPrevChild();

		bool isFirstChild();
		bool isLastChild();

	
	private:

		static char int2char(int);
		
		TagNode* rootNode, *currentNode;
		std::list<node*>::iterator currentChild;
		std::stack<WriterFrame> parentStore;
	
		TagNode* makeNewTagNode();
		bool DataLast(); // a simple predicate that checks whether the current child of the current TagNode is a DataNode
	
	
	};


// This default constructor is the normal constructor, expected to be used
// whenever a user wants to write a LinBox object in the XML format.  Most
// importantly, it sets the parent flag to true, enabling destruction and
// writing
//
	Writer::Writer() {
		rootNode = makeNewTagNode();
		currentNode = rootNode;
		currentChild = currentNode->children.begin();
	}

// write - Takes in an std::ostream and writes the tag structure currntly contained
// if writing is allowed
// pre-condition: std::ostream ready to be written
// post-conditon: tag structure represented by Writer is written to std::ostream,
// if parent flag is set, otherwise false is returned
//
	bool Writer::write(std::ostream & o) {

		o << "<?xml version=\"1.0\"?>" << std::endl;
		rootNode->xml_out(o);
		return true;
	}

// setTagName - Takes a const char* and assigns to tag at the current level
// the name given by the char*
	void Writer::setTagName(const char* newName) {
		currentNode->tag.assign(newName); 
		return;
	}

// setTagName - Takes a std::string and assigns to the tag at the current level
// the name given by the std::string
	void Writer::setTagName(const std::string &newName) {
		currentNode->tag = newName;
		return;
	}

// setAttribute - Takes two char*, the the attribute name and the second it's
// value, and either 1) overwrites the attribute value if the attribute
// is already defined or 2) adds the new attribute if it isn't already defined
//
	void Writer::setAttribute(const char* name, const char* value) {
		
		std::list<std::string>::iterator atPtr;

		for(atPtr = currentNode->attrib.begin(); 
		    atPtr != currentNode->attrib.end() && *atPtr != name; 
		    ++atPtr) ++atPtr;

		if(atPtr == currentNode->attrib.end() ) {
			currentNode->attrib.push_back(std::string(name));
			currentNode->attrib.push_back(std::string(value));
		}
		else {
			++atPtr;
			atPtr->assign(value);
		}

		return;
	}
	
// setAttribute - Takes two strings, the attribute name and the second it's value,
// and either 1) overwrite the attribute value if the attribute is already defined
// or 2) adds the new attribute if it isn't already defined
//
	void Writer::setAttribute(const std::string &name, const std::string & value) {
		
		std::list<std::string>::iterator atPtr;

		for(atPtr = currentNode->attrib.begin();
		    atPtr != currentNode->attrib.end() && *atPtr != name;
		    ++atPtr) ++atPtr;

		if(atPtr == currentNode->attrib.end() ) {
			currentNode->attrib.push_back(name);
			currentNode->attrib.push_back(value);
		}
		else {
			++atPtr;
			*atPtr = value;
		}

		return;
	}


// numToString<int> - A utility for LinBox users.  Takes in an integer and returns a std::string representation
// of that number.  Uses an internal character buffer to convert the int to a std::string using the sprintf
// function from standard C, then returns a std::string w/ this data internally.  Note, uses a 100
// char buffer, even though this is more than necessary, as even a 64 bit int is only about 20 digits
// long. 
	template<>
	std::string &Writer::numToString(std::string &source, const int &theNumber) {

		char buffer[100];
		sprintf(buffer, "%d", theNumber);
		source.assign(buffer);

		return source;
	}

	template<>
	std::string &Writer::numToString(std::string &source, const long &theNumber) {
		
		char buffer[100];
		sprintf(buffer, "%ld", theNumber);
		source.assign(buffer);

		return source;
	}

	template<>
	std::string &Writer::numToString(std::string &source, const size_t &theNumber) {

		char buffer[100];
		sprintf(buffer, "%zd", theNumber);
		source.assign(buffer);

		return source;
	}

	// A generic method for when we have a peice of data represented
	// by a template parameter.  This code is slower but will work
	// with any numerical type that supports integral / & % and comparisons
	// against integers
	//
	// Uses a stringstream and relies on the type in question having
	// a C++ stream insertion operator available (which I can then
	// use.  If no such operator is available, the operation fails
	// Sorry
	//
	template<class Num>
	std::string &Writer::numToString(std::string &source, const Num &theNum) {
		std::ostringstream oss;

		oss << theNum;
		source = oss.str();

		return source;
	}



		//		Num buff;
		//		int x;
		//		if(theNum == 0) {
		//			source = "0";
		//		}
		//		else {
		//			source.assign("");
		//			if(theNum < 0) {
		//				buff = -1 * theNum;
		//			}
		//			else
		//				buff = theNum;
		//			while(buff > 0) {
		//				x = buff % 10;
		//				source.insert(0, 1, Writer::int2char(x));
		//				buff /= 10;
		//			}
		//			if(theNum < 0)
		//				source.insert(0, 1, '-');
		//		}
		//		return source;
		//	}
		

	//	template<class Num>
        //std::string &Writer::numToString(std::string &source, const Num &theNum) {
	//	num2String(source, theNum);
	//}



// intToCString - A utility for LinBox users.  Takes in an int and a char* which is assumed to
// point to a location with the space neccessary to hold the int.  So if something is overwritten,
// it is assumed to be the user's fault :-)
//
	template<>
	char* Writer::numToCString<int>(char* buffer, const int &theNumber ) {
		
		sprintf(buffer, "%d", theNumber);
		return buffer;
	}

	template<>
	char* Writer::numToCString<long>(char* buffer, const long &theNumber) {
		
		sprintf(buffer, "%ld", theNumber);
		return buffer;
	}

	template<>
	char* Writer::numToCString<size_t>(char* buffer, const size_t &theNumber) {
		
		sprintf(buffer, "%zd", theNumber);
		return buffer;
	}

	// GMPToString - Takes a LinBox integer and converts it to an
	// STL std::string using the built in conversion function
	//
	template<>
	std::string &Writer::numToString<integer>(std::string &source, const integer& num) {
		return Integer2string(source, num);
	}

	// GMPTOCString - Takes a char* (which is assumed to be
	// initalized with enough space) and writes the contents of the
	// of GMP integer there.  
	//
	template<>
	char* Writer::numToCString<integer>(char* buffer, const integer& num) {
		std::string helper;
		Integer2string(helper, num);
		strcpy(buffer, helper.c_str());
		return buffer;
	}


	// generic version of the above functions, for when we
	// don't know the type of the number.  Note that it is assumed
	// that the char* given has been initalized to have enough space
	//
	template<class Num>
	char* Writer::numToCString(char* buffer,  const Num& theNum) {

		std::ostringstream oss;

		oss << theNum;
		strcpy(buffer, oss.str().c_str());

		return buffer;
	}


		//		std::list<int> L;
		//		size_t i;
		//		Num buff;
		//
		//		if(theNum == 0) {
		//			strcpy(buffer, "0");
		//		}
		//		else {
		//			if(theNum < 0) {
		//				buff = -1 * theNum;
		//				i = 1;
		//				buffer[0] = '-';
		//			}
		//			else {
		//				buff = theNum;
		//				i = 0;
		//			}
		//			while(buff > 0) {
		//				L.push_back(buff % 10);
		//				buff /= 10;
		//			}
		//			
		//			for(std::list<int>::iterator it = L.begin(); it != L.end(); ++it, ++i) {
		//				buffer[i] = Writer::int2char(*it);
		//			}
		//			buffer[i] = "\0";
		//		}
		//		return buffer;
		//	}
		
				


	// Takes a char* and adds a new DataChild to the XML structure.  Note if a DataChild was the
	// last thing added to the class, this peice of data is just concatenated onto the end of it
	//
	bool Writer::addDataChild(const char* dataToAdd) {
	
		node* lastPtr;
		DataNode* DataNodePtr;
		if( DataLast() ) {
     
		// get the last child and typecast it to a DataNode
			lastPtr = *currentChild;
			DataNodePtr = dynamic_cast<DataNode*>(lastPtr);
			DataNodePtr->data += dataToAdd; // appends this new char data onto the end of the DataNode's data
		}
		else {
			DataNodePtr = new DataNode();
			DataNodePtr->data.assign(dataToAdd);
			DataNodePtr->valid = true;
			++currentChild;
			currentNode->children.insert(currentChild, DataNodePtr);
			--currentChild;
		}

		return true;
	}

	bool Writer::insertDataChild(const char* dataToAdd) {
		node* lastPtr;
		DataNode* DataNodePtr;
		if( DataLast() ) {
			lastPtr = *currentChild;
			DataNodePtr = dynamic_cast<DataNode*>(lastPtr);
			DataNodePtr->data = dataToAdd + DataNodePtr->data;
		}
		else {
			DataNodePtr = new DataNode();
			DataNodePtr->data.assign(dataToAdd);
			DataNodePtr->valid = true;
			currentNode->children.insert(currentChild, DataNodePtr);
			--currentChild;
		}
		return true;
	}



	// Takes a std::string and adds a new DataChild to the XML structure. Note if a DataChild was the
	// last thing added to the class, this peice of data is just concatenated onto the end of it
	//
	bool Writer::addDataChild(const std::string &dataToAdd) {

		node* lastPtr;
		DataNode* DataNodePtr;
		if(DataLast() ) {

		// get the last child and typecase it to a DataNode
			lastPtr = *currentChild;
			DataNodePtr = dynamic_cast<DataNode*>(lastPtr);
			DataNodePtr->data += dataToAdd; // appends this new char data onto the end of the std::string
		}
		else {
			DataNodePtr = new DataNode();
			DataNodePtr->data = dataToAdd;
			DataNodePtr->valid = true;
			++currentChild;
			currentNode->children.insert(currentChild, DataNodePtr);
		        --currentChild;
		}
  
		return true;
	}

	bool Writer::insertDataChild(const std::string &dataToAdd) {
		node* lastPtr;
		DataNode* DataNodePtr;
		if( DataLast() ) {
			lastPtr = *currentChild;
			DataNodePtr = dynamic_cast<DataNode*>(lastPtr);
			DataNodePtr->data = dataToAdd + DataNodePtr->data;
		}
		else {
			DataNodePtr = new DataNode();
			DataNodePtr->data = dataToAdd;
			DataNodePtr->valid = true;
			currentNode->children.insert(currentChild, DataNodePtr);
			--currentChild;
		}
		return true;
	}

	template<class Num>
	bool Writer::addNum(const Num &N) {
		std::string s;

		addTagChild();
		setTagName("cn");
		Writer::numToString(s, N);
		addDataChild(s);
		upToParent();

		return true;
	}

	template<class Num>
	bool Writer::insertNum(const Num &N) {
		std::string s;

		insertTagChild();
		setTagName("cn");
		Writer::numToString(s, N);
		addDataChild(s);
		upToParent();

		return true;
	}
		



// Takes a vector of ints and adds a data child which represents the vector in the following way:
// <1 2 3 4 5 6> -> "1 2 3 4 5 6".  Used for things like "<elements>1 2 3 4. . .</elements>" or
// <index>1 2 3 4. . .</index>
//
	template<class Num>
	bool Writer::addNumericalList(const vector<Num> &numVect, bool subOne) {
		
		typename vector<Num>::const_iterator iter;
		std::string s;

		for(iter = numVect.begin(); iter != numVect.end(); ++iter) {
			addTagChild();
			setTagName("cn");
			if(subOne) 
				Writer::numToString(s, *iter - 1);
			else
				Writer::numToString(s, *iter);

			addDataChild(s);
			upToParent();
		}

		return true;
	}

	template<class Num>
	bool Writer::insertNumericalList(const vector<Num> &numVector, bool subOne) {

		typename vector<Num>::const_iterator iter;
		std::string s;
		
		if(numVector.empty()) return true;

		iter = numVector.begin();
		insertTagChild();
		setTagName("cn");
		if(subOne)
			Writer::numToString(s, *iter - 1);
		else
			Writer::numToString(s, *iter);
		addDataChild(s);
		upToParent();

		++iter;
		while( iter != numVector.end()) {
			addTagChild();
			setTagName("cn");
			if(subOne)
				Writer::numToString(s, *iter - 1);
			else
				Writer::numToString(s, *iter);
			addDataChild(s);
			upToParent();
			++iter;
		}

		return true;
	}




	template<class Field>
	bool Writer::addNumericalList(const Field &F, const vector<typename Field::Element> &v) {

		typedef typename Field::Element Element;
		typename vector<Element>::const_iterator iter;

		for(iter = v.begin(); iter != v.end(); ++iter) {
			addTagChild();
			F.toTag(*this, *iter);
			upToParent();
		}

		return true;
	}

	template<class Field>
	bool Writer::insertNumericalList(const Field &F, const vector<typename Field::Element> &v) {

		typedef typename Field::Element Element;
		typename vector<Element>::const_iterator iter = v.begin();

		if(v.empty()) return;

		insertTagChild();
		F.toTag(*this, *iter);
		upToParent();
		++iter;

		while(iter != v.end()) {
			addTagChild();
			F.toTag(*this, *iter);
			upToParent();
			++iter;
		}

		return true;
	}
			

	// addSequencePromise - A templatized memberfunction, templatized
	// on the Container you are using (it must be a linear sequence
	// of something), takes in one of these sequence types, and
	// adds a PromiseNode.  Note, if you add one of these things,
	// BE REALLY, REALLY, REALLY, REALLY SURE that you don't delete
	//  your data structure before outputing the XML, else the thing
	// won't come out properly (in fact, it could be very likely that
	// you'll get a horrific mess, or an error!)
	// Note, usage of addSequencePromise implies 3 things about the
	// data type you are handing it:  1)  it supports the sequence
	// container interface (more specifically, it has empty(), 
	// forward iterators, and begin() & end()), 2) The elements
	// of this sequence container have a proper operator<< 
	// setup (so if you try to print numbers w/ a call like
	// o << *it, digits are written to the stream), and
	// 3) The container you input persists past whenever you
	// actually call the write function for the writer

	template<class ContainerType>
	bool Writer::addSequencePromise(const ContainerType &cont) {

		PromiseNode<ContainerType>* pn = new PromiseNode<ContainerType>(cont);

		++currentChild;
		currentNode->children.insert(currentChild, pn);
		--currentChild;

		return true;
	}





	// insertSequencePromise - A templatized memberfunction, templatized
	// on the Container you are using (it must be a linear sequence
	// of something), takes in one of these sequence types, and
	// adds a PromiseNode.  Note, if you add one of these things,
	// BE REALLY, REALLY, REALLY, REALLY SURE that you don't delete
	//  your data structure before outputing the XML, else the thing
	// won't come out properly (in fact, it could be very likely that
	// you'll get a horrific mess, or an error!)
	// Note, usage of addSequencePromise implies 3 things about the
	// data type you are handing it:  1)  it supports the sequence
	// container interface (more specifically, it has empty(), 
	// forward iterators, and begin() & end()), 2) The elements
	// of this sequence container have a proper operator<< 
	// setup (so if you try to print numbers w/ a call like
	// o << *it, digits are written to the stream), and
	// 3) The container you input persists past whenever you
	// actually call the write function for the writer

	template<class ContainerType>
	bool Writer::insertSequencePromise(const ContainerType &cont) {

		PromiseNode<ContainerType>* pn = new PromiseNode<ContainerType>(cont);

		currentNode->children.insert(currentChild, pn);
		--currentChild;

		return true;
	}





// addTagChild - Adds a new TagNode as a child to the current TagNode, then makes this new TagNode
// the current TagNode (so the pattern of usage is at any particular level, you say 
//
//  Writer W;
//  W.addTagChild();
//  W.setTagName("name");
//  W.upToParent();
//  W.addTagChild();
//  W.setTagName("super");
//  W.upToParent - Creates a writer w/ an empty tag that has two tag children, "name" and "super", so as is
// this thing will contain:
//
//  <blank>
//     <name></name>
//     <super></super>
//  </blank>
// as blank is the default name of a tag (and we didn't give the root a name)
//
	bool Writer::addTagChild() {
		TagNode* newNode = makeNewTagNode();
		WriterFrame wf;

		++currentChild;
		currentNode->children.insert(currentChild, newNode);
		--currentChild;

		wf.tnode = currentNode;
		wf.iter = currentChild;
		
		parentStore.push(wf);
		currentNode = newNode;
		currentChild = currentNode->children.begin();
		
		return true;
	}


	bool Writer::insertTagChild() {
		TagNode *newNode = makeNewTagNode();
		WriterFrame wf;
		
		currentNode->children.insert(currentChild, newNode);
		--currentChild;

		wf.tnode = currentNode;
		wf.iter = currentChild;

		parentStore.push(wf);
		currentNode = newNode;
		currentChild = currentNode->children.begin();

		return true;
	}



// upToParent - Goes up one level to the parent of this node.  Simply picks the parent off the std::stack
// (if there is a parent to pick).
//

	bool Writer::upToParent() {
		WriterFrame wf;

	
		if(parentStore.empty() ) {
			return false;
		}
		else {
			wf = parentStore.top();
			parentStore.pop();
			currentNode = wf.tnode;
			currentChild = wf.iter;
		}


		return true;
	}

	bool Writer::isFirstChild() {
		if( currentChild == currentNode->children.begin()) 
			return true;
		else
			return false;
	}

	bool Writer::isLastChild() {
		bool check;

		++currentChild;
		check = currentChild == currentNode->children.end();
		--currentChild;

		return check;
	}

	bool Writer::getNextChild() {
		if(!isLastChild()) {
			++currentChild;
			return true;
		}
		else
			return false;
	}

	bool Writer::getPrevChild() {
		if(!isFirstChild()) {
			--currentChild;
			return true;
		}
		else
			return false;
	}




// makeNewTagNode - a utility function that creates a new TagNode pointer which is initalized to
// "blank" w/ no attributes or children
//
	TagNode* Writer::makeNewTagNode() {
		TagNode* newNode = new TagNode();
		newNode->tag.assign("blank");

		return newNode;
	}

// DataLast - a utility predicate that checks whether the last child of the currentNode is in fact
// a DataNode by checking the Ident of that Node
	bool Writer::DataLast() {
		node* lastPtr;
		if(currentNode->children.empty() ) {
			return false;
		}
  
		lastPtr = currentNode->children.back();
		return lastPtr->Ident == Data;
	}

	// converts a until digit to it's ASCII character value
	//
	char Writer::int2char(int x) {
		return (char) x + 48;
	}


}


#endif // #ifndef __LINBOX_WRITER
