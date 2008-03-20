/* Rich Seagraves - xml-tree.h.  A DOM like XML parsing class that wraps
 * around the expat XML parser.  This file contains the declarations and 
 * definitions of many helper classes and types.  The type of interest is 
 * XMLTree.  This is used to parse and std::istream containing XML text into
 * a tree structure, and write a tree structure in XML formatting to an std::ostream
 */

/*
Copyright (c) 2003 Rich Seagraves

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
 
#ifndef XML_TREE
#define XML_TREE

// Standard libs
#include <list> 
#include <map>
#include <stack>
#include <cstring>
#include <iostream>
#include <string>

// Parts of the C++ standard library that are used

// #include <cstdlib> <- In the file describing LinBox::Reader and
// #include <cctype>  <- LinBox::Writer

// #include "linbox/integer.h" <- Place these in the file describing
// using LinBox::integer;      <- the LinBox::Reader and LinBox::Writer

// For the expat parser
#include "expat.h"

// nodeType is an enumeration used to identify the different parts of the node
// heirarchy described below.  These are Node (for node), Tag (for TagNode), 
// Data (for DataNode) and Root (for RootNode).  The node class has an
// Ident member of type nodeType to tell the different peices of the class 
// apart
enum nodeType { Node, Tag, Data, Promise};

// Class node - A mostly virtual base class for the various types of nodes
// that will go into a parsed XML tree.  Defines that every node must have
// a clone function (virtual copy constructor), killself function 
// (destructor-like function that frees all dynamic memory), and xml_out 
// function (the XML output writer for each individual node).

// Also has an Ident member inhereited by every class so that pointers
// to the various node types can be distinguished (if you have a TagNode versus
// a DataNode, which one do you have?)

class node {
 public:
  nodeType Ident;
  node() {}
  virtual node* clone() const { return new node(); }
  virtual void killself() {}
  virtual void xml_out(std::ostream &) const {}
  virtual ~node() {}
};


// Class DataNode - a node type that holds XML PCDATA (or plain character
// data).  Implements the interface given by node.  Essentially a simple
// std::string class that takes expat formatted XML data (a char * and a len)
// and turns this into a null terminated std::string on the heap.

class DataNode : public node {
 public:
  std::string data;

  DataNode();
  DataNode(const char*, int);
  DataNode(const DataNode &);
  ~DataNode() { killself(); }

  const DataNode & operator=(const DataNode &rhs);
  bool operator==(const DataNode &rhs);
  bool operator!=(const DataNode &rhs);

  // Overloaded methods from the node class
  DataNode* clone() const;
  void killself() {}
  void xml_out(std::ostream &) const;

  // The method used for setting the data.  Note, if the int is 0, it means
  // the std::string is null terminated
  void set(const char*, int = 0);

  bool valid;
};


// class PromiseNode - A special node for writing XML vectors that
// respresents a contract between the Writer & the rest of the LinBox
// objects.  When used, the PromiseNode doesn't contain any XML
// text or tags, but instead contains a pointer to a C++ style
// sequence container that support STL iterators.  When writing
// XML, the PromiseNodes use this interface to write the contents of
// the container as a space seperated text std::string.  It is assumed
// that the data contained w/n the container know how to properly 
// print themselves using the operator<<()
//
template<class ContainerType>
class PromiseNode : public node {

 public:
	PromiseNode();
	PromiseNode(const ContainerType &Ct);
	~PromiseNode();
	PromiseNode(const PromiseNode<ContainerType> &);

	const PromiseNode<ContainerType> &operator=(const PromiseNode<ContainerType>&);
	bool operator==(const PromiseNode<ContainerType> &) const;
	bool operator!=(const PromiseNode<ContainerType> &) const;

	bool valid() const;
	
	// overload methods of Node class
	PromiseNode<ContainerType>* clone() const;
	void killself();
	void xml_out(std::ostream &) const;

 private:
	const ContainerType* containerPointer;
	bool initFlag;
};





// Class TagNode - A type of node representing a tag in XML format
// Has a tag, std::list of attributes, and linked std::list of children
// Each TagNode is essentially a subtree.  Fulfills the node functionality,
// but adds some additional helper functions, such as addTagChild and 
// addDataChild (which add a node of the specified type to the TagNodes
// child std::list)

class TagNode : public node {
 public:

  std::string tag;
  std::list<std::string> attrib;
  std::list<node*> children;

  TagNode();
  TagNode(const char*, const char**);
  TagNode(const TagNode &);
  ~TagNode();

  const TagNode & operator=(const TagNode &);
  bool operator==(const TagNode &);
  bool operator!=(const TagNode &);
  
  TagNode* clone() const;
  void killself();
  void addDataChild(const char*, int);
  TagNode* addTagChild(const char*, const char**);
  void xml_out(std::ostream &) const;
  void init(const char*, const char**);

};


/*
// Class RootNode - A special type of TagNode representing the root
// element of the XML document.  No current use yet, but there are
// some features allowing the RootNode to easily be upgraded to give
// it some details about the entire XML document.  The stacks defined below
// are not used.  In addition, none of the node interface functions are
// defined here (so that RootNode uses TagNode's implentations for these).
// This may change


class RootNode : public TagNode {
 public:
  std::stack<node *> entries;
  std::stack<int> attribs;
  
  RootNode();
  RootNode(const char*, const char**);

};
*/

// Handler functions used by the expat parser when parsing XML documents
// start is called when expat encounters an open tag, data when expat
// encounters PCDATA (raw character data inside a tag), and finish when
// expat encounters a close tag.  Are declared as friends to the 
// XMLTree class below.  

void start(void*, const char*, const char**);
void data(void*, const char*, int);
void finish(void*, const char*);
// int convert_int(const char*); <- In the file containing LinBox::Reader and
// integer convert_integer(const char*); <- LinBox::Writer

// ErrorType is an enumeration of the possible Errors that can crop up
// when reading/writing XML data.  Either a)  an XML document couldn't be
// parsed 
// for lack of well-formedness b) We had trouble converting from character
// data to integral data or c) the tree structure you are trying to write
// to XML is broken.  This ErrorType is used in the Error struct below

enum ErrorType {XML_ParseError, Invalid_XML, Data_ConvError, UnInit, Undef_Action, NoError};

// Error is a struct thrown when an error occurs in the function of
// the parser.  User code using the XMLTree class should try to catch
// the Error class, as it can be thrown.

struct Error
{
  const char* error_descrip;
  ErrorType EType;
  int param;
  Error(const char* desc, ErrorType e, int para = 0) 
  {
    error_descrip = desc;
    EType = e;
    param = para;
  }
};


// Class XMLTree - a DOM like XML reading/writing utility that wraps
// the expat parser.  Declares the start, data and finish classes as friends
// To read in an XML document and parse it, call parse with an std::istream.
// Note that it will read from the beginning of the stream until it encounters
// an EOF line.  It also ignores all new-line and tab characters (to ensure
// they don't get in the way).  It can be accessed using the Tree element,
// which is the Root XML element of the document.  It also has a LabelTable
// for reference labels, which haven't been implemented yet.

// To write an XMLTree into an XML document, call the write method with
// the std::ostream you want written to.  Won't write unless the initFlag is set,
// which is only set if the XMLTree has a fully parsed XML document inside.

class XMLTree {
  friend void start(void*, const char*, const char**);
  friend void data(void*, const char*, int);
  friend void finish(void*, const char*);

 public:
  XMLTree();
  XMLTree(std::istream &in, const char* encoding = "US-ASCII");
  ~XMLTree();
  void parse(std::istream &in, const char* encoding = "US-ASCII");
  void write(std::ostream &);

  TagNode Tree;
  std::map<std::string, TagNode*> LabelTable;
  bool initFlag;

 private:
  XML_Parser p; // the actual expat parser
  std::stack<TagNode*> currLeaf; // a std::stack useful in parsing an xml document

};

// The main constructor for the TagNode class
// Most of it's functionality was export to the init function so that
// I could pull a trick with the RootNode of the XMLTree.
// This function is normally called when a new TagNode is created
// from an opened XML tag.  The init function performs all the memory
// allocation and std::string copying.  

TagNode::TagNode()
{
  Ident = Tag;
}

TagNode::TagNode(const char *TagName, const char **Attrib) 
{
  init(TagName, Attrib);
  Ident = Tag;
}


// Real copy constructor.  Neat, huh :-).  Won't suffer from the siamese twin
// problem.  Note, I don't use init here because the init function is made
// to be used with data recieved straight from the expat parser.  As the
// data is saved in the TagNode class in a slightly different manner, for
// these copying functions an "own rolled" version seemed better

TagNode::TagNode(const TagNode &In) {

  // Copy the tagname
  tag = In.tag;

  // Copy the attributes.  Calls the std::list equal operator, which I believe calls
  // the = operator on each of the data members.  As DataNode's operator= are
  // already overloaded, this is okay
  attrib = In.attrib;

  // Make a copy of the source TagNode's children std::list
  // Must specifically copy their data, not just their pointer location, so
  // call the clone function.
  // Note, clone puts these nodes on the heap, so these will be cleaned up
  // by calling each child's killself method

  for(std::list<node*>::const_iterator vi = In.children.begin(); vi != In.children.end(); ++vi) {
    children.push_back( (*vi)->clone() );
  }

  Ident = Tag;

}

// assignment operator.  Makes a copy of a TagNode, without suffering from
// the siamese twin problem.  Calls the clone function of each of rhs's
// children for it's own children std::list

const TagNode & TagNode::operator=( const TagNode &rhs) {

  if( *this != rhs) { // To prevent bone-headed uses like X = X
    killself(); // Deletes all previously saved data for this TagNode

    tag = rhs.tag;
    attrib = rhs.attrib;

    
    // Clone the sources children.  
    for( std::list<node*>::const_iterator vi = rhs.children.begin(); vi != rhs.children.end(); ++vi) {
      children.push_back( (*vi)->clone() );
    }
  }
  return *this; // For X = Y = Z chaining
}

// The init function.  Should be called only on TagNodes whose tag and attrib
// fields have not been initalized, otherwise this function creates a memory
// leak.

void TagNode::init(const char* TagName, const char** Attrib)
{
  int i;

  tag.assign(TagName);

 
  // Pushes back as many DataNodes as needed to hold the attributes
  for(i= 0; Attrib[i] ; ++i) {
    attrib.push_back(std::string(Attrib[i]));
  }
  

}  



/* Note, this function is incomplete (ie, a hack).  
 * Test before continued use
 */

// An equality operator.  Checks whether two TagNodes are the same

bool TagNode::operator==(const TagNode &rhs) {
  return this == &rhs;
}

/* Note, this function is incomplete (ie, a hack).
 * Test before continued use
 */

// TagNode inequality operator.  Checks whether two TagNodes are not
// the same

bool TagNode::operator!=(const TagNode &rhs) {
   return this != &rhs;
}


// Virtual destructor.  Frees all memory used by this tag and calls the
// virtual destructor of all children below this point.

void TagNode::killself() 
{

  // Nothing needs to be done for the attributes, the std::list class destructor
  // will call the destructor on each DataNode, which takes care of internal
  // memory

  // Ensures all children are deleted
  for(std::list<node*>::iterator it = children.begin(); it != children.end(); ++it) {
    (*it)->killself();
  }
  children.clear(); // clear the child std::list

}

// Real destructor.  Calls the virtual destructor.
// Note that this real destructor / virtual destructor business is a misnomer
// To stay standard, the destructor for the node class is virtual, even though
// I like to think of killself as the "virtual" desctructor

TagNode::~TagNode() 
{
  killself();
}

// addDataChild - Utility function for Tag class.  Creates a new data
// node and puts it at the end of the children std::list for this TagNode.
// Created from data supplied by the expat "CharacterDataHandler" function
// The new data child is on the heap and will be freed by this Tag's
// killself function

void TagNode::addDataChild(const char* text, int len) {
  DataNode* newData = new DataNode(text, len);

  if(newData->valid) children.push_back(newData); // if the new DataNode
						// is valid, add it as a child
						// Otherwise, just kill it
  else delete newData;
  return;
}

// addTagChild - Utility function for Tag class.  Creates a new tag node
// and puts it at the end of the children std::list.  Note this lower tag node
// is on the heap and will be deleted later by this tag node's killself
// function
// Note, takes input straight from XML parser

TagNode* TagNode::addTagChild( const char* Tag, const char** Attrib) {
  TagNode* newTag = new TagNode(Tag, Attrib);

  children.push_back( newTag );
  return newTag;
}

// Virtual copy constructor.  Returns a copy of the tagnode on the heap
// Is a dynamic copy, so if you use it, you've got to clean it up

TagNode* TagNode::clone() const 
{
  return new TagNode(*this);
}

// Virtual xml_output function.  Called recursively by the tag above
// Is used to write out this portion of XML data.  Calls the xml_output
// function on all of it's children, if there are any
// Oh yeah, outputs this tag as well-formed XML, indented according to depth

void TagNode::xml_out(std::ostream & out) const
{
  std::list<node*>::const_iterator li;
  std::list<std::string>::const_iterator di;

  //  for(i = 0; i < depth; ++i) out << "\t"; // write the correct indentation
  out << "<" << tag; // Write the start of the tag

  for(di = attrib.begin(); di != attrib.end(); di++) { // Write attributes
    out << " " << *di << " = \"";
    di++;
    out << *di << "\"";
  }
  
  if( !children.empty() ) { // If there are children, close the tag
    out << ">";
    for(li = children.begin(); li != children.end(); ++li) {
      (*li)->xml_out(out); // and write each child at a greater depth
    }
    //    for(i = 0; i < depth; ++i) out << "\t"; // Correct indentation
    out << "</" << tag << ">";
  }
  else out << "/>";
  return;
}

// Default constructor for DataNode.  Ensures element is properly
// identified
DataNode::DataNode()
{
  Ident = Data;
  valid = false;
}

   
// Main constructor for the DataNode class.  Takes plain character data
// provided by the expat parser and transforms this into a null terminated
// std::string

// Note that the validity of a DataNode is determined by whether or not
// it is solely constructed of whitespace (in this case '\t' and '\n'
// characters).  If it is, it is pretty much worthless, so don't copy it
// into the XMLTree 

DataNode::DataNode(const char* text, int len) {

  Ident = Data;
  valid = false;

  // If we've maybe got something to set, otherwise, don't bother
  if(len > 0) set(text, len);

}

// This function will destroy any data already inside, just to let you know
// Notice that if len = 0, it is assumed that not only is the std::string 
// null-terminated, but that it is also valid (ie - not just a std::string of
// white-space). Also, chop and eliminate all leading white space by
// using the std::string's find & substring functions to get a substring
// that eliminates leading whitespace
//
void DataNode::set(const char* text, int len)
{
  int i, j;
  killself(); // First delete all old data

  if(len == 0) { // uh-oh, we've got a null-terminated std::string
    data.assign(text);
    valid = true;
  }
  else {
    // now check the std::string for validity (it isn't all white-space)
    j = 0;
    for(i = 0; i < len; ++i) {
      if( text[i] != '\n' && text[i] != '\t') {
	++j;
      }
    }
    
    // Now check if the std::string is valid, if so, write it
    if(j == 0) {
      valid = false;
      return;
    }
    
    // otherwise, write the thing.  Overwrite the std::string w/
    // a substring that starts at the end of leading whitespace
    data.assign(text, len);
    i = data.find_first_not_of(" \t\n");
    j = data.find_last_not_of(" \t\n");
    data = data.substr(i, j - i + 1);
    valid = true;

  }
  return;
}


// Copy constructor.  Just copies another DataNode

DataNode::DataNode(const DataNode &In) {

  Ident = Data; 
  if(!In.valid) {
    valid = false;
  }
  else {
    valid = true;
    data = In.data;
  }
}

/*
// "virtual" destructor.  Deletes all data used by the class
void DataNode::killself() {}

// The real destructor.  Do I have to be any more explicit
DataNode::~DataNode() { killself() }
*/

// Assignment operator.  Destroys previous data and 
// copies source data
const DataNode & DataNode::operator=(const DataNode &rhs) {

  if( this != &rhs) {
    if(rhs.valid) data = rhs.data; // if rhs is a valid std::string
    else {
      
      valid = false; // and invalidate it
    }
  }

  return *this;
}

// equality operator
bool DataNode::operator==(const DataNode &rhs) {
  return data == rhs.data;
}

// inequality operator
bool DataNode::operator!=(const DataNode &rhs) {
  return data != rhs.data;
}

// Virtual copy constructor.  Returns a copy of this DataNode on the heap
// Any users of this function are required to clean up after it
DataNode* DataNode::clone() const
{
  return new DataNode(*this);
}

// Virtual xml_output function.  Writes the character data in a well-formed
// idented manner.  Indents according to depth.  Writes to the std::ostream supplied

void DataNode::xml_out(std::ostream &out) const
{
	//  int i;
  if( !valid) return; // If nothings here, then return

  //  for(i = 0; i < depth; ++i) out << "\t";
  out << data;// << endl;
  return;
}


// PromiseNode - default constructor, essentially set the pointer to
// null, the flag to false
//
template<class ContainerType>
PromiseNode<ContainerType>::PromiseNode() : containerPointer(NULL), initFlag(false) {
	Ident = Promise;

}

// PromiseNode - Main constructor, take in a container we will
// print to XML and save a pointer to that container.  That's it
//
template<class ContainerType>
PromiseNode<ContainerType>::PromiseNode(const ContainerType &Ct) 
	: containerPointer(&Ct), initFlag(true) {
	Ident = Promise;
}

// PromiseNode - Copy constructor.  Just copy the pointer.  I know,
// this creates two objects pointing to the same data structure, but
// ultimately, this is what a PromiseNode is supposed to do, it just points
// to a data structure that we will eventually convert to XML.  The
// direct data copy is okay, as the state of one will automatically be
// propogated to the other
//
template<class ContainerType>
PromiseNode<ContainerType>::PromiseNode(const PromiseNode<ContainerType> &other) 
	: containerPointer(other.containerPointer), initFlag(other.initFlag) {
	Ident = Promise;

} 


// Destructor.  This one just calls killself, which doesn't do anything,
// so the destructor doesn't really do anything
//
template<class ContainerType>
PromiseNode<ContainerType>::~PromiseNode() { killself(); }

// operator=.  Once again, copy the state of one PromiseNode to the
// other so that the two point to the same object
//
template<class ContainerType>
const PromiseNode<ContainerType> &PromiseNode<ContainerType>::operator=(const PromiseNode<ContainerType> &rhs) {

	containerPointer = rhs.containerPointer;
	initFlag = rhs.initFlag;

	return *this;
}

// operator==.  This one is fairly tricky.  You could have two PromiseNodes
// pointing to the same object, or you could have two PromiseNodes pointing
// to two objects that are the same.  For now, I'm going to assume that
// two PromiseNodes are equal if they point to the same object (as 
// PromiseNodes are going to be used mainly for storing HUGE objects, so
// doing an element by element equality test would be a HUGE time consuming
// operation
//
template<class ContainerType>
bool PromiseNode<ContainerType>::operator==(const PromiseNode<ContainerType> &rhs) const {
	
	return containerPointer == rhs.containerPointer;
}

// operator!= - same story as above
//
template<class ContainerType>
bool PromiseNode<ContainerType>::operator!=(const PromiseNode<ContainerType> &rhs) const {
	
	return containerPointer != rhs.containerPointer;
}

// on to the Node functions . . .

// clone - Returns a dynamically allocated copy of the current node.
// This one just copies the pointer, nothing else.  As w/ all clones, you
// have to clean up for them yourself, that isn't our job
//
template<class ContainerType>
PromiseNode<ContainerType>* PromiseNode<ContainerType>::clone() const {

	return new PromiseNode<ContainerType>(*this);
}

// killself - Deletes any dynamically allocated data and return the node
// to a pre-allocation state.  The best this class can do is merely
// return the containerPointer to null
//
template<class ContainerType>
void PromiseNode<ContainerType>::killself() {
	containerPointer = NULL;
	initFlag = false;
}

// xml_out - the interesting one, the whole reason this class is here.
// This function uses the iterator interface we are assuming ContainerType
// has, and am using this interface to print each element of the
// ContainerType item by item, in a space seperated fashion
//
template<class ContainerType>
void PromiseNode<ContainerType>::xml_out(std::ostream &o) const {

	if(!valid() || containerPointer->empty()) return;
	
	typename ContainerType::const_iterator it = containerPointer->begin();

	o << *it;
	++it;

	while(it != containerPointer->end()) {
		o << " " << *it;
		++it;
	}

	return;
}

// valid - just checks whether the initFlag is set or not
//
template<class ContainerType>
bool PromiseNode<ContainerType>::valid() const {
	return initFlag;
}
							      

/*
// Ensures root has the proper identity

RootNode::RootNode() {
  Ident = Root;
  // Insert other code about special Root related data structure

}

// RootNode constructor for when the RootNode is created during the
// XMl parsing session

RootNode::RootNode(const char* TagName, const char** Attribs)
{
  init(TagName, Attribs);
  Ident = Root;
  // Insert other lines about special Root related data structure
}
*/

// Default constructor for XMLTree.  Creates the parser and ensures the
// initalization flag is false
XMLTree::XMLTree()
{
  p = XML_ParserCreate(NULL);
  initFlag = false;
}

XMLTree::XMLTree(std::istream &in, const char* encoding) 
{
  p = XML_ParserCreate(NULL);
  initFlag = false;
  parse(in, encoding);
}


// Destructor.  Frees all parser memory.

XMLTree::~XMLTree()
{
  Tree.killself();
  XML_ParserFree(p);
}  
  

/** XMLTree::parse - Reads characters off the instream one at a time and parses
 *  until it reachs the End of File mark.  Throws an Error
 *  if there is an error in parsing.  If not, sets the initFlag flag
 */

// Rich Seagraves: 7-2-2003
// Note to Self (or later maintainer, how are you doing?  Hope my code isn't
// giving you too much trouble :-) ): This method has been extended to
// allow multiple XML documents to occur in the same std::istream.  The
// code starting at line 823 replaced this:

/* 
   c = In.get();
   while(c != EOF) {
      buffer.push_back(c):
      c = In.get();
   }
*/


void XMLTree::parse(std::istream &In, const char* encoding)
{
  std::string buffer;
  char c;
  size_t tcount = 0;
  bool done = false;

  // First dump and reset all the data in each data structure
  while( !currLeaf.empty()) currLeaf.pop();
  LabelTable.clear();
  Tree.killself();
  
  XML_ParserReset(p, encoding);
  XML_SetElementHandler(p, start, finish);
  XML_SetCharacterDataHandler(p, data);
  XML_SetUserData(p, this);
 

  c = In.get();
  // skip inital whitespace to get parser off my back about
  // <?xml version="1.0"?> tag
  while(c == ' ' || c == '\n' || c == '\t') c = In.get();
  while(c != EOF && !done) {

	  buffer.push_back(c);
	  if(c == '<') {
		  c = In.get();
		  // processing instruction, do nothing
		  if(c == '?') {
			  buffer.push_back(c);
		  }

		  // CDATA or comment
		  else if(c == '!') {
			  c = In.get();
			  
			  // CDATA
			  if(c == '[') {
				  buffer.push_back(c);
				  while(c != '>') {
					  c = In.get();
					  // prevent infinite loops
					  if(c == EOF) break;
					  else  buffer.push_back(c);
				  }
			  }
			  // comment
			  else {
				  if(c != EOF) buffer.push_back(c);
			  }
		  }
		  // close tag
		  else if(c == '/') {
			  buffer.push_back(c);
			  --tcount;
			  if(tcount == 0) {
				  done = true;
			  }
			  while(c != '>') {
				  c = In.get();
				  if(c == EOF) break;
				  buffer.push_back(c);
			  }
		  }
		  // open tag
		  else {
			  buffer.push_back(c);
			  ++tcount;
			  while(c != '>' && c != '/') {
				  c = In.get();
				  if(c == EOF) break;
				  buffer.push_back(c);
			  }
			  if(c == '/') {
				  --tcount;
				  if(tcount == 0) {
					  done = true;
				  }
				  c = In.get();
				  if(c != '>') {
					  done = true; // parse error
				  }
				  else
					  buffer.push_back(c);
			  }
		  }
	  }
	  if(!done) {
		  c = In.get();
	  }
  }


  if(!XML_Parse(p, buffer.c_str(), buffer.length(), 1)) { 
	  initFlag = false;
	  throw Error(XML_ErrorString(XML_GetErrorCode(p)), XML_ParseError, XML_GetCurrentLineNumber(p));
  } 


  initFlag = true;
  return;
}

// XMLTree write function.  Takes a reference to an std::ostream, and writes 
// a properly formatted XML block


void XMLTree::write(std::ostream &out)
{
  // Throw an error if trying to write to an un-initalized block
  if( !initFlag) throw Error("Trying to write un-initialized data structure.  Nothing written.", UnInit);
  Tree.xml_out(out); //, 0); // Call the xml_out method to write
                        // valid XML to the std::ostream
}



/*
  The following block of code is an iterative version of the code that 
  is currently implemented above.  I opted for the recursive version because
  it eliminates the reliance on RTTI to work properly.  The write
  function for each individual class in the representation tree is a
  virtual function.  That way, if we add more classes into the node
  hierarchy, I don't have to do a complicated re-write of this function.

  I'm saving this in case I have to go back (and because I thought it was
  cool :-) ).  


bool XMLTree::write(std::ostream &out) {
  std::stack< std::pair<node*, std::list<node*>::iterator> > theStack;
  std::list<node*>::iterator li;
  std::pair<node*, std::list<node*>::iterator> holder;
  node * currPtr;
  TagNode* tagPtr;
  DataNode* dataPtr;
  int i, depth;
  bool newFlag = false;

  if( !initFlag) return false;
  holder.first = &Tree;
  holder.second = Tree.children.begin();
  theStack.push(holder);
  currPtr = &Tree;
  depth = -1;

  while( !theStack.empty() ) { // Go until there's nothing left on the std::stack

    if(newFlag) { // If we've encountered a child
      switch(currPtr->Ident) { // Check what the child is
      
      case Data: // We have data
	dataPtr = dynamic_cast<DataNode*>(currPtr);
	for(i = 0; i < depth; ++i) out << "\t";
	out << dataPtr->data << endl; // Just print the data
	newFlag = false;
	depth = depth - 1; // re-adjust depth
	break;
      
      case Tag: // We have a new tag
	tagPtr = dynamic_cast<TagNode*>(currPtr);
	for(i = 0; i < depth; i++) out << "\t"; // Proper indent
	out << "<" << tagPtr->tag; // print the tag opener
	for(std::list<DataNode>::iterator di = tagPtr->attrib.begin(); di != tagPtr->attrib.end(); ++di) {
	   out << " " << di->data << " = \"";
	   ++di;
	   out << di->data << "\""; 
	}
        if( !tagPtr->children.empty() ) { // check whether the tag has children
          out << ">" << endl; // If so, print the > tag ender
          newFlag = true;     // Signal next time that there are children
          depth = depth + 1;  // Adjust depth
          li = tagPtr->children.begin(); // Get the std::list iterator
          holder.first = currPtr;
          currPtr = *li;
          holder.second = ++li;  // Save the current child to the std::stack, and set the new child
          theStack.push(holder);  // as the current pointer
        }
        else {  // There are no children, this is an open and close tag
          newFlag = false; // No more children
          depth = depth - 1; // subtract the depth
          out << "/>" << endl; // print the start and end closer for tag
        }
        break;
      } // close the switch
    } // close the up-most if
    else {
      holder = theStack.top(); // Get the outer most tag
      theStack.pop(); // Get this off so we can update/remove it
      tagPtr = dynamic_cast<TagNode*>(holder.first); // Cast it properly
      li = holder.second;
      if( li != tagPtr->children.end() ) { // There are more children to process
        depth = depth + 1; // Reset the variables for a new child
        newFlag = true;    
        currPtr = *li;     // Get the pointer to the new child
        holder.second = ++li; // Save the new iterator position
        theStack.push(holder);
      }
      else {
	if(tagPtr->Ident != Root) { // If we aren't at the rootnode
	  for(i = 0; i < depth; i++) out << "\t"; // print the end tag
	  out << "</" << tagPtr->tag << ">" << endl;
	}
	depth = depth - 1; // Close up shop
	newFlag = false;
      }
    }
  }

  return true;
}
  
*/


// Handler function for open tags.  Called by expat when it opens valid
// XML tags.  Checks to see whether the std::stack is empty (so we need to setup
// the root node and put it on the std::stack), otherwise pull the top node off
// the std::stack, add the current tag as a child to that TagNode, and put
// the new TagNode at the top of the std::stack.  Do other Root update stuff
  
void start(void* dataforme, const char *tagname, const char** Attributes)
{
  XMLTree *ParseData = (XMLTree*) dataforme;
  TagNode *currNode, *newNode;

  if(ParseData->currLeaf.empty() ) { // We are at the root element
    ParseData->Tree.init(tagname, Attributes); // Setup root's TagNode data
    newNode = & ParseData->Tree; // prepare to put it on the std::stack
    // Other RootNode setup data here
  }
  else {
    currNode = ParseData->currLeaf.top(); // Get the topmost node
    newNode = currNode->addTagChild(tagname, Attributes); // add the new node

    for(std::list<std::string>::iterator di = newNode->attrib.begin(); di != newNode->attrib.end(); ++di) {
      if( (*di) == "label") {
	++di;
	ParseData->LabelTable.insert(std::pair<std::string,TagNode*>(*di, newNode));
      }
      else ++di;
    }
  }
 

  ParseData->currLeaf.push(newNode);
  return;
}

// Handler for character data.  Called by expat when it encounters character
// data in an XML document.  In this case, adds a DataNode to the current
// TagNode's children.

void data(void* dataforme, const char *text, int len)
{

  XMLTree* ParseData = (XMLTree*) dataforme;
  TagNode* currNode = ParseData->currLeaf.top();
  currNode->addDataChild(text, len);

  /* Do stuff with RootNode data here
   *
   *
   *
   */

 return;
}

// Handler function for when a tag closes.  Just pops the current
// TagNode off the std::stack.  Also does something with the RootNode, when
// the RootNode data structure is more focused.

void finish(void* dataforme, const char *tagname)
{

  XMLTree* ParseData = (XMLTree*) dataforme;
  
  /* Do stuff with RootNode data here
   *
   *
   *
   *
   *
   */
  ParseData->currLeaf.pop(); // Pop off the top-most tag, it's done for now
  return;
}

#endif



