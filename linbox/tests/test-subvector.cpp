#include <iostream>
#include <vector>
#include "tests/test-common.h"
#include "linbox/field/unparametric.h"

#include "subvector.h"

int main(void)
{
	typedef LinBox::UnparametricField<int> Field;
	typedef Field::Element Element;
	typedef std::vector<Element>	Vector;
	typedef Vector::iterator	Iterator;
	typedef LinBox::Subvector<Vector>	Subvector;
	typedef Subvector::iterator	Subiterator;
	typedef Subvector::const_iterator	ConstSubiterator;

	typedef Vector::reverse_iterator	ReverseIterator;
	typedef Subvector::reverse_iterator	ReverseSubiterator;

	Field F;

	int n = 8;
	Vector v(8);
	for (int i = 0; i < n; i++) v[i] = i;

	printVector(F, cout, v);

	int start = 1;
	int stride = 2;
	int length = 3;

	Subvector w(v, start, stride, length);
	
	cout << "start = " << w._start << endl;
	cout << "stride = " << w._stride << endl;
	cout << "length = " << w._length << endl;

#if 1
	cout << endl << "*** Testing forward iterator" << endl << endl;
	
	Iterator i(v.begin());
	Subiterator j(w.begin());

//	cout << "stride = " << j._stride << endl;

	cout << "*i = " << *i << endl;
	cout << "*j = " << *j << endl;

	cout << "*i++ = " << *i++ << endl;
	cout << "*j++ = " << *j++ << endl;

	cout << "*i = " << *i << endl;
	cout << "*j = " << *j << endl;

	cout << "*++i = " << *++i << endl;
	cout << "*++j = " << *++j << endl;

	cout << "*i = " << *i << endl;
	cout << "*j = " << *j << endl;

	cout << "*i-- = " << *i-- << endl;
	cout << "*j-- = " << *j-- << endl;

	cout << "*i = " << *i << endl;
	cout << "*j = " << *j << endl;

	cout << "*--i = " << *--i << endl;
	cout << "*--j = " << *--j << endl;

	cout << "*i = " << *i << endl;
	cout << "*j = " << *j << endl;

	if (j == w.begin()) 
		cout << "at begining" << endl;
	else
		cout << "not at beginning" << endl;

	cout << "Iterating through vector: (";
	for (j = w.begin(); j != w.end(); j++)
	{
		cout << *j;
		if ( j < (w.end() - 1) ) cout << ", ";
	}
	cout << ')' << endl;

	j = w.begin();	
	cout << "Random access through vector: (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		cout << j[i];
		if ( i < (w.size() - 1) ) cout << ", ";
	}
	cout << ')' << endl;

#endif

#if 1
	cout << endl << "*** Testing reverse iterator" << endl << endl;
	
	ReverseIterator ir(v.rbegin());
	ReverseSubiterator jr(w.rbegin());

	cout << "*ir = " << *ir << endl;
	cout << "*jr = " << *jr << endl;

	cout << "*ir++ = " << *ir++ << endl;
	cout << "*jr++ = " << *jr++ << endl;

	cout << "*ir = " << *ir << endl;
	cout << "*jr = " << *jr << endl;

	cout << "*++ir = " << *++ir << endl;
	cout << "*++jr = " << *++jr << endl;

	cout << "*ir = " << *ir << endl;
	cout << "*jr = " << *jr << endl;

	cout << "*ir-- = " << *ir-- << endl;
	cout << "*jr-- = " << *jr-- << endl;

	cout << "*ir = " << *ir << endl;
	cout << "*jr = " << *jr << endl;

	cout << "*--ir = " << *--ir << endl;
	cout << "*--jr = " << *--jr << endl;

	cout << "*ir = " << *ir << endl;
	cout << "*jr = " << *jr << endl;

	if (jr == w.rbegin()) 
		cout << "at begining" << endl;
	else
		cout << "not at beginning" << endl;

	cout << "Iterating through vector: (";
	for (jr = w.rbegin(); jr != w.rend(); jr++)
	{
		cout << *jr;
		if ( jr < (w.rend() - 1) ) cout << ", ";
	}
	cout << ')' << endl;

	jr = w.rbegin();	
	cout << "Random access through vector: (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		cout << jr[i];
		if ( i < (w.size() - 1) ) cout << ", ";
	}
	cout << ')' << endl;

#endif

	cout << endl << "*** Testing vector" << endl << endl;

	// Need to check on const_iterator variants of functions.
	
	cout << "w.size() = " << w.size() << endl;
	cout << "w.max_size() = " << w.size() << endl;
	cout << "w.capacity() = " << w.size() << endl;

	cout << "w.front() = " << w.front() << endl;
	cout << "w.back() = " << w.back() << endl;

	printVector(F, cout, w);

	cout << "Printing using operator[]: (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		cout << w[i];
		if ( i < (w.size() - 1) ) cout << ", ";
	}
	cout << ')' << endl;

#if 0
	cout << "Printing using at(): (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		cout << w.at(i);
		if ( i < (w.size() - 1) ) cout << ", ";
	}
	cout << ')' << endl;
#endif

	cout << "Copying subvector: ";
	Subvector ww(w);
	printVector(F, cout, ww);

	cout << "Constructing subvector from iterators: ";
	Subvector www(w.begin(), w.end());
	printVector(F, cout, www);

	
	
	
} // int main(void)
