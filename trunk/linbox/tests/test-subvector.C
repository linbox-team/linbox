/* -*- mode: c; style: linux -*- */

/* tests/test-subvector.C
 * Evolved from Will Turner's test-subvector.cpp  -bds
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream.h>
#include <fstream>
#include <vector>
#include "test-common.h"

#include <linbox/field/unparametric.h>

#include <linbox/vector/subvector.h>

/* Test: Subvector class - has the vector interface less those that
 * can invalidate iterators.
 *
 * Return true on success and false on failure
 */

// for testing purposes this vector class is defined
namespace LinBox {
template <class Element>
class CopyReportingVector: public std::vector<Element> {
    public:
	CopyReportingVector(int n): std::vector<Element>(n), copy(false){};
	CopyReportingVector(CopyReportingVector& V): 
		std::vector<Element>(V), copy(true){};
	bool copy;
};
// Vector traits for CopyReportingVector wrapper
template <class Element> 
struct VectorTraits< CopyReportingVector<Element> >
{
	typedef typename VectorTraits<std::vector<Element> >::VectorCategory 
		VectorCategory;
};
}// namespace LinBox

using namespace LinBox;

template <class Field>
static bool testSubvector(Field &F, size_t n, ostream &report) 
{
	report << "Testing subvector class:" << endl;

	bool ret = true;

	typedef typename Field::Element Element;
	typedef std::vector<Element>	Vector;
	//typedef CopyReportingVector<Element>	Vector;
	typedef typename Vector::iterator       Iter;
	typedef Subiterator<typename Vector::iterator>	Subiter;
	typedef typename LinBox::Subvector<Vector, Subiter>	Subvect;
	typedef typename Subvect::const_iterator	ConstSubiterator;

	typedef typename Subvect::reverse_iterator	ReverseIterator;
	typedef typename Subvect::reverse_iterator	ReverseSubiterator;

	Vector v(n);
	for (int z = 0; z < n; z++) v[z] = z;

	printVector(F, report, v);

	int start = 1;
	int stride = 2;
	int length = 3;

	//if (v.copy) report << "vector copy problem" << endl;
	Subvect w(v, start, stride, length);
	//if (v.copy) {report << "vector copy problem" << endl;
	//	     ret = false; return ret;}
	
	//report << "start = " << w._start << endl;
	//report << "stride = " << w._stride << endl;
	//report << "length = " << w._length << endl;

	report << endl << "*** Testing forward iterator" << endl << endl;
	
	Iter i = v.begin();

	Subiter j = w.begin();

//	report << "stride = " << j._stride << endl;

	report << "*j = 1 = " << *j << endl;
	ret = ret && 1 == *j;

	report << "*j++ = 1 = " << *j++ << endl;
	j--;
	ret = ret && 1 == *j++;

	report << "*j = 3 = " << *j << endl;
	ret = ret && 3 == *j;

	report << "*++j = 5 = " << *++j << endl;
	--j;
	ret = ret && 5 == *++j;

	report << "*j = 5 = " << *j << endl;
	ret = ret && 5 == *j;

	report << "*j-- = 5 = " << *j-- << endl;
	++j;
	ret = ret && 5 == *j--;

	report << "*j = 3 = " << *j << endl;
	ret = ret && 3 == *j;

	report << "*--j = 1 = " << *--j << endl;
	++j;
	ret = ret && 1 == *--j;

	report << "*j = 1 = " << *j << endl;
	ret = ret && 1 == *j;

	ret = ret && (j == w.begin());
	if (j != w.begin())
		report << "not at beginning" << endl;
        //j = j + 3;
        j += 3;
	if (j != w.end())
		report << "not at end" << endl;

	report << "Iterating through vector: (";
	for (j = w.begin(); j != w.end(); j++)
	{
		report << *j;
		if ( j < (w.end() - 1) ) report << ", ";
	}
	report << ')' << endl;

	j = w.begin();	
	report << "Random access through vector: (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		report << j[i];
		if ( i < (w.size() - 1) ) report << ", ";
	}
	report << ')' << endl;

#if 1
	report << endl << "*** Testing reverse iterator" << endl << endl;
	
	ReverseSubiterator jr(w.rbegin());

	report << "*jr = 5 = " << *jr << endl;
	ret = ret && 5 == *jr;

	report << "*jr++ = 5 = " << *jr++ << endl;
	--jr;
	ret = ret && 5 == *jr++;

	report << "*jr = 3 = " << *jr << endl;
	ret = ret && 3 == *jr;

	report << "*++jr = 1 = " << *++jr << endl;
	--jr;
	ret = ret && 1 == *++jr;

	report << "*jr = 1 = " << *jr << endl;
	ret = ret && 1 == *jr;

	report << "*jr-- = 1 = " << *jr-- << endl;
	++jr;
	ret = ret && 1 == *jr--;

	report << "*jr = 3 = " << *jr << endl;
	ret = ret && 3 == *jr;

	report << "*--jr 5 = = " << *--jr << endl;
	++jr;
	ret = ret && 5 == *--jr;

	report << "*jr = 5 = " << *jr << endl;
	ret = ret && 5 == *jr;

	if (jr != w.rbegin()) 
		report << "not at beginning" << endl;
	ret = ret && 5 == *jr;

	report << "Iterating through vector: (";
	for (jr = w.rbegin(); jr != w.rend(); jr++)
	{
		report << *jr;
		if ( jr < (w.rend() - 1) ) report << ", ";
	}
	report << ')' << endl;
	ret = ret && jr == w.rend();

	jr = w.rbegin();	
	report << "Random access through vector: (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		report << jr[i];
		if ( i < (w.size() - 1) ) report << ", ";
	}
	report << ')' << endl;

#endif

	report << endl << "*** Testing vector" << endl << endl;

	// Need to check on const_iterator variants of functions.
	
	report << "w.size() = 3 = " << w.size() << endl;
	ret = ret && w.size() == 3;
	report << "w.max_size() = 3 = " << w.max_size() << endl;
	ret = ret && w.max_size() == 3;

	report << "w.front() = 1 = " << w.front() << endl;
	ret = ret && w.front() == 1;
	report << "w.back() = 5 = " << w.back() << endl;
	ret = ret && w.back() == 5;

	report << "w.empty() = false = " << w.empty() << endl;
	ret = ret && !w.empty();

	printVector(F, report, w);

	report << "Printing using operator[]: (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		report << w[i];
		if ( i < (w.size() - 1) ) report << ", ";
	}
	report << ')' << endl;

#if 0
	report << "Printing using at(): (";
	for (unsigned long i = 0; i < w.size(); i++)
	{
		report << w.at(i);
		if ( i < (w.size() - 1) ) report << ", ";
	}
	report << ')' << endl;
#endif

	report << "Assigning subvector: ";
	Vector vv(n);
	Subvect ww(vv, 0, 1, 3);
	//vector<int> ww(3, 77);
	w = ww;
	printVector(F, report, ww);
#if 0
	report << "Constructing subvector from iterators: ";
	Subvect www(w.begin(), w.end());
	printVector(F, report, www);
#endif

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();
	return ret;
}

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 8;

	static Argument args[] = {
		{ 'n', "-n N", "Set size of vector to N (default 8)"}
	};

	parseArguments (argc, argv, args);
	typedef LinBox::UnparametricField<int> Field;
	Field F;

	cout << "Subvector test suite" << endl << endl;

	pass = pass && testSubvector<Field> (F, n, report);

	return pass ? 0 : -1;
}
