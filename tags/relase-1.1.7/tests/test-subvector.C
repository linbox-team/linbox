/* Copyright (C) LinBox
 *
 * Evolved from Will Turner's test-subvector.cpp  -bds
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

#include "linbox/linbox-config.h"

#include "linbox/util/commentator.h"
#include "linbox/field/unparametric.h"
#include "linbox/vector/subvector.h"

#include "test-common.h"

using namespace LinBox;

template <class Field>
static bool testSubvector(Field &F, size_t n); 

int main(int argc, char** argv)
{	
    // set up command line options
    static size_t n = 8;
    static Argument args[] = 
    {
 		{ 'n', "-n N", "Set size of vector to N.", TYPE_INT, &n},
		{ '\0' }
    };
    parseArguments (argc, argv, args);
 
    // start testing
	commentator.start("Subvector test suite", "Subvector");
    bool pass = true;
 
    // call tests
    typedef LinBox::UnparametricField<int> Field;
    Field F;
    pass = testSubvector<Field> (F, n);
 
    // finish
	commentator.stop("Subvector test suite");
    return pass? 0 : -1;
}

/* Test Subvector class 
 * Subvector has the vector interface less those that
 * can invalidate iterators.
 */

using namespace LinBox;

template <class Field>
static bool testSubvector(Field &F, size_t n) 
{
	// commentator setup
	const char *  title = "Subvector test";
	commentator.start(title, title, 1);
	ostream &report = commentator.report 
		(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "This test currently neglects several members including constructors." << endl;
	bool ret = true;

	typedef typename Field::Element Element;
	typedef std::vector<Element>	Vector;
	typedef typename Vector::iterator       Iter;
	typedef Subiterator<typename Vector::iterator>	Subiter;
	typedef typename LinBox::Subvector<Subiter>	Subvect;
	//typedef typename LinBox::Subvector<Vector, Subiter>	Subvect;
	typedef typename Subvect::const_iterator	ConstSubiterator;

	typedef typename Subvect::reverse_iterator	ReverseIterator;
	typedef typename Subvect::reverse_iterator	ReverseSubiterator;

	Vector v(n);
	for (size_t z = 0; z < n; z++) v[z] = z;

	printVector(F, report, v);

	int start = 1;
	int stride = 2;
	int length = 3;
	Subiter sb(v.begin()+start, stride);
	Subiter se(sb+length);

//	Subvect w(v, start, stride, length);
//	Subvect w(sb, length);
	Subvect w(sb, se);

	// implicit (not stored) stride of 1
	Subvector<typename Vector::iterator> 
		z(v.begin(), v.end()); 
	// fixme: at least constructor compiles.

	// explicit (stored) stride of 1
	Subvect zz(v.begin(), v.end()); 
	// fixme: at least constructor compiles.
	
	//report << "start = " << w._start << endl;
	//report << "stride = " << w._stride << endl;
	//report << "length = " << w._length << endl;

	report << endl << "*** Testing forward iterator" << endl << endl;
	
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

	// finish
	commentator.stop (MSG_STATUS (ret), (const char *) 0, title);
	return ret;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
