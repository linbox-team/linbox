/* -*- mode: c; style: linux -*- */

/* linbox/tests/test-field-common.h
 * Copyright (C) 2001, 2002 David Saunders
 * Use subject to GNU LGPL.  See linbox/COPYING for details.
 */

#ifndef __TEST_FIELD_COMMON_H
#define __TEST_FIELD_COMMON_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include "test-common.h"

using namespace LinBox;


/* Check the field functionality */

bool isPower(integer n, integer m) // true if n = m^e for some e.
{	return (n == 1) || ((n % m) == 0) && isPower(n/m, m); } 

template<class Field>
bool test_field(Field& F, size_t k, ostream& report, int iters) 
{
	typename Field::element zero, one, two, three;
	typename Field::element a, b, c, d, e, f;

	integer n, m;
	bool pass = true;

	F.characteristic(n);
	F.cardinality(m);
	if ( n > 0 && !isPower(m, n) )
	{	pass = false; 
	report << "characteristic, cardinality mismatch" << endl;
	}

/* tests for presence of members with minimal check of semantics */
// these checks need improvement 
	F.init(zero, 0);
	F.init(one, 1);
	if ( !F.isZero(zero) || F.isZero(one) )
	{	pass = false; report << "isZero problem" << endl; }
	if ( F.isOne(zero) || !F.isOne(one) )
	{	pass = false; report << "isOne problem" << endl; }

	F.init(two, 2);
	F.init(three, 3);

	F.init(a, 0);
	F.init(b, n-2);
	F.convert(m, b);
	if ( m != n - 2 )
	{	pass = false; report << "init/convert problem" << endl; }

	F.init(d, n-2);
	F.init(e, 3);

	F.add(a, three, two); 
	F.assign(d, three); F.addin(d, two);
	if ( !F.areEqual(a, F.init(f, 5)) || !F.areEqual(d, a) )
	{	pass = false; report << "add problem" << endl; }

	F.neg(a, two); 
	F.assign(d, two); F.negin(d);
	if ( !F.areEqual(a, F.init(f, -2)) || !F.areEqual(d, a) )
	{	pass = false; report << "neg problem" << endl; }

	F.sub(a, three, two); F.init(d, 3); F.subin(d, two);
	if ( !F.areEqual(a, one) || !F.areEqual(d, a) )
	{	pass = false; report << "sub problem" << endl; }

	F.mul(a, two, three); F.assign(d, two); F.mulin(d, three);
	F.init(f, 6);
	if ( !F.areEqual(a, f) || !F.areEqual(d, a) )
	{	pass = false; report << "mul problem" << endl; }

	F.inv(a, one); F.assign(d, one); F.invin(d);
	if ( !F.areEqual(a, one) || !F.areEqual(d, a) )
	{	pass = false; report << "inv problem" << endl; }

	F.div(a, two, two); F.assign(d, three); F.divin(d, three);
	if ( !F.areEqual(a, one) || !F.areEqual(d, a) )
	{	pass = false; report << "div problem" << endl; }

	F.axpy(a, two, three, two); 
	F.assign(d, two); F.axpyin(d, two, three);
	if ( !F.areEqual(a, F.init(f, 8)) || !F.areEqual(d, a) )
	{	pass = false; report << "axpy problem" << endl; }

	//,..
	// 2^101 - 1 vs 1 + 2 + 4 + ... + 2^100
	F.init(a, 1);
	F.init(b, 2);
	F.init(c, 0);

	for(int i = 1; i <= 101; ++i)
	{	F.addin(c, a);
	F.mulin(a, b);
	}
	F.subin(a, F.init(f, 1));
	if ( !F.areEqual(a, c) )
	{	pass = false; report << "power problem" << endl; }

	F.assign(d, one);
	for(int i = 1; i < 101; ++i)
		F.axpy(d, two, d, one);
	if ( !F.areEqual(a, d) )
	{	pass = false; report << "axpy power problem" << endl; }

	// there is an extra char in the output - bds 3/02
	report << "field self description: " << F.write(report) << endl;
	report << "field element 2: " << F.write(report, two) << endl;

	/* untested so far
	   ostream &write (ostream &os) const 
	   istream &read (istream &is)
	   ostream &write (ostream &os, const element &x) const 
	   istream &read (istream &is, element &x) const
	   Field_archetype (Field_abstract*, Element_abstract*, RandIter_abstract* = 0)
	*/

	return pass;
}

/****************
class Field_archetype
{
    public:

	typedef Element_archetype element;
	typedef RandIter_archetype RandIter;
    
	Field_archetype (const Field_archetype &F) 
	~Field_archetype (void) 
	Field_archetype &operator=(const Field_archetype &F)

	element &init (element &x, const integer &y = 0 ) const
	integer &convert (integer &x, const element &y = 0) const
	element &assign (element &x, const element &y) const

	integer &cardinality (integer &c) const 
	integer &characteristic (integer &c) const

	bool areEqual (const element &x, const element &y) const
	bool isZero (const element &x) const 
	bool isOne (const element &x) const 

	element &add (element &x, const element &y, const element &z) const
	element &sub (element &x, const element &y, const element &z) const
	element &mul (element &x, const element &y, const element &z) const
	element &div (element &x, const element &y, const element &z) const
	element &neg (element &x, const element &y) const
	element &inv (element &x, const element &y) const
	element &axpy (element       &r, 

	element &addin (element &x, const element &y) const
	element &subin (element &x, const element &y) const
	element &mulin (element &x, const element &y) const
	element &divin (element &x, const element &y) const
	element &negin (element &x) const
	element &invin (element &x) const
	element &axpyin (element &r, const element &a, const element &x) const
	ostream &write (ostream &os) const 
	istream &read (istream &is)
	ostream &write (ostream &os, const element &x) const 
	istream &read (istream &is, element &x) const
	Field_archetype (Field_abstract*, Element_abstract*, RandIter_abstract* = 0)
	template<class Field_qcq> Field_archetype (Field_qcq *f) 
****************/

#endif // __TEST_FIELD_COMMON_H
