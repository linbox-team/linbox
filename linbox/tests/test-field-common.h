/* -*- mode: c; style: linux -*- */

/* linbox/tests/test-field-common.h
 * Copyright (C) 2001, 2002 David Saunders
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Dave Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __TEST_FIELD_COMMON_H
#define __TEST_FIELD_COMMON_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include "linbox/util/commentator.h"
#include "test-common.h"

using namespace LinBox;

bool isPower (integer n, integer m)
{
	return (n == 1) || ((n % m) == 0) && isPower (n/m, m);
} 

/* Test 1: Test of field operations
 *
 * Test various field oeprations
 *
 * F - Field over which to perform computations
 * title - String to use as the descriptive title of this test
 *
 * Return true on success and false on failure
 */

template<class Field>
bool testField (Field &F, char *title) 
{
	typename Field::element zero, one, two, three;
	typename Field::element a, b, c, d, e, f;

	commentator.start (title, "testField", 5);

	// there is an extra char in the output - bds 3/02
	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Field self description: " << F.write (report) << endl;
	commentator.indent (report);
	report << "field element 2: " << F.write (report, two) << endl;

	integer n, m;
	bool pass = true, part_pass;

	commentator.start ("Testing characteristic/cardinality match");
	part_pass = true;

	F.characteristic (n);
	F.cardinality (m);

	if (n > 0 && !isPower (m, n)) {
		pass = part_pass = false; 
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Characteristic, cardinality mismatch" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	/* tests for presence of members with minimal check of semantics */
	// these checks need improvement 

	commentator.start ("Testing correctness of 0 and 1");
	part_pass = true;

	F.init (zero, 0);
	F.init (one, 1);

	if (!F.isZero (zero)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isZero (0) is false" << endl;
	}

	if (F.isZero (one)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isZero (1) is true" << endl;
	}

	if (F.isOne (zero)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isOne (0) is true" << endl;
	}

	if (!F.isOne (one)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isOne (1) is false" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing init/convert");
	part_pass = true;

	F.init (a, n);
	F.convert (m, a);

	if (m != n) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: F.convert (m, F.init (a, n)) != n" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing field arithmetic");
	part_pass = true;

	F.init (two, 2);
	F.init (three, 3);

	F.init (b, n-2);
	F.init (d, n-2);
	F.init (e, 3);

	F.add (a, three, two); 
	F.assign (d, three);
	F.addin (d, two);

	if (!F.areEqual (a, F.init (f, 5)) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of add are incorrect" << endl;
	}

	F.neg (a, two); 
	F.assign (d, two);
	F.negin (d);

	if (!F.areEqual (a, F.init (f, -2)) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of neg are incorrect" << endl;
	}

	F.sub (a, three, two);
	F.init (d, 3);
	F.subin (d, two);

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of sub are incorrect" << endl;
	}

	F.mul (a, two, three);
	F.assign (d, two);
	F.mulin (d, three);
	F.init (f, 6);

	if (!F.areEqual (a, f) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of mul are incorrect" << endl;
	}

	F.inv (a, one);
	F.assign (d, one);
	F.invin (d);

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of inv are incorrect" << endl;
	}

	F.div (a, two, two);
	F.assign (d, three);
	F.divin (d, three);

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of div are incorrect" << endl;
	}

	F.axpy (a, two, three, two); 
	F.assign (d, two);
	F.axpyin (d, two, three);

	if ( !F.areEqual (a, F.init (f, 8)) || !F.areEqual (d, a) ) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Resutls of axpy are incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing summation of powers of 2");

	//,..
	// 2^101 - 1 vs 1 + 2 + 4 + ... + 2^100

	F.init (a, 1);
	F.init (b, 2);
	F.init (c, 0);

	for (int i = 1; i <= 101; ++i) {
		F.addin (c, a);
		F.mulin (a, b);
	}

	F.subin (a, F.init (f, 1));

	if (!F.areEqual (a, c)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results are incorrect using mul/add" << endl;
	}

	F.assign (d, one);
	for (int i = 1; i < 101; ++i)
		F.axpy (d, two, d, one);

	if (!F.areEqual (a, d)) {
		pass = part_pass = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results are incorrect using axpy" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	/* untested so far
	   ostream &write (ostream &os) const 
	   istream &read (istream &is)
	   ostream &write (ostream &os, const element &x) const 
	   istream &read (istream &is, element &x) const
	   Field_archetype (Field_abstract*, Element_abstract*, RandIter_abstract* = 0)
	*/

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "testField");

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
