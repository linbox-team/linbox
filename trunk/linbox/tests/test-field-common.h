/* -*- mode: c; style: linux -*- */

/* linbox/tests/test-common.C
 * Copyright (C) 2001, 2002 David Saunders
 * Use subject to GNU LGPL.  See linbox/COPYRIGHT for details.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include "test-common.h"

using namespace LinBox;


/* Check the field functionality */

template<class Field>
bool test_field(Field& F, size_t k, ostream& report, int iters) 
{
/* tests for mere presence of members */
	typename Field::element a, b, c, d, e, f;

	integer n, m;
	bool pass = true;

	F.characteristic(n);
	F.cardinality(m);
	if (n > 0 && (m < n || m % n != 0)) 
	{	pass = false; 
		report << "characteristic, cardinality mismatch" << endl;
	}

	F.init(a, 0);
	F.init(b, 2);
	F.init(c, 3);
	F.init(d, 2);
	F.init(e, 3);
	F.init(f, 5);
	F.add(a, b, c); F.addin(d, e);
	if ( !F.areEqual(a, f) || !F.areEqual(d, f) )
	      {pass = false; report << "add problem" << endl;}

	F.neg(a, b); F.negin(d);
	if ( !F.areEqual(a, F.init(f, -2)) || !F.areEqual(d, F.init(f, -5)) )
		{pass = false; report << "neg problem" << endl;}

	F.sub(a, b, c); F.init(d, 2); F.subin(d, e);
	if ( !F.areEqual(a, F.init(f, -1)) || !F.areEqual(d, F.init(f, -1)) )
		{pass = false; report << "sub problem" << endl;}

	F.mul(a, b, c); F.init(d, 2); F.mulin(d, e);
	if ( !F.areEqual(a, F.init(f, 6)) || !F.areEqual(d, F.init(f, 6)) )
		{pass = false; report << "mul problem" << endl;}

	F.inv(a, F.init(f, 1)); F.init(d, 1); F.invin(d);
	if ( !F.areEqual(a, F.init(f, 1)) || !F.areEqual(d, F.init(f, 1)) )
		{pass = false; report << "inv problem" << endl;}

	F.div(a, b, b); F.init(d, 3); F.divin(d, e);
	if ( !F.areEqual(a, F.init(f, 1)) || !F.areEqual(d, F.init(f, 1)) )
		{pass = false; report << "div problem" << endl;}

	// axpy

	//,..
	return pass;
}

