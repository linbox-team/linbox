
/*
 * examples/rr.C
 *
 * Copyright (C) 2008, 2010 A. Urbanska
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file examples/rr.C
 * @example  examples/rr.C
 * @ingroup examples
 * @brief Rational Reconstruction.
 */

#include <iostream>

#include <linbox/integer.h>
#include <linbox/field/PID-integer.h>
#include <linbox/algorithms/rational-reconstruction-base.h>

#include <linbox/util/matrix-stream.h>
#include <linbox/util/timer.h>

using namespace LinBox;
using namespace std;

typedef PID_integer Ints;
//typedef Ints::Element Integer;

int main (int argc, char **argv)
{
	srand48( time(NULL) );

	integer a,b;
	integer x(1);
	integer m(1);

	x = x*rand()*rand()*rand()*rand()*rand()*rand();
	m = m*rand()*rand()*rand()*rand()*rand()*rand()*rand()*rand();

	if (x > m) {integer t=x; x=m; m=t;}

	cout << "Searching a,b: a =" << x << "b mod " << m << endl << flush;
	cout << "size of m" << m.bitsize() << endl << flush;

	Ints Z;
	WangClassicRationalReconstruction<Ints> RRB(Z,true,false);
	//RationalReconstruction<Ints, WangClassicRationalReconstruction<Ints> > RR(RRB);
	//RationalReconstruction<Ints, MaxQClassicRationalReconstruction<Ints> > RR(Z);
	//RationalReconstruction<Ints, WangFastRationalReconstruction<Ints> > RR(Z);
	RationalReconstruction<Ints, MaxQFastRationalReconstruction<Ints> > RR(Z);

	LinBox::UserTimer t;
	t.clear();
	t.start();

	for (int i=0; i < 1 ; ++i) {
		if (RR.reconstructRational(a,b,x,m,5)) {
			cout << "Found a,b: "<< a <<"=" << x << "x" << b << " mod " << m << endl << flush;
			cout << "Does agree with bounds\n";
		}
		else {
			cout << "Found a,b: "<< a <<"=" << x << "x" << b << " mod " << m << endl << flush;
			cout << "Does not agree with bounds\n";
		}
	}
	t.stop();
	cout << "Time:";
	t.print(cout);
	cout << endl;

	return 0 ;

}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

