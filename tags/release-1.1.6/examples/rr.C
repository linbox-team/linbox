/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#include <iostream>

#include "linbox/integer.h"
#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/rational-reconstruction-base.h"

#include "linbox/util/matrix-stream.h"
#include "linbox/util/timer.h"

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

	UserTimer t;
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
		

} 
