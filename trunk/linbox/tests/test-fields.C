/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-fields.C
 * Written by Dan Roche
 * Copyright (C) June 2004 Dan Roche
 *
 * See COPYING for license information
 */

// to compile: g++ -I../ -I./ -I../.. -I/usr/local/algebra/include -L/usr/local/algebra/lib -lgmp -O3 -Wall -o test-fields test-fields.C ../linbox/.libs/liblinbox.a -lntl

#include "linbox/util/timer.h"
#include "linbox/field/double-fmod.h"
#include "linbox/field/ntl-zz_p.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-int.h"
#include "linbox/field/modular-double.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setw;

// Namespace in which all LinBox code resides
namespace LinBox {

template< class Element >
Element& noop( Element& a, Element& b, Element& c ) { return a; }

template< class Element >
Element& noop( Element& a, Element& b ) { return a; }

template< class Element >
Element& noop( Element& a, Element& b, Element&c, Element& d ) { return a; }

/* fieldTest is a template function to test out the performance of a given field on a
 * machine.  Taken are three arguments.  The first is a field class object.  The second
 * is an array, declared but not necessarily initialized, of seven doubles.  These
 * values will be filled with timings for add, sub, neg, mul, int, div, and axpy,
 * respectively.  The third argument is optional and specifies how many loop
 * iterations to use.
 */

template< class Field >
void fieldTest( const Field& f, double* array, int iter = 1000000 ) {

	// initialize and fill array of random elements.
	typename Field::RandIter r(f);
	typename Field::Element *elements;
	elements = new typename Field::Element[ iter * 3 ];
	for( int i = 0; i < iter*3; i++ ) {
		do { r.random( elements[i] ); }
			while( f.isZero( elements[i] ) );
	}
	typename Field::Element returnValue;

	UserTimer timer;

	// Compute overhead time for 2-argument functions.
	timer.clear();
	timer.start();
	for( int i = 0; i < iter; i++ ) {
		noop( returnValue, elements[ i*3 ] );
	}
	timer.stop();
	double overHeadTime2 = timer.time();

        // Compute overhead time for 3-argument functions.
        timer.clear();
        timer.start();
        for( int i = 0; i < iter; i++ ) {
                noop( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
        }
        timer.stop();
        double overHeadTime3 = timer.time();

	// Compute overhead time for 4-argument functions.
	timer.clear();
	timer.start();
	for( int i = 0; i < iter; i++ ) {
		noop( returnValue, elements[ i*3 ], elements[ i*3 + 1 ], elements[ i*3 + 2 ] );
	}
	timer.stop();
	double overHeadTime4 = timer.time();


	// add
	timer.clear();
	timer.start();
	for( int i = 0; i < iter; i++ ) {
		f.add( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
	}
	timer.stop();
	array[0] = timer.time() - overHeadTime3;

        // sub
        timer.clear();
        timer.start();
        for( int i = 0; i < iter; i++ ) {
                f.sub( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
        }
        timer.stop();
        array[1] = timer.time() - overHeadTime3;

        // neg
        timer.clear();
        timer.start();
        for( int i = 0; i < iter; i++ ) {
                f.neg( returnValue, elements[ i*3 ] );
        }
        timer.stop();
        array[2] = timer.time() - overHeadTime2;

        // mul
        timer.clear();
        timer.start();
        for( int i = 0; i < iter; i++ ) {
                f.mul( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
        }
        timer.stop();
        array[3] = timer.time() - overHeadTime3;

        // inv
        timer.clear();
        timer.start();
        for( int i = 0; i < iter; i++ ) {
                f.inv( returnValue, elements[ i*3 ] );
        }
        timer.stop();
        array[4] = timer.time() - overHeadTime2;

        // div
        timer.clear();
        timer.start();
        for( int i = 0; i < iter; i++ ) {
                f.div( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
        }
        timer.stop();
        array[5] = timer.time() - overHeadTime3;

        // axpy
        timer.clear();
        timer.start();
        for( int i = 0; i < iter; i++ ) {
                f.axpy( returnValue, elements[ i*3 ], elements[ i*3 + 1 ], elements[ i*3 + 2 ] );
        }
        timer.stop();
        array[6] = timer.time() - overHeadTime4;
}

}

using namespace LinBox;

void printTimings( double* timings ) {
	cout << setw(8) << timings[0]
	     << setw(8) << timings[1]
	     << setw(8) << timings[2]
	     << setw(8) << timings[3]
	     << setw(8) << timings[4]
	     << setw(8) << timings[5]
	     << setw(8) << timings[6];
}

int main(int argc, char** argv) {
	double timings[7];
	long prime = 268435399;
	cout << setw(20) << "Field Name"
	     << setw(8) << "add"
	     << setw(8) << "sub"
	     << setw(8) << "neg"
	     << setw(8) << "mul"
	     << setw(8) << "inv"
	     << setw(8) << "div"
	     << setw(8) << "axpy" << endl;

	DoubleFmod field1( prime );
	fieldTest( field1, timings );
	cout << setw(20) << "DoubleFmod";
	printTimings( timings );
	cout << endl;

	Modular<double> field15( prime );
	fieldTest( field15, timings );
	cout << setw(20) << "Modular<double>";
	printTimings( timings );
	cout << endl;

	NTL_zz_p field2( prime );
	fieldTest( field2, timings );
	cout << setw(20) << "NTL_zz_p";
	printTimings( timings );
	cout << endl;

	Modular<int> field3( prime );
	fieldTest( field3, timings );
	cout << setw(20) << "Modular<int>";
	printTimings( timings );
	cout << endl;

	return 0;
}
