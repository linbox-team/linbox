/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-fields.C
 * Written by Dan Roche
 * Copyright (C) June 2004 Dan Roche
 *
 * See COPYING for license information
 */

#include "linbox/util/timer.h"
#include "linbox/field/givaro-gfq.h"
#include "linbox/field/ntl-zz_p.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-int.h"
#include "linbox/field/modular-double.h"
#include "linbox/vector/stream.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setw;

// Namespace in which all LinBox code resides
namespace LinBox {

template< class Element >
Element& noop( Element& a, const Element& b ) { return a = b; }

/* fieldTest is a template function to test out the performance of a given field on a
 * machine.  Taken are three arguments.  The first is a field class object.  The second
 * is an array, declared but not necessarily initialized, of ten doubles.  The
 * first nine values will be filled with mops for add, sub, neg, mul, int, div,
 * axpy, dot1, and dot2, respectively.  (Dot1 is dense*dense, Dot2 is dense*sparse).
 * The last value is filled with mops for walking through an array of size iter.
 * The third argument is optional and specifies how many loop iterations to use.
 */

template< class Field >
void fieldTest( const Field& f, double* array, long iter = 1000000 ) {

	long vectorSize = 10000;
	float sparsity = .01;
	int i;

	// initialize and fill array of random elements.
	typename Field::RandIter r(f);
	typename Field::Element *elements;
	elements = new typename Field::Element[ iter * 3 ];
	for( int i = 0; i < iter*3; i++ ) {
		do { r.random( elements[i] ); }
			while( f.isZero( elements[i] ) );
	}
	register typename Field::Element returnValue;

	// initialize random vector streams
	RandomDenseStream<Field> dense1( f, vectorSize, iter/vectorSize );
	RandomDenseStream<Field> dense2( f, vectorSize, iter/vectorSize );
	RandomSparseStream<Field> sparse( f, sparsity, vectorSize, 0 );

	// initialize individual vectors to hold results
	typename RandomDenseStream<Field>::Vector dv1;
	typename RandomDenseStream<Field>::Vector dv2;
	typename RandomSparseStream<Field>::Vector sv;
	VectorWrapper::ensureDim (dv1,vectorSize);
	VectorWrapper::ensureDim (dv2,vectorSize);
	VectorWrapper::ensureDim (sv,vectorSize);

	VectorDomain<Field> VD( f );

	UserTimer timer;

	// Compute overhead time for non-vector functions.
	timer.clear();
	timer.start();
	for( i = 0; i < iter; i++ ) {
		noop( returnValue, elements[ i*3 ] );
		
	}
	timer.stop();
	double overHeadTime = timer.time();

	// Compute overhead time for DotProduct1( dense*dense )
	timer.clear();
	timer.start();
	while( dense1 && dense2 ) {
		dense1.get( dv1 );
		dense2.get( dv2 );
	}
	timer.stop();
	dense1.reset();
	dense2.reset();
	double overHeadTimeDD = timer.time();

	// Compute overhead time for DotProduct2( dense*sparse )
	timer.clear();
	timer.start();
        while( dense1 ) {
                dense1.get( dv1 );
                for( float fi = 0; fi < 1; fi += sparsity ) sparse.get( sv );
        }
	timer.stop();
	dense1.reset();
	sparse.reset();
	double overHeadTimeDS = timer.time();

	// add
	timer.clear();
	timer.start();
	for( i = 0; i < iter; i++ ) {
		f.add( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
	}
	timer.stop();
	array[0] = timer.time() - overHeadTime;

        // sub
        timer.clear();
        timer.start();
        for( i = 0; i < iter; i++ ) {
                f.sub( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
        }
        timer.stop();
        array[1] = timer.time() - overHeadTime;

        // neg
        timer.clear();
        timer.start();
        for( i = 0; i < iter; i++ ) {
                f.neg( returnValue, elements[ i*3 ] );
        }
        timer.stop();
        array[2] = timer.time() - overHeadTime;

        // mul
        timer.clear();
        timer.start();
        for( i = 0; i < iter; i++ ) {
                f.mul( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
        }
        timer.stop();
        array[3] = timer.time() - overHeadTime;

        // inv
        timer.clear();
        timer.start();
        for( i = 0; i < iter; i++ ) {
                f.inv( returnValue, elements[ i*3 ] );
        }
        timer.stop();
        array[4] = timer.time() - overHeadTime;

        // div
        timer.clear();
        timer.start();
        for( i = 0; i < iter; i++ ) {
                f.div( returnValue, elements[ i*3 ], elements[ i*3 + 1 ] );
        }
        timer.stop();
        array[5] = timer.time() - overHeadTime;

        // axpy
        timer.clear();
        timer.start();
        for( i = 0; i < iter; i++ ) {
                f.axpy( returnValue, elements[ i*3 ], elements[ i*3 + 1 ],
			elements[ i*3 + 2 ] );
        }
        timer.stop();
        array[6] = timer.time() - overHeadTime;

	// DotProduct1 ( dense * dense )
	timer.clear();
	timer.start();
        while( dense1 && dense2 ) {
                dense1.get( dv1 );
                dense2.get( dv2 );
		VD.dot( returnValue, dv1, dv2 );
        }
	timer.stop();
	dense1.reset();
	dense2.reset();
	array[7] = timer.time() - overHeadTimeDD;

	// DotProduct2 ( dense * sparse )
	timer.clear();
	timer.start();
        while( dense1 && sparse ) {
                dense1.get( dv1 );
                for (float fi = 0; fi < 1; fi += sparsity) {
			sparse.get( sv );
			VD.dot( returnValue, dv1, sv );
		}
	}
	timer.stop();
	dense1.reset();
	sparse.reset();
	array[8] = timer.time() - overHeadTimeDS;

	// Convert timings to mops (million operations per second)
	array[9] = overHeadTime;
	for( i = 0; i < 10; i++ )
		array[i] = iter / (array[i] * 1000000);
}

}

using namespace LinBox;

void printTimings( double* timings ) {
	cout << setw(8) << timings[0] << ' '
	     << setw(8) << timings[1] << ' '
	     << setw(8) << timings[2] << ' '
	     << setw(8) << timings[3] << ' '
	     << setw(8) << timings[4] << ' '
	     << setw(8) << timings[5] << ' '
	     << setw(8) << timings[6] << ' '
	     << setw(8) << timings[7] << ' '
	     << setw(8) << timings[8] << ' '
	     << setw(8) << timings[9];
}

int main(int argc, char** argv) {
	double timings[10];
	long prime = 65521;
	cout << setw(15) << "Field Name"
	     << setw(8) << "add"
	     << setw(9) << "sub"
	     << setw(9) << "neg"
	     << setw(9) << "mul"
	     << setw(9) << "inv"
	     << setw(9) << "div"
	     << setw(9) << "axpy"
	     << setw(9) << "dot d*d"
	     << setw(9) << "dot d*s"
	     << setw(9) << "array" << endl;

	Modular<double> field1( prime );
	fieldTest( field1, timings );
	cout << setw(15) << "Modular<double>";
	printTimings( timings );
	cout << endl;

	NTL_zz_p field2( prime );
	fieldTest( field2, timings );
	cout << setw(15) << "NTL_zz_p";
	printTimings( timings );
	cout << endl;

	Modular<int> field3( prime );
	fieldTest( field3, timings );
	cout << setw(15) << "Modular<int>";
	printTimings( timings );
	cout << endl;

	return 0;
}
