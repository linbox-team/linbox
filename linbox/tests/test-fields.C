/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-fields.C
 * Written by Dan Roche
 * Copyright (C) June 2004 Dan Roche
 *
 * See COPYING for license information
 */

#include "linbox-config.h"
#include "linbox/util/timer.h"
// #include "linbox/field/givaro-gfq.h"

#ifdef __LINBOX_HAVE_NTL
#include "linbox/field/ntl-lzz_p.h"
#include "linbox/field/ntl-ZZ.h"
#include "linbox/field/ntl-ZZ_p.h"
#include "linbox/field/ntl-pid-lzz_p.h"
#include "linbox/field/PIR-ntl-ZZ_p.h"
#endif

#include "linbox/field/modular.h"
#include "linbox/field/modular-int32.h"
#include "linbox/field/modular-double.h"
#include "linbox/field/field-traits.h"
#include "linbox/vector/stream.h"
#include "linbox/integer.h"
#include "linbox/field/PIR-modular-int32.h"
// #include "linbox/field/gf2.h"
#include "linbox/field/gmp-rational.h"
#include "linbox/field/local2_32.h"
#include "linbox/field/modular-byte.h"
#include "linbox/field/modular-short.h"

#ifdef __LINBOX_HAVE_LIDIA
#include "linbox/field/lidia.h"
#endif

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
		array[i] = iter / (array[i] > 0 ? (array[i] * 1000000) : 0) ;
}

/* This simple test takes and int and a float as arguments, only to make
 * sure the compiler does not optimize too much to make the test useless.
 * The number returned is the number of times per second the inner loop
 * (one floating-point and one int operation) can be executed on the current
 * machine.
 */
int64 getOps(int& a, float& b) {
	int64 ops = 1;
	int64 i = 0;
	a = 13;
	b = 1.3;
	UserTimer opsClock;
	opsClock.clear();
	while( opsClock.time() < 1 ) {
		ops *= 2;
		i = 0;
		opsClock.start();
		while( ++i < ops ) {
			a *= a;
			b *= b;
		}
		opsClock.stop();
	}
	return ops;
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

template <class Field>
void doTest(integer& p, integer& exp, int64& iter) {
	static double mops[10];
	if( FieldTraits<Field>::goodModulus( p ) &&
	    FieldTraits<Field>::goodExponent( exp ) ) {
		Field fld( p, exp );
		fieldTest( fld, mops, iter );
		// print name
		printTimings( mops );
		// cout << endl;
	}
}

int main(int argc, char** argv) {
	int a; float b;
	int64 ops = getOps(a,b);
	cout << "Ops per sec, roughly: " << ops << endl;
	int64 iterations = ops / (1<<8); // should be ops / 32
	integer prime(2), exp(1);
	if( argc >= 2 ) prime = integer( argv[1] );
	if( argc == 3 ) exp = integer( argv[2] );
	if( argc > 3 ) exit(1);

	cout << setw(20) << "Field Name"
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

        cout << setw(20) << "Modular<int32>";
        doTest< Modular<int32> >( prime, exp, iterations );
        cout << endl;

        cout << setw(20) << "Modular<int16>";
        doTest< Modular<int16> >( prime, exp, iterations );
        cout << endl;

        cout << setw(20) << "Modular<int8>";
        doTest< Modular<int8> >( prime, exp, iterations );
        cout << endl;

        cout << setw(20) << "Modular<double>";
        doTest< Modular<double> >( prime, exp, iterations );
        cout << endl;

        cout << setw(20) << "PIRModular<int32>";
        doTest< PIRModular<int32> >( prime, exp, iterations );
        cout << endl;

#ifdef __LINBOX_HAVE_NTL
        cout << setw(20) << "NTL_zz_p";
        doTest< NTL_zz_p >( prime, exp, iterations );
        cout << endl;

        cout << setw(20) << "NTL_PID_zz_p";
        doTest< NTL_PID_zz_p >( prime, exp, iterations );
        cout << endl;

/*
        cout << setw(20) << "NTL_ZZ_p";
        doTest< NTL_ZZ_p >( prime, exp, iterations );
        cout << endl;
*/

        cout << setw(20) << "PIR_ntl_ZZ_p";
        doTest< PIR_ntl_ZZ_p >( prime, exp, iterations );
        cout << endl;

        cout << setw(20) << "NTL_ZZ";
        doTest< NTL_ZZ >( prime, exp, iterations );
        cout << endl;
#endif

#ifdef __LINBOX_HAVE_LIDIA
        cout << setw(20) << "LidiaGfq";
        doTest< LidiaGfq >( prime, exp, iterations );
        cout << endl;
#endif

/*
        cout << setw(20) << "GF2";
        doTest< GF2 >( prime, exp, iterations );
        cout << endl;
*/

        cout << setw(20) << "GMPRationalField";
        doTest< GMPRationalField >( prime, exp, iterations );
        cout << endl;

        cout << setw(20) << "Local2_32";
        doTest< Local2_32 >( prime, exp, iterations );
        cout << endl;

	return 0;
}
