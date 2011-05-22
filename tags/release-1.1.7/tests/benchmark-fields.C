
/* tests/test-fields.C
 * Written by Dan Roche
 * Copyright (C) June 2004 Dan Roche, part of LinBox, GNU LGPL. See COPYING for license.
 */

#include "linbox/linbox-config.h"
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
#include "linbox/field/modular-int.h"
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
#include <vector>

using namespace LinBox;

/* fieldTest is a template function to test out the performance of a given field on a
 * machine.  Taken are three arguments.  The first is a field class object.  The second
 * is an array, declared but not necessarily initialized, of ten doubles.  The
 * first nine values will be filled with mops for add, sub, neg, mul, int, div,
 * axpy, dot1, and dot2, respectively.  (Dot1 is dense*dense, Dot2 is dense*sparse).
 * The last value is filled with mops for walking through an array of size iter.
 * The third argument is optional and specifies how many loop iterations to use.
 */

template< class Field >
void fieldTest( const Field& f, double* array, long iter = 1000000, bool fulltest = false ) {

	long vectorSize = 10000;
	float sparsity = .01;
	long sparsity_inv = 100;
	int i;

	// initialize a few field elements,
	typedef typename Field::Element Element;
	register Element returnValue; f.init(returnValue, 1);
	register Element s; f.init(s, 0); 

	register Element a, b, c;
	typename Field::RandIter r(f);
	r.random( a ); r.random( b ); r.random( c ); 
	std::vector<Element> dv1( vectorSize ), dv2( vectorSize );
	for (i = 0; i < vectorSize; ++i ) {
		r.random( dv1[i] );
		r.random( dv2[i] );
	}
	RandomSparseStream<Field> sparse( f, sparsity, vectorSize ); 
	typename RandomSparseStream<Field>::Vector sv; sparse.get( sv );

/*
	// initialize and fill array of random elements.
	typename Field::RandIter r(f);
	typename Field::Element *elements;
	elements = new typename Field::Element[ iter * 3 ];
	for( int i = 0; i < iter*3; i++ ) {
		do { r.random( elements[i] ); }
			while( f.isZero( elements[i] ) );
	}

	// initialize random vector streams
	RandomDenseStream<Field> dense( f, vectorSize, 2);
	typename RandomDenseStream<Field>::Vector dv1; dense.get( dv1 );
	typename RandomDenseStream<Field>::Vector dv2; dense.get( dv2 );
	RandomSparseStream<Field> sparse( f, sparsity, vectorSize ); 
	typename RandomSparseStream<Field>::Vector sv; sparse.get( sv );

	RandomDenseStream<Field> dense1( f, vectorSize, iter/vectorSize );
	RandomDenseStream<Field> dense2( f, vectorSize, iter/vectorSize );
	RandomSparseStream<Field> sparse( f, sparsity, vectorSize ); 

	// initialize individual vectors to hold results
	typename RandomDenseStream<Field>::Vector dv1;
	typename RandomDenseStream<Field>::Vector dv2;
	typename RandomSparseStream<Field>::Vector sv;
*/
	VectorWrapper::ensureDim (dv1,vectorSize);
	VectorWrapper::ensureDim (dv2,vectorSize);
	VectorWrapper::ensureDim (sv,vectorSize);

	VectorDomain<Field> VD( f );

	UserTimer timer;
	double overHeadTime;

	timer.clear(); timer.start();
	f.init(s, 0); 
	for( i = 0; i < iter; i++ ) { f.init(returnValue, i); f.addin(s, returnValue); }
	timer.stop(); overHeadTime = timer.time();

	// add
	timer.clear(); timer.start();
	for( i = 0; i < iter; i++ ) {
		f.init(a, i);
		f.add( returnValue, a, b);
		f.addin(s, returnValue);
	}
	timer.stop(); array[0] = timer.time() - overHeadTime;
//std::cout << iter << " add done " << array[0] << std::endl;

if (fulltest) {
    // sub
    timer.clear(); timer.start();
    for( i = 0; i < iter; i++ ) {
		f.init(a, i);
		f.sub( returnValue, a, b);
		f.addin(s, returnValue);
    }
    timer.stop(); array[1] = timer.time() - overHeadTime;

	// neg
	timer.clear(); timer.start();
	for( i = 0; i < iter; i++ ) {
		f.init(a, i);
		f.neg( returnValue, a);
		f.addin(s, returnValue);
	}
	timer.stop(); array[2] = timer.time() - overHeadTime;
} // end if (fulltest)

	// mul
	timer.clear(); timer.start();
	for( i = 0; i < iter; i++ ) {
		f.init(a, i);
		f.mul( returnValue, a, b);
		f.addin(s, returnValue);
	}
	timer.stop(); array[3] = timer.time() - overHeadTime;
//std::cout << iter << " mul done " << array[3] << std::endl;

if (fulltest) {
	// inv
	timer.clear(); timer.start();
	for( i = 0; i < iter; i++ ) {
		f.init(a, i);  if (f.isZero(a)) f.init(a, 1);
		f.inv( returnValue, a);
		f.addin(s, returnValue);
	}
	timer.stop(); array[4] = timer.time() - overHeadTime;

	// div
	timer.clear(); timer.start();
	for( i = 0; i < iter; i++ ) {
		f.init(a, i); 
		f.div( returnValue, a, b);
		f.addin(s, returnValue);
	}
	timer.stop(); array[5] = timer.time() - overHeadTime;
} // end if (fulltest)

	// axpy
	timer.clear(); timer.start();
	for( i = 0; i < iter; i++ ) {
		f.init(a, i);
		f.axpy( returnValue, a, b, c);
		f.addin(s, returnValue);
	}
	timer.stop(); array[6] = timer.time() - overHeadTime;
	//std::cout << timer.time() << "  - " << overHeadTime << " = " << array[6] << std::endl;;

	// DotProduct1 ( dense * dense )
	timer.clear(); timer.start();
	for( i = 0; i < iter/vectorSize; i++ ) {
		f.init(dv1.back(), i);
		VD.dot( returnValue, dv1, dv2 );
		f.addin(s, returnValue);
	}
	timer.stop(); array[7] = timer.time();
//	std::cout << (iter/vectorSize) << " dd " << timer.time() << std::endl;;

if (fulltest) {
	// DotProduct2 ( dense * sparse )
	timer.clear(); timer.start();
	for( i = 0; i < iter/vectorSize; i++ ) {
		f.init(dv1.back(), i);
    	for ( int j = 0; j < sparsity_inv; ++ j ) {
			f.init(dv1.front(), j);
			VD.dot( returnValue, dv1, sv );
			f.addin(s, returnValue);
		}
	}
	timer.stop(); array[8] = timer.time();
	//std::cout << "ds " << timer.time() << std::endl;;
} // end if (fulltest)

	// Convert timings to mops (million operations per second)
	for( i = 0; i < 9; i++ ) {	
		double t = array[i];
		array[i] = iter / (t > 0 ? (t * 1000000) : 0) ;
	}
	// use s (just in case compiler cares)
	if (f.isZero(s)) std::cout << "zero sum" << std::endl;
}

/* This simple test takes and int and a float as arguments, only to make
 * sure the compiler does not optimize too much to make the test useless.
 * The number returned is the number of times per second the inner loop
 * (one floating-point and one int operation) can be executed on the current
 * machine.
 */
int64 getOps(int unit) {
	int64 ops = 1;
	int64 i = 0;
	int a = 13;
	double b = 1.3;
	UserTimer opsClock;
	opsClock.clear();
	long double c;
	while( opsClock.time() < unit ) {
		ops *= 2;
		i = 0;
		opsClock.start();
		while( ++i < ops ) {
			a *= a;
			b *= b;
		}
		opsClock.stop();
		// random code to prevent optimization of the loop
		if (a<b)
			b=a;
		else
			b = 2*a;
		c = a+b;

	}
	return ops;
}

void printTimings( double* timings, bool fulltest = false ) {
	if (fulltest){ std::cout 
	     << std::setw(11) << timings[0] << ' '
	     << std::setw(11) << timings[1] << ' '
	     << std::setw(11) << timings[2] << ' '
	     << std::setw(11) << timings[3] << ' '
	     << std::setw(11) << timings[4] << ' '
	     << std::setw(11) << timings[5] << ' '
	;} std::cout 
	     << std::setw(11) << timings[6] << ' '
	     << std::setw(11) << timings[7] << ' '
	; if (fulltest){ std::cout 
	     << std::setw(11) << timings[8] << ' '
	;} std::cout 
	     << std::setw(11) << timings[6]/(1/(1/timings[0] + 1/timings[3])); // axpy/(mul+add) ratio
}

template <class Field>
void doTest(const char* name, integer& p, integer& exp, int64& iter, bool fulltest = false) {
	static double mops[11];
	if( FieldTraits<Field>::goodModulus( p ) &&
	    FieldTraits<Field>::goodExponent( exp ) ) {
		Field fld( p, exp );
		fieldTest( fld, mops, iter, fulltest);
		// print name
		std::cout << std::setw(20) << name;
		printTimings( mops, fulltest);
		std::cout << std::endl;
	} 
	else {
		std::cout << std::setw(20) << name << ": " << p << "^" << exp << " is out of range" << std::endl;
	}
}

int main(int argc, char** argv) {
	int64 ops = getOps(1);
	std::cout << "timings recorded in mops.  Bigger is better." << std::endl;
	std::cout << "Ops per sec, roughly: " << ops << std::endl;
	//int64 iterations = ops/16;
	int64 iterations = ops;
	integer prime(101), exp(1);
	if( argc >= 2 ) prime = integer( argv[1] );
	if( argc >= 3 ) exp = integer( argv[2] );
	//bool fulltest = true;
	bool fulltest = false;
	if( argc > 3 ) fulltest = ( argv[3][0] == 1 ); 
	if( argc > 4 ) exit(1);

	std::cout << std::setw(20) << "Field Name";
	if (fulltest) { std::cout 
	     << std::setw(12) << "add "
	     << std::setw(12) << "sub "
	     << std::setw(12) << "neg "
	     << std::setw(12) << "mul "
	     << std::setw(12) << "inv "
	     << std::setw(12) << "div "
	;} std::cout 
	     << std::setw(12) << "axpy"
	     << std::setw(12) << "dot d*d "
	; if (fulltest) { std::cout 
	     << std::setw(12) << "dot d*s "
	;} std::cout 
	     << std::setw(12) << "axpy/(mul+add)"
		 << std::endl;

    doTest< Modular<int8> >( "Modular<int8>", prime, exp, iterations, fulltest );
    doTest< Modular<int16> >( "Modular<int16>", prime, exp, iterations, fulltest );
    doTest< Modular<int32> >( "Modular<int32>", prime, exp, iterations, fulltest );
    doTest< Modular<int> >( "Modular<int>", prime, exp, iterations, fulltest );
    doTest< Modular<double> >( "Modular<double>", prime, exp, iterations, fulltest );
    doTest< Modular<float> >( "Modular<float>", prime, exp, iterations, fulltest );

#ifdef __LINBOX_HAVE_NTL
    doTest< NTL_zz_p >( "NTL_zz_p", prime, exp, iterations, fulltest );
    doTest< NTL_PID_zz_p >( "NTL_PID_zz_p", prime, exp, iterations, fulltest ); 
    doTest< NTL_ZZ_p >( "NTL_ZZ_p", prime, exp, iterations, fulltest );
    doTest< PIR_ntl_ZZ_p >( "PIR_ntl_ZZ_p", prime, exp, iterations, fulltest );
    doTest< NTL_ZZ >( "NTL_ZZ", prime, exp, iterations, fulltest );
#endif
#ifdef __LINBOX_HAVE_LIDIA
    doTest< LidiaGfq >( "LidiaGfq", prime, exp, iterations, fulltest );
#endif
//	doTest< GF2 >( "GF2", prime, exp, iterations, fulltest );
    doTest< GMPRationalField >( "GMPRationalField", prime, exp, iterations, fulltest ); 
	//if (prime == 2)
    	doTest< PIRModular<int32> >( "PIRModular<int32>", prime, exp, iterations, fulltest );
    doTest< Local2_32 >( "Local2_32", prime, exp, iterations, fulltest );

	return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
