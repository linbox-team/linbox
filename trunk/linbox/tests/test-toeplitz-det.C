#include <iostream>
#include <vector>
#include <linbox/blackbox/toeplitz.h>
#ifdef __LINBOX_HAVE_NTL
#include <linbox/field/ntl-lzz_pX.h>
#include <linbox/field/ntl-lzz_p.h>
#endif
#include <linbox/solutions/det.h>
#include <linbox/blackbox/dense.h>
#include <linbox/randiter/random-prime.h>
#include <cstdlib>
#include <ctime>
#include <linbox/integer.h>

using namespace std;
using namespace LinBox;

const int N_BOUND = 1000;
const int N_REPS = 100;

int main() {
	ostream& report = cout;
	report << "======> Testing Toeplitz Determinant" << endl;
	bool pass = true;
#ifdef __LINBOX_HAVE_NTL
	srandom(time(0));
	RandomPrimeIterator rp;
	NTL_zz_p::RandIter rand;
	report << "\tUsing random primes and square matrices of size 2 to " << N_BOUND << endl;
	int toDel = 0;
	report << '\t';
	for( int i = 0; pass && i < N_REPS; ++i ) {
		if( !(i % (N_REPS/100)) ) {
			int percentDone = (100*i)/N_REPS;
			for( int j = 0; j < toDel; ++j ) report.put(8);
			report << percentDone << "%";
			report.flush();
			toDel = 3;
			if( percentDone < 10 ) --toDel;
		}
		size_t n;
		do { n = random() % N_BOUND; } while( n < 2 );

		NTL_zz_p CF( *rp );
		NTL_zz_pX PF(CF);
		
		DenseMatrix<NTL_zz_p> A(CF,n,n);
		
		NTL_zz_p::Element temp;
		NTL_zz_pX::Element poly;
		PF.init(poly,0);
		size_t r,c;

		for( int diff = 1 - ((int)n); diff <= ((int)n) - 1; ++diff ) {
			rand.random(temp);
			PF.setCoeff(poly,(size_t)(diff + n - 1), temp );
			r = c = 0;
			if( diff < 0 ) c = (size_t)(diff*-1);
			else r = (size_t)diff;
			for( ; r < n && c < n; ++r, ++c )
				A.setEntry(r,c,temp);
		}

		Toeplitz<NTL_zz_p,NTL_zz_pX> T( PF, poly, n );

		NTL_zz_p::Element res1, res2;
		det(res1,T);
		det(res2,T);
		
		if( res1 != res2 ) pass = false;
	}
#endif
	report << endl;
	if( pass ) report << "<====== Passed!" << endl;
	else report << "<====== Failed!" << endl;
	return (pass ? 0 : 1);

#if 0
	NTL_ZZ_pX F( 65521 );
	NTL_ZZ_pX::Element a,b;
	F.init(a,1);
	F.init(b,4);
	F.mulin(a,b);
	F.write(cout,a) << endl;
	F.write(cout) << endl;
	NTL_ZZ_pX::Element T;
	Toeplitz<NTL_ZZ_p> Q(F.getCoeffField());
	vector<int> v;
	v.push_back( -3 );
	v.push_back( 50 );
	v.push_back( 71 );
	v.push_back( -43 );
	v.push_back( 91 );
	v.push_back( 16 );
	v.push_back( 78 );
	F.init(T,v);
	Toeplitz<NTL_ZZ_p,NTL_ZZ_pX> mat( F, T, 4 );
	Toeplitz<NTL_ZZ_p,NTL_ZZ_pX>::Element res;
	//mat.det( res );
	det(res,mat);

	/*
	F.init( T, v );
	NTL_ZZ_pX::Coeff res;
	toeplitz_determinant( F, res, T, 4 );
	*/
	cout << res << endl;
#endif
}
