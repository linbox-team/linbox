/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=g0,t0,\:0

/** @file tests/checker.C
@brief script to run LinBox tests

Checker is compiled and run by the check macro invoked by "make check" in the top source dir or in tests/.  It may be run with "make test" or by "make checker; checker" in /tests, avoiding spurrious output of the autotool's check macro.  Run without args, checker prints an explanation of command line options, then proceeds with default test behaviour.

Each LinBox test gets a line below.  A LinBox unit/regresssion test has a default behaviour -- when there are no command line file names -- in which nothing is written to any output stream and the return value is 0 for a passing test, nonzero for any failure.  

The current convention is that (1) linbox' check macro, checker.C, runs the tests with no command line parameters at all, and (2) if there is a command line file name, verbose diagnostic output is written to that file and more terse output may be written to standard output streams.  The second feature is intended to assist debugging with individual tests.  

In future, command line arguments could be used in checking to vary matrix sizes and other parameters in check runs.

There should be a 1-1 correspondence between files tests/test-*.C and calls
here to build_n_run() or no_build_n_run() (possibly commented out).
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
//#include <iomanip>
using namespace std;

#include "../linbox/linbox-config.h"

struct counts {
	int pass; int buildfail; int runfail; int skipped;
	counts() : pass(0), buildfail(0), runfail(0), skipped(0) {}
	int total() { return pass + buildfail + runfail; }
};

void build_n_run(string s, counts& cnt, int flag = 1) ;
void no_build_n_run(string s, string r, counts& cnt, int flag) ;

int main(int argc, char* argv[]){
	int flag = 1; // default verbosity, no force
	if (argc < 2) {
		cout << "usage: " << argv[0] << " [summary,default,errors,verbose] [c]" << endl;
		cout << "  summary: only 4 lines printed." << endl;
		cout << "  default: also one line per test." << endl;
		cout << "  errors: also any build and run output for error cases." << endl;
		cout << "  verbose: also any build and run output for each test." << endl;
		cout << "  2nd arg, if present, forces rebuild of all tests." << endl;
	} else { 
		if (argv[1][0] == 's') flag = 0; // summary
		if (argv[1][0] == 'd') flag = 1; // default
		if (argv[1][0] == 'e') flag = 2; // errors
		if (argv[1][0] == 'v') flag = 3; // verbose
		if (argc > 2) flag += 4; // force rebuilds
	}

#if 0
	// the setup
	system("rm -f checkdata");
	string compiler;
	//system("make -n "); 
	system("make -n checkerdata > checkdata");
	ifstream fin("checkdata", ifstream::in);
	do {getline(fin, compiler);}
	while (fin && compiler.find("something") == compiler.npos);
	cout << compiler << endl << endl;
	system("rm -f checkdata");
#endif

	// the tests
	counts counter;

/*
Each test gets a line below.  A LinBox unit/regresssion test has a default behaviour -- when there are no command line file names -- in which nothing is written to any output stream and the return value is 0 for a passing test, nonzero for any failure.  

The current convention is that (1) linbox' check macro, checker.C, runs the tests with no command line parameters at all, and (2) if there is a command line file name, verbose diagnostic output is written to that file and more terse output may be written to standard output streams.  The second feature is intended to assist debugging with individual tests.  

In future, command line arguments could be used in checking to vary matrix sizes and other parameters in check runs.

There should be a 1-1 correspondence between files tests/test-*.C and calls
here to build_n_run or no_build_n_run (possibly commented out).
Thus "ls test-*.C |wc" and "grep test- checker.C |grep \
build |wc" should yield the same number of lines.
*/

//BASIC_TESTS
	build_n_run("test-bitonic-sort", counter, flag);
	build_n_run("test-blas-domain", counter, flag);
	build_n_run("test-butterfly", counter, flag);
	build_n_run("test-commentator", counter, flag);
	build_n_run("test-companion", counter, flag);
	build_n_run("test-dense", counter, flag);
	build_n_run("test-det", counter, flag);
	build_n_run("test-diagonal", counter, flag);
	build_n_run("test-dif", counter, flag);
	build_n_run("test-direct-sum", counter, flag);
	build_n_run("test-ffpack", counter, flag);
	build_n_run("test-frobenius", counter, flag);
	build_n_run("test-getentry", counter, flag);
	build_n_run("test-gf2", counter, flag);
	build_n_run("test-gmp-rational", counter, flag);
	build_n_run("test-hilbert", counter, flag);
	build_n_run("test-hom", counter, flag);
	build_n_run("test-inverse", counter, flag);
	build_n_run("test-isposdef", counter, flag);
	build_n_run("test-ispossemidef", counter, flag);
	build_n_run("test-matrix-domain", counter, flag);
	build_n_run("test-matrix-stream", counter, flag);
	build_n_run("test-minpoly", counter, flag);
	build_n_run("test-modular", counter, flag);
	build_n_run("test-modular-balanced-int", counter, flag);
	build_n_run("test-modular-byte", counter, flag);
	build_n_run("test-modular-double", counter, flag);
	build_n_run("test-modular-int", counter, flag);
	build_n_run("test-modular-short", counter, flag);
	build_n_run("test-block-ring", counter, flag);
	build_n_run("test-moore-penrose", counter, flag);
	build_n_run("test-randiter-nonzero", counter, flag);
	build_n_run("test-rank", counter, flag);
	build_n_run("test-rational-matrix-factory ", counter, flag);
	build_n_run("test-scalar-matrix", counter, flag);
	build_n_run("test-sparse", counter, flag);
	build_n_run("test-subiterator", counter, flag);
	build_n_run("test-submatrix", counter, flag);
	build_n_run("test-subvector", counter, flag);
	build_n_run("test-trace", counter, flag);
	build_n_run("test-triplesbb", counter, flag);
	build_n_run("test-vector-domain", counter, flag);
	build_n_run("test-qlup", counter, flag);
	build_n_run("test-charpoly", //"intermittent inf loop, bb or cp responsible?", 
	counter, flag);
	build_n_run("test-zero-one", counter, flag);
	build_n_run("test-solve", counter, flag);
	build_n_run("test-nullspace", counter, flag);
	build_n_run("test-smith-form-binary", counter, flag);
	build_n_run("test-rational-solver", counter, flag);
	build_n_run("test-dyadic-to-rational", counter, flag);
	build_n_run("test-last-invariant-factor", counter, flag);
	build_n_run("test-cra", counter, flag);
	build_n_run("test-cradomain", counter, flag);
	build_n_run("test-rational-reconstruction-base", counter, flag);
#if __LINBOX_HAVE_LAPACK
	if (flag > 0) cout << "	Lapack dependent tests" << endl;
	build_n_run("test-rational-solver-adaptive", counter, flag);
	// needs output cleanup.  Resolve whether a benchmark or a test.
	build_n_run("test-solve-nonsingular", counter, flag);
#endif

#if __LINBOX_HAVE_GIVARO
	if (flag > 0) cout << "	Givaro tests" << endl;
	build_n_run("test-givaro-zpz", counter, flag);
	build_n_run("test-givaro-zpzuns", counter, flag);
	build_n_run("test-rat-solve", counter, flag);
	build_n_run("test-rat-minpoly", counter, flag);
#else 
	if (flag > 0) cout << "	not doing Givaro dependent tests" << endl;
	cout << "Configuration problem?  __LINBOX_HAVE_GIVARO is not set, but LinBox requires Givaro" << endl;
#endif

#if __LINBOX_HAVE_NTL
	if (flag > 0) cout << "	NTL dependent tests" << endl;
	build_n_run("test-ntl-lzz_p", counter, flag);
	build_n_run("test-ntl-ZZ_p", counter, flag);
	build_n_run("test-ntl-toeplitz", counter, flag);
	build_n_run("test-ntl-hankel", counter, flag);
	build_n_run("test-ntl-sylvester", counter, flag);
	build_n_run("test-sum", counter, flag);
	build_n_run("test-smith-form-iliopoulos", counter, flag);
	build_n_run("test-toeplitz-det", counter, flag);
#else 
	if (flag > 0) cout << "	not doing NTL dependent tests" << endl;
#endif

#if __LINBOX_HAVE_LIDIA
	if (flag > 0) cout << "	Lidia dependent tests" << endl;
	build_n_run("test-lidia-gfq", counter, flag);
#else 
	if (flag > 0) cout << "	not doing Lidia dependent test" << endl;
#endif

#if __LINBOX_HAVE_ATLAS
	if (flag > 0) cout << "	Atlas dependent tests" << endl;
	build_n_run("test-optimization", counter, flag);
#else 
	if (flag > 0) cout << "	not doing Atlas dependent test" << endl;
#endif

#if 1
// tests of components not adhering to LinBox requirements.
//NONCOMPLIANT_TESTS = 
	if (flag > 0) cout << "	Noncompliant tests" << endl;
	build_n_run("test-unparametric-field", counter, flag);
	no_build_n_run("test-ntl-RR", "make fails, nobody responsible", counter, flag);
	build_n_run("test-param-fuzzy", counter, flag);
//NONCOMPILING_TESTS =
	if (flag > 0) cout << "	Noncompiling tests" << endl;
	no_build_n_run("test-dense-zero-one", "make fails, bds responsible", counter, flag);
	no_build_n_run("test-smith-form", "make fails, bds responsible", counter, flag);
	no_build_n_run("test-smith-form-local", "make fails, bds responsible", counter, flag);
	no_build_n_run("test-smith-form-adaptive", "make fails, bds responsible", counter, flag);
//FAILING_TESTS = 
	if (flag > 0) cout << "	Failing tests" << endl;
	no_build_n_run("test-rat-charpoly", "infinite loop, cp responsible?", counter, flag);
//DUBIOUS_TESTS =     
	if (flag > 0) cout << "	Immature tests" << endl;
	no_build_n_run("test-ftrmm", "half baked, Boyer's?",  counter, flag);
	// test-integer-tools -- there is no test-integer-tools.C file
	no_build_n_run("test-quad-matrix", "half baked, make fails, bds responsible", counter, flag);
	no_build_n_run("test-zo", "half baked, make fails, bds responsible", counter, flag);
	// no one has taken these on.
	//no_build_n_run("test-mg-block-lanczos", "make fails, nobody responsible", counter, flag);
	//no_build_n_run("test-la-block-lanczos", "make fails, nobody responsible", counter, flag);
#endif

#if 0
	if (flag > 0) cout << "	Non tests" << endl;
	no_build_n_run("test-template", "yup, model for tests, not an actual test", counter, flag);
	// test-common.C is included in test-common.h, 
	//  but could be separately compiled (no templates).
	no_build_n_run("test-common", "not a test", counter, flag);
	no_build_n_run("test-fields", "deprecated", counter, flag);
	no_build_n_run("test-matrix", "deprecated?", counter, flag);
	no_build_n_run("test-image-field", "deprecated", counter, flag);

// a benchmarker, not to be included in check.
//BENCHMARKS =
	build_n_run("benchmark-fields", counter, flag);
	build_n_run("benchmark-blas-domain", counter, flag);

# I put test-gssv_rank in hmrg:~saunders/gen_superlu.  It is not linbox until and unless it is made to work
# test-gssv is an ntl test

# no explicit test needed, I guess...
#        Transpose is tested in test-triplesbb
#        Compose is tested in test-butterfly

#endif

	// the summary
	cout << "--------------------------------------------------------------------" << endl;
	if (counter.buildfail || counter.runfail)
		cout << "Of " << counter.total() << " tests, " 
			<< counter.buildfail << " failed to compile and " 
			<< counter.runfail << " failed at runtime." 
			<< "  (" << counter.skipped << " skipped tests.)" << endl;
	else
		cout << endl << "All " << counter.total() << " tests pass." 
			<< "  (" << counter.skipped << " skipped tests.)" << endl;
	cout << "  File tests/checkdata contains any output from the builds and runs." << endl;
	cout << "--------------------------------------------------------------------" << endl;

	return counter.buildfail || counter.runfail ? -1 : 0;
}

void no_build_n_run(string s, string r, counts& cnt, int flag) {
	if (flag >= 4) flag -= 4;
	if (flag > 0) {
		cout.width(35); cout << left << s;
		cout << "not tried, " << r << "." << endl;
	}
	cnt.skipped++;
}

void build_n_run(string s, counts& cnt, int flag) {
	string cmd;
	if (flag >= 4) { // force build
		flag -= 4;
		cmd = "touch " + s + ".C";
		system(cmd.c_str());
	}
	if (flag > 0) {
		cout.width(35); cout << left << s;
		cout << "make"; cout.flush();
	}
	cmd = "make " + s + " 2>> checkdata >> checkdata"; 
	int status = system(cmd.c_str());
	if (status != 0) {
		if (flag > 0) cout << " FAILS" << endl;
		if (flag > 1) system ("cat checkdata");
		cnt.buildfail++;
		//cerr << "compile fails " << endl;
		//system("cat checkdata");
	} else {
		if (flag > 0) { cout << "\b\b\b\brun "; cout.flush(); }
		status = system(s.c_str());
		if (status != 0) { 
			if (flag > 0) cout << "\b\b\b\bFAILS" << endl;
			if (flag > 1) system ("cat checkdata");
			cnt.runfail++;
			//cerr << "run fails " << endl;
		} else {
			if (flag > 0) cout << "\b\b\b\bOK  " << endl;
			if (flag > 2) system ("cat checkdata");
			cnt.pass++;
			//cerr << "ok" << endl;
		}
	}
	if (flag > 1) system("rm -f checkdata");
}

