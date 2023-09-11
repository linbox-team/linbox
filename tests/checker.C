/*
 * Copyright (c) 2015 LinBox
 * This file is part of the library LinBox. See COPYING for license info.
 * Written by bds.
 */
/** @file tests/checker.C
    @brief script to run LinBox tests

    Checker is compiled and run by the check macro invoked by "make fullcheck" in the top source dir or in tests/.

    Run without args,
    each test is build and run, with a one line report of success or failure.
    A summary at the end reports number of test failing to compile, failing at runtime, and skipped.
    There should be a 1-1 correspondence between files tests/test-*.C and report lines

    A LinBox unit/regresssion test has a default behaviour -- when there are no command line file names -- in which nothing is written to any output stream and the return value is 0 for a passing test, nonzero for any failure.

    The current convention is that (1) linbox' checker.C, runs the tests with no command line parameters at all, and (2) if there is a command line file name, verbose diagnostic output is written to that file and more terse output may be written to standard output streams.  The second feature is intended to assist debugging with individual tests.

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
//#include <iomanip>
using namespace std;
//#include "linbox-config.h"
#include "linbox/linbox-config.h"
//#include "fflas-ffpack/fflas-ffpack.h"

// globals
//map< string, string> skip_note;
map< string, string> warn_note;
set< string> skips;
set< string> hides;
void warn(const string& t, const string& w) { warn_note[t] = w; }
void skip(const string& t, const string& w) { warn_note[t] = w; skips.insert(t); }
void hide(const string& t, const string& w) { warn_note[t] = w; hides.insert(t); }

int main(int argc, char* argv[]) {
	int arg;
	bool force_build = false, reporting = true, honor_skips = true;
	for (arg = 1; arg < argc; ++arg)
		for (char * i = argv[arg]; *i != 0; ++i){
			if (*i == 'r') force_build = true;
			if (*i == 's') reporting = false;
			if (*i == 'a') honor_skips = false;
		}
    if (argc > 1 and not force_build and reporting and honor_skips){ // bogus arg
		cout << "usage: " << argv[0] << " [-][r][s][a]" << endl;
		cout << "  -r Force recompilation of each test."<< endl;
		cout << "  -s Summary only: 3 lines printed. (Default is one line per test.)" << endl;
		cout << "  -a Build and run all. (Ignore skip and hide commands.)" << endl;
		return 0;
	}

// fullcheck customization: warnings, skips, hides, optional package dependencies

//// notes section (customize fullcheck here) ////

    if (honor_skips) {
        hide("test-block-ring", "non commutative rings not supported");
//hide("test-ffpack", "testTURBO fails, move to ffpack tests?");
        hide("test-ftrmm", "should move to attic");
//         skip("test-givaro-fields", "may fail on small fields because of supposed non-randomness or failure to find a non trivial element");
        hide("test-image-field", "deprecated");
// skip("test-invariant-factors", "not unit/regression test conforming");
//skip("test-isposdef", "intermittent inf loop");
//skip("test-ispossemidef", "intermittent inf loop");
        hide("test-la-block-lanczos", "not maintained. operator >> missing");
        hide("test-mg-block-lanczos", "not maintained");
        hide("test-modular", "deprecated");
        hide("test-modular-byte",  "deprecated");
        hide("test-modular-short",  "deprecated");
//skip("test-modular-balanced-int",  "test and modular-balanced disagree on init/convert");
//skip("test-modular-balanced-double",  "test and modular-balanced disagree on init/convert");
//skip("test-moore-penrose", "inf loop");
        skip("test-optimization", "not unit/regression test conforming");
//skip("test-rational-reconstruction-base", "inf loop");
        skip("test-rat-charpoly", "inf loop");
        skip("test-rat-minpoly", "stale test. solns over QQ need fresh tests"); // "intermittent failures")
        skip("test-rat-solve", "stale test. solns over QQ need fresh tests"); // "infinite loop")
        skip("test-poly-det", "incomplete test (if still relevant)");
        skip("test-sparse-map-map", "const issue in givranditer, curious use of nonexistant next() in Extension");
//Tests requiring further development
    }

//warn("test-echelon-form", "new");
    warn("test-fibb",  "incomplete");
//warn("test-matrix-domain", "intermittent row permutation failure");
    warn("test-param-fuzzy", "Noncompliant field");
//		template <class F2Field> BitVector (const F2Field&) {}
//warn("test-qlup", "GF2 fails to compile");
//warn("test-rank-u32", "intermittent failure"/*, "vector (bb) responsible"*/);
    warn("test-rat-charpoly", "stale test. solns over QQ need fresh tests");//, "infinite loop, cp responsible?")
    warn("test-rat-solve", "infinite loop");
//warn("test-smith-form-local", "bds, intermittent failures");
    warn("test-solve", "most of the tests are commented out");
    warn("test-toom-cook", " One method does not work");
//warn("test-transpose", "sometimes fails on Sparsematrix/getEntry");
/* Quad matrix is a dormant project. Eventual revival is expected. Test works.
    warn("test-quad-matrix", "half baked, bds responsible");
*/

    warn("test-one-invariant-factor", " Probabilistic algorithm, sometimes fails");


//// optional package dependency section ////
	set< string> ntl_tests;
	ntl_tests.insert("test-ntl-hankel");
	ntl_tests.insert("test-ntl-lzz_p");
	ntl_tests.insert("test-ntl-lzz_pe");
	ntl_tests.insert("test-ntl-lzz_px");
	ntl_tests.insert("test-ntl-lzz_pex");
	warn("test-toeplitz", "we should have a non NTL version.");
	ntl_tests.insert("test-toeplitz");
	warn("test-toeplitz-det", "we should have a non NTL version.");
	ntl_tests.insert("test-toeplitz-det");
	warn("test-ntl-rr", "floating point equality");
	ntl_tests.insert("test-ntl-rr");
	ntl_tests.insert("test-ntl-sylvester");
	ntl_tests.insert("test-ntl-toeplitz");
	ntl_tests.insert("test-ntl-zz_p");

        // are these really ntl dependent?
//	ntl_tests.insert("test-smith-form");
//	ntl_tests.insert("test-smith-form-adaptive");
//	ntl_tests.insert("test-smith-form-iliopoulos");
	ntl_tests.insert("test-polynomial-local-x");
	ntl_tests.insert("test-weak-popov-form");
	ntl_tests.insert("test-frobenius-leading-invariants");
	ntl_tests.insert("test-frobenius-large");
	ntl_tests.insert("test-invariant-factors");
	ntl_tests.insert("test-frobenius-small");
	ntl_tests.insert("test-poly-smith-form");

	set< string> ocl_tests;
	ocl_tests.insert("test-opencl-domain");

    set<string> mpi_tests;
    mpi_tests.insert("test-mpi-comm");
//// Things are automatic from here onward. ////

        // process optional dependencies
#ifndef LINBOX_HAVE_OCL
    for (set< string>::iterator i = ocl_tests.begin(); i != ocl_tests.end(); ++i)
		skip(*i, "OpenCL not present");
#endif
#ifndef __LINBOX_HAVE_NTL
    for (set< string>::iterator i = ntl_tests.begin(); i != ntl_tests.end(); ++i)
		skip(*i, "NTL not present");
#endif
#ifndef LINBOX_HAVE_MPI
    for (set< string>::iterator i = mpi_tests.begin(); i != mpi_tests.end(); ++i)
		skip(*i, "MPI not present");
#endif
// build and run the tests section
	string t, cmd;
	if ( system("ls test-*C >.tmp-tests") ) {
              std::cerr << "FAIL. Problem finding test sources.\n";
		exit(-1);
	}
	ifstream tests(".tmp-tests");
	int buildfail = 0, runfail = 0, skipped = 0, pass = 0; // counts
	vector<string> all_tests;
	while (tests >> t) {t.resize(t.size()-2); all_tests.push_back(t);}
//#ifdef __LINBOX_USE_OPENMP
//#pragma omp parallel for
//PARFOR1D(i, 0, all_tests.size(), SPLITTER(),
//#endif
	for (size_t i = 0; i < all_tests.size(); ++i)
	{
		t = all_tests[i];
		if (hides.count(t) > 0) continue; // ignore entirely
		stringstream report;
		report.width(35);
		report << left << t;
		if (skips.count(t) > 0) { // print skip message
			report << "skipped, " + warn_note[t];
			skipped++;
		} else { // do build
                        int status(0);
			if (force_build) {
				cmd = "touch " + t + ".C";
				status += system(cmd.c_str());
			}
                // build
			cmd = "make " + t + " 2> /dev/null > /dev/null";
			report << "build ";
			status += system(cmd.c_str());
//			#ifndef __LINBOX_USE_OPENMP
//			if (status == 2) break;
//			#endif
			if (status != 0) { // build failure
				buildfail++;
				report << "FAIL. ";
			} else {// do run
				report << "OK, run ";
				std::ostringstream prog ;
				prog << "./" << t ;
				status = system(prog.str().c_str());
//				#ifndef __LINBOX_USE_OPENMP
//				if (status == 2) break;
//				#endif
				if (status == 0) {
					report << "OK.  " + warn_note[t];
					pass++;
				} else {
					runfail++;
					report << "FAIL. " + warn_note[t];
				}
			}
		}
		if (reporting) cout << report.str() << endl;
	} // for i
        //); // parfor

        //// summary ////
	cout << "--------------------------------------------------------------------" << endl;
	int total = pass + buildfail + runfail + skipped;
	if (buildfail || runfail)
		cout << "Of " << total << " tests, "
             << pass << " passed, "
             << buildfail << " failed to compile, "
             << runfail << " failed at runtime, "
             << " and " << skipped << " were skipped." << endl;
	else
		cout << endl << "All " << total - skipped << " tests pass."
             << "  (There were " << skipped << " skipped tests.)" << endl;
	cout << "--------------------------------------------------------------------" << endl;

	return buildfail || runfail ? -1 : 0;
} // main

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
