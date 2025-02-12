/* Copyright (C) LinBox
 *
 *  Author: bds
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file tests/test-smith-form.h
 * @ingroup tests
 * @brief tools for making matrix with known SNF.
 */

#ifndef __TEST_SMITH_FORM_H
#define __TEST_SMITH_FORM_H

#include "linbox/util/commentator.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/solutions/smith-form.h"

using std::endl;
using namespace LinBox;

template <class PIR> // This is for PIR = Z or Z/nZ.
BlasVector<PIR> & makeBumps(BlasVector<PIR> & b, int choice) {
	const PIR & R = b.field();
	typename PIR::Element two, three, nine, x;
	R.init(two,2);
	R.init(three,3);
	R.init(nine,9);
	R.init(x,202);
	// b is a single row
	size_t n = b.size();
	switch (choice) {
		case 0: // all zero
				for(size_t i = 0; i < n; ++i) b.setEntry(i,R.zero);
					break;
		case 1: // identity
				for(size_t i = 0; i < n; ++i) b.setEntry(i,R.one);
					break;
		case 2: // powers of 2
				for(size_t i = 0; i < n; ++i) b.setEntry(i,two);
					break;
		case 3: // random followed by 202,0.  Random part is largely 1's.
				for(size_t i = 0; i < n-2; ++i) {
					size_t r = rand()%20;
					if (r < 17) b.setEntry(i,R.one);
					if (r == 17) b.setEntry(i,two);
					if (r == 18) b.setEntry(i,three);
					if (r == 19) b.setEntry(i,nine);
				}
				b.setEntry(n-2,x);
				b.setEntry(n-1,R.zero);
	}
	// negate a few
	for (size_t k = rand()%4; k > 0; --k){
		size_t i = rand()%n;
		b.getEntry(x,i);
		b.setEntry(i,R.negin(x));
	}

	return b;
}


// For any PIR, build an increasing sequence of smith invariants d from "bumps" b.
template <class PIR>
BlasVector<PIR> & prefixProduct (BlasVector<PIR> & d, const BlasVector<PIR> & b) {
	const PIR& R = d.field();
	typename PIR::Element x,y; R.init(x); R.init(y);
	d.setEntry(0,b.getEntry(x,0));
	for (size_t i = 1; i < d.size(); ++i){
		d.getEntry(x,i-1);
		b.getEntry(y,i);
		d.setEntry(i,R.mulin(x, y));
	}
	return d;
}

// Generate A with snf = diag(d) (up to sign), based on the bumps.
// Think of bumps[i] as s_i/s_{i-1}, quotient of smith invariants.
// The lumps are used for off diagonal entries in L,U (triangular scramblers).
template <class PIR>
void makeSNFExample(DenseMatrix<PIR>& A,
					BlasVector<PIR> & d,
			  const BlasVector<PIR> & bumps,
			  const BlasVector<PIR> & lumps) {
	//LinBox::VectorWrapper::ensureDim (d, bumps.size());
	//LinBox::VectorWrapper::ensureDim (d, std::min(A.rowdim(), A.coldim()));
	prefixProduct(d, bumps);

    //    std::cout<<"d="<<d<<std::endl;
	// make A = UDL for random unimodular L,U
	const PIR & R = A.field();
	DenseMatrix<PIR> L(R, A.coldim(), A.coldim()),
					U(R, A.rowdim(), A.rowdim());
	typename PIR::Element x; R.init(x);
	size_t i, j, k;
	k = lumps.size();
	A.zero();
	for(i = 0; i < d.size(); ++i) A.setEntry(i,i,d.getEntry(x,i));

    //A.write(std::cout);
    //std::cout<<"PPP\n";
	L.zero();
	for(i = 0; i < L.rowdim(); ++i) L.setEntry(i,i,R.one);
	for (i = 0; i < L.rowdim(); ++ i)
		for (j = 0; j < i; ++ j) L.setEntry(i,j,lumps[rand()%k]);

	U.zero();
	for(i = 0; i < U.rowdim(); ++i) U.setEntry(i,i,R.one);
	for (i = 0; i < U.rowdim(); ++ i)
		for (j = i+1; j < U.coldim(); ++ j) U.setEntry(i,j,lumps[rand()%k]);



	// A <- UAL
	BlasMatrixDomain<PIR> MD(R);
	MD.mulin_left(A,L);
	MD.mulin_right(U,A);

	for (i = 0; i < d.size(); ++ i)
		d.setEntry(i,R.abs(x, d.getEntry(x,i)));
	// Now A is matrix equivalent to diag prefix product of bumps.
	// Now d is SNF diagonal (vector of invariants) for A.

	std::ostream & report =
#ifndef DISABLE_COMMENTATOR
        commentator().report()
#else
        std::clog
#endif
;
    A.write(report << "Created: ", Tag::FileFormat::linalg) << std::endl;
}

	std::ostream & report =
#ifndef DISABLE_COMMENTATOR
        commentator().report(LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
#else
        std::clog;
#endif
;

template <class PIR>
bool checkSNFExample( const BlasVector<PIR>& d, const BlasVector<PIR>& x ){

	VectorDomain<PIR> VD(d.field());

	report << "Expected smith form: ";
	VD.write (report, d) << endl;

	report << "Computed Smith form: ";
	VD. write (report, x) << endl;

	if (not VD.areEqual (d, x)) {
		report << "ERROR: Computed not as Expected." << endl;
		return false;
	}

    report << "PASSED." << endl;

    return true;
}

template <class PIR>
bool checkSNFExample( const LinBox::SmithList<PIR>& d, const LinBox::SmithList<PIR>& x, const PIR& R){

/*
#ifndef DISABLE_COMMENTATOR
	std::ostream & report =
      commentator().report(LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
#else
	std::ostream & report = std::clog;
#endif
*/
	report << "Expected smith form SL: " << '{';
    for(auto const & sit: d) report << '{' << sit.first << ',' << sit.second << '}';
    report << '}' << std::endl;

	report << "Computed Smith form SL: " << '{';
    for(auto const & sit: x) report << '{' << sit.first << ',' << sit.second << '}';
    report << '}' << std::endl;


    auto dit=d.begin();
    auto xit=x.begin();
    bool pass=true;
    for( ; dit != d.end(); ++dit, ++xit) {
        if (!R.areEqual(dit->first,xit->first)
            || dit->second!=xit->second) {
            R.write( R.write(
                report << "ERROR: Computed not as Expected. "
                << '(', dit->first) << ',' << dit->second << ')'
                       << " != "
                       << '(', xit->first) << ',' << xit->second << ')'
                       << endl;
            pass = false;
        }
    }

    report << "PASSED." << endl;

    return pass;
}
#endif // __TEST_SMITH_FORM_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
