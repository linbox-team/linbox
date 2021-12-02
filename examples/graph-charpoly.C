/*
 * examples/graph-charpoly.C
 *
 * Copyright (C) 2005, 2007, 2010 C. Pernet
  ========LICENCE========
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

/*! @file examples/graph-charpoly.C
  * @example examples/graph-charpoly.C
 * @ingroup examples
 * @brief Compute the integer characteristic polynomial of symetric matrices with (0,1) coefficients.
 * Example for Pr G. Royle.
 */

#include <linbox/linbox-config.h>

#include <iostream>

#include <linbox/blackbox/zero-one.h>
#include <givaro/zring.h>
#include <linbox/solutions/charpoly.h>
#include <linbox/solutions/methods.h>

using namespace std;
using namespace LinBox;

template <class Field, class Polynomial>
void printPolynomial(const Field& F, const Polynomial& P)
{
	int n= (int) P.size()-1;
	if (n==1){
		cout<<"X";
		if ( P[0] != 0)
			F.write(cout<<((P[0]>0)?"+":""),P[0]);
	}
	else{
		cout<<"X^"<<n;
		for ( int i=n-1; i>1; --i)
			if (!F.isZero(P[(size_t)i]))
				F.write(cout<<((P[(size_t)i]>0)?"+":""),P[(size_t)i])<<"*X^"<<i;
		if ( P[1] != 0)
			F.write(cout<<((P[1]>0)?"+":""),P[1])<<"*X";
		if ( P[0] != 0)
			F.write(cout<<((P[0]>0)?"+":""),P[0]);
	}
}


int main (int argc, char **argv)
{
	commentator().getMessageClass (BRIEF_REPORT).setMaxDepth (2);
	commentator().getMessageClass (BRIEF_REPORT).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);


	if (argc != 2) {
		cerr << "Usage: graph-charpoly <0/1-symmetric-matrix-file-in-SMS-format>" <<endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) {
		cerr << "Error opening matrix file " << argv[1] << endl;
		return -1;
	}

    typedef Givaro::ZRing<Integer> IRing_t;

	IRing_t ZZ;
	ZeroOne<IRing_t > A(ZZ);
	A.read (input);
	commentator().report(1, BRIEF_REPORT)<< "A is " << A.rowdim() << " by " << A.coldim() << endl;

    DensePolynomial<IRing_t> c_A(ZZ);

	charpoly (c_A, A, Method::Blackbox(Method::Wiedemann(Shape::Symmetric)));

	PolynomialRing<IRing_t>(ZZ,'X').write(
        cout<< "Characteristic Polynomial is ", c_A) << std::endl;


	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
