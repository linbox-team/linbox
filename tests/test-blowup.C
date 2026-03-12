/* tests/test-blowup.C
 * Copyright (C) 2023 BDS
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file  tests/test-blowup.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <cstdio>

#include "linbox/blackbox/blowup.h"
#include "linbox/util/commentator.h"
//#include "linbox/field/archetype.h"
#include "linbox/ring/modular.h"
#include <givaro/givranditer.h>
//***
#include "linbox/matrix/sparse-matrix.h"
//#include "linbox/ring/ntl.h"
#include "linbox/ring/polynomial-ring.h"
#include "linbox/solutions/rank.h"
//***

#include "test-blackbox.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 7;
	static size_t d = 3;
	static integer q = 65521U;
	static int iterations = 1; // was 100

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set deg of modulus to D.", TYPE_INT,     &d },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
//cout << "n " << n << ", d " << d << ", q " << q << endl;

   using Field = Givaro::Modular<double>;
	Field F (q);
   VectorDomain<Field> VD(F);

	srand ((unsigned)time (NULL));
   Field::RandIter iter(F);

	commentator().start("Blowup matrix black box test suite", "sparsemat");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

//* A is D
   SparseMatrix<Field> D(F, n, n);
	for (size_t i = 0; i < n; i++) D.setEntry(i, i, F.one);
   D.setEntry(1,0,  0 ); 
   D.setEntry(2,1,  0 ); 
   D.setEntry(3,2,  1 ); 
   D.setEntry(4,3,  0 );
   D.setEntry(5,4,  1 );
   D.setEntry(6,5,  1 );

cout << "D:"; double tmp;
for (size_t i = 0; i < n; i++) cout << D.getEntry(tmp, i, i) << ", ";
cout << endl;
for (size_t i = 1; i < n; i++) cout << D.getEntry(tmp, i, i-1) << ", ";
cout << endl;

//* f
/* f = x + 1                : (1,1)
   f^2 = x^2 + 2x + 1       : (1,2,1)
   f^3 = x^3 + 3x^2 + 3x + 1: (1,3,3,1)
   f^4 = (x+1)^4            : (1,4,6,4,1)
   f^5 = (x+1)^5            : (1,5,10,5,1)
*/

for (size_t d = 1; d < 5; ++d) {
   BlasVector<Field> f(F,d+1), s(F, d+2); // polynomial
   f[0] = 1; // i = 0 case
   s[1] = 1;
	for (size_t i = 1; i <= d; ++i) { // build (x+1)^i
	   for (size_t j = 1; j <= i; ++j) f[j] += s[j];
	   for (size_t j = 1; j <= i+1; ++j) s[j] = f[j-1];
	}
   //f[d] = F.one;

cout << "f:";
for (size_t i = 0; i < d+1; i++) cout << f[i] << ", ";
cout << endl;
//*

	LinBox::Blowup<SparseMatrix<Field> > B(D, f);

   size_t r;
   LinBox::rank< Blowup<SparseMatrix<Field> >(r, B);
   cout << "n " << n << ". rank " << r << ", de " << d << ", d 1" << ", q " << q << endl;
} // for d
#if 0
// my initial tests
   BlasVector<Field> u(F,n*d), v(F, n*d), Av(F, n*d), uAT(F, n*d);
	for (size_t i = 0; i < n*d; i++) {
		iter.random(u[i]);
		iter.random(v[i]);
	}
	for (size_t i = 0; i < n; i += 1) {
      //u[i] = F.one;
      u[i*d] = i;
	   for (size_t j = 1; j < n; j += 1) 
         u[i*d + j] = 1;
	}

cout << "  u:     ";
for (size_t i = 0; i < n*d; i++) cout << u[i] << ", "; cout << endl;
   v.clear();
   B.apply(v, u);
cout << "  v:     ";
for (size_t i = 0; i < n*d; i++) cout << v[i] << ", "; cout << endl;
#endif
#if 0  // transpose test
cout << "  u:     ";
for (size_t i = 0; i < n*d; i++) cout << u[i] << ", "; cout << endl;
cout << "  v:     ";
for (size_t i = 0; i < n*d; i++) cout << v[i] << ", "; cout << endl;
   Av.clear();
   B.apply(Av, v);

   uAT.clear();
   B.applyTranspose(uAT, u);

   Field::Element r1, r2;
   VD.dot (r1, uAT, v); 
   VD.dot(r2, u, Av); 
   cout << "(uA^T)^Tv " << r1 << ", u(Av) " << r2;
   if (not F.areEqual(r1, r2)) { pass = false; cout << " fail" << endl; }
   else cout << " pass " << endl;
if (not pass) {
cout << "uAT:     ";
for (size_t i = 0; i < n*d; i++) cout << uAT[i] << ", "; cout << endl;
cout << " Av:     ";
for (size_t i = 0; i < n*d; i++) cout << Av[i] << ", "; cout << endl;
}

/*
	for (size_t i = 0; i < n*d; ++i) {
      if (v[i] != w[i]) {pass = false; break;}
   }
cout << "u:     ";
for (size_t i = 0; i < n*d; i++) cout << u[i] << ", ";
cout << endl;
cout << "v = Au: ";
for (size_t i = 0; i < n*d; i++) cout << v[i] << ", ";
cout << endl;
cout << "w=A^Tu: ";
for (size_t i = 0; i < n*d; i++) cout << w[i] << ", ";
cout << endl;
	for (size_t i = 0; i < n*d; i++) 
      if (w[i]!= v[i]) {pass = false; break;}

*/
   if (not pass) cout << "initial test fail" << endl;
// end initial tests
#endif

//   pass = pass && testBlackboxNoRW(B);

	commentator().stop (MSG_STATUS (pass));

   cout << (pass ? "pass" : "fail") << endl;
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
