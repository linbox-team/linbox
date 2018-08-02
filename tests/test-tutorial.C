/* Copyright 2013 (C) the LinBox group
 *
 * updated in 2018 by bds
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

/*! @file  tests/test-tutorial.C
 * @ingroup tests
 *
 * @brief no doc
 *
 * @test no doc.
 */
#include <fstream>
struct NullBuffer : public std::streambuf
{ int overflow(int c) { return c; } };
NullBuffer null_buffer;

std::ostream nullout(&null_buffer);
std::ifstream smsin("data/sms.matrix"); 

/* Each time the tutorial is changed, this test of it should be updated by 
   including the tutorial example below (the include's and all of main).
   Make two changes:

   2. replace "std::cin" with "smsin".
   2. replace "std::cout" with "nullout"
*/

// Begin inclusion of tutorial 
#include <givaro/modular.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/det.h>

using namespace LinBox;

int main()
{
	typedef Givaro::Modular<double> Field;
	Field F(101);

	SparseMatrix<Field> A(F);
	//A.read(std::cin);
	A.read(smsin);

	Field::Element d;
	Method::SparseElimination M;

	det(d, A, M);

	//F.write(std::cout << "the determinant is ", d) << std::endl;
	F.write(nullout << "the determinant is ", d) << std::endl;
	return 0; // if it runs, it passes.
}
// End inclusion of tutorial 
   

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
