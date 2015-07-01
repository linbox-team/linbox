/* Copyright 2013 (C) the LinBox group
 *
 * Written by :
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


#include <cstdlib>
#include <stdio.h>

int main(){
	bool pass = true;
	int d;
	FILE * infile = fopen("tutorial-1.in", "w");
	fprintf(infile, "2 2\n");
	fprintf(infile, "3 1\n");
	fprintf(infile, "1 2\n");
	fclose(infile);

	system("./tutorial-1 <tutorial-1.in >tutorial-1.out");

	FILE * outfile = fopen("tutorial-1.out", "r");
	fscanf(outfile, "the determinant is %d", &d);
	fclose(outfile);

	system("rm tutorial-1.in tutorial-1.out");

	pass = pass and d == 5;
	return pass? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

