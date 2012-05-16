
/*
 * examples/echelon.C
 *
 * Copyright (C) 2007, 2010 C. Pernet
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

/** \file examples/echelon.C
 * @example examples/echelon.C
  \brief Echelon form of matrix over  Zp.
  \ingroup examples
  */

#include <iostream>

#include "linbox/field/modular.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/algorithms/echelon-form.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
	if (argc < 3 || argc > 3) {
		cerr << "Usage: echelon <matrix-file-in-SMS-format> <p>" << endl;
		return -1;
	}
	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	typedef Modular<double> Field;
	double q = atof(argv[2]);
	Field F(q);

	MatrixStream<Field> ms(F,input);
	BlasMatrix<Field> A(ms);
	BlasMatrix<Field> E(F,A.rowdim(),A.coldim());
	cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

	EchelonFormDomain<Modular<double> > EFD (F);

	EFD.rowReducedEchelon(E,A);

	if (E.coldim() <20)
		E.write(cerr<<"Echelon = "<<endl)<<endl;

	return 0;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

