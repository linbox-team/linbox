/*
 * examples/checksolve.C
 *
 * Copyright (C) 2007 C. Pernet
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
/**\file examples/checksolve.C
 * @example  examples/checksolve.C
  \brief Solving of sparse matrix over Z or Zp.
  \ingroup examples
  */
#include <iostream>

#include <linbox/field/modular.h>
#include <linbox/solutions/solve.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/solutions/methods.h>

using namespace LinBox;
// using namespace std;


int main (int argc, char **argv)
{

	// 	commentator().setMaxDetailLevel (-1);
	// 	commentator().setMaxDepth (-1);
	// 	commentator().setReportStream (std::cerr);


	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: checksolve <matrix-file-in-supported-format> <dense-vector-file> <p>" << std::endl;
		return 0;
	}


	std::ifstream input (argv[1]);
	if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }

	std::ifstream invect(argv[2]);
	if (!input) { std::cerr << "Error opening vector file " << argv[2] << std::endl; return -1; }


	typedef Modular<double> Field;
	double q = atof(argv[3]);
	Field F(q);
	MatrixStream< Field > ms ( F, input );
	BlasMatrix<Field> A (ms); // A.write(std::cout);
	std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

	BlasVector<Field> X(F, A.coldim()),B(F, A.rowdim());

	for(BlasVector<Field>::iterator it=B.begin();
	    it != B.end(); ++it)
		invect >> *it;

	std::cout << "B is [ "<<B<< "]" << std::endl;

	solve (X, A, B, Method::BlasElimination());

	std::cout << "(BlasElimination) Solution is [ "<<X<< "]" << std::endl;
	BlasVector<Field> r(F, A.rowdim());
	BlasMatrixDomain<Field> BMD(F);
	BMD.mul(r, static_cast<BlasMatrix<Field>& >(A), X);
	//A.apply (r,X);
	VectorDomain<Field> VD(F);
	if (VD.areEqual (r,B))
		std::cout<<"CHECK"<<std::endl;
	else{
		std::cout<<"FAIL"<<std::endl;
		std::cout<<"r = "<<r<<std::endl;
		std::cout<<"B = "<<B<<std::endl;
	}

	return 0;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

