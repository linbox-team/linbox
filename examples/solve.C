/*
 * examples/solve.C
 *
 * Copyright (C) 2005, 2010 J-G Dumas, D. Saunders, P. Giorgi
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

/** @file examples/solve.C
 * @ingroup examples
 * @brief Blackbox solvers.
 * @warning some are commented out...
 * @example  examples/solve.C
 */

#include <linbox/linbox-config.h>

#include <iostream>

#include <givaro/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/algorithms/gauss.h>
#include <linbox/util/timer.h>

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{


	if (argc != 4) {
		cerr << "Usage: solve <matrix-file-in-supported-format> <output-file> <p>" << endl;
		return 0;
	}
	srand48( BaseTimer::seed() );

	std::ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	std::ofstream nsb;
	nsb.open (argv[2], std::ofstream::out);
	if (!nsb) { cerr << "Error opening nullspace output file " << argv[2] << endl; return -1; }


        cout<<"Computation is done over Z/("<<atoi(argv[3])<<")"<<endl;
		typedef Givaro::Modular<int64_t> Field;
		double q = atof(argv[3]);
		typedef DenseVector<Field> DenseVector ;
		Field F(q);
		MatrixStream< Field > ms ( F, input );
		SparseMatrix<Field> A (ms);  // A.write(std::cout);
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
                if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cerr << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
		DenseMatrix<Field> N(F, A.rowdim(), 15);
		Timer chrono;

		// Sparse Elimination
		chrono.clear();
		chrono.start();
		GaussDomain<Field> GD ( A.field() );
		GD.nullspacebasisin(N, A);

		chrono.stop();

		N.write(nsb) << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<<std::endl;;

#if 0
		// Wiedemann
		std::cout << "Blackbox" << std::endl;
		chrono.clear();
		chrono.start();
		solve (X, A, B, Method::Blackbox());
		chrono.stop();

		std::cout << "(Wiedemann) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			F.write(cout, *it) << " ";
		std::cout << "]" << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<<std::endl;;
#endif
#if 0
		// Lanczos
		std::cout << "Lanczos" << std::endl;
		chrono.clear();
		chrono.start();
		solve (X, A, B, Method::Lanczos());
		chrono.stop();

		std::cout << "(Lanczos) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			F.write(cout, *it) << " ";
		std::cout << "]" << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<< std::endl;


		// Block Lanczos
		std::cout << "Block Lanczos" << std::endl;
		Method::BlockLanczos MBL;
		MBL.preconditioner(Specifier::FULL_DIAGONAL);
		chrono.clear();
		chrono.start();
		solve (X, A, B, MBL);
		chrono.stop();

		std::cout << "(Block Lanczos) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			F.write(cout, *it) << " ";
		std::cout << "]" << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<< std::endl;
#endif

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
