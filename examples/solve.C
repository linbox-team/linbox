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
#include <givaro/zring.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/solve.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/solutions/methods.h>

using namespace LinBox;

int main (int argc, char **argv)
{

	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::cerr);


	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: solve <matrix-file-in-supported-format> [<dense-vector-file>] [<p>]" << std::endl;
		return 0;
	}
	srand48( BaseTimer::seed() );

	std::ifstream input (argv[1]);
	if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }
	std::ifstream invect;

	bool createB = false;
	int ModComp = 0;
	if (argc == 2) {
		createB = true;
		ModComp = 0;
	}

	if (argc == 3) {
		invect.open (argv[2], std::ifstream::in);
		if (!invect) {
			createB = true;
			ModComp = 2;
		}
		else {
			createB = false;
			ModComp = 0;
		}
	}

	if (argc == 4) {
		ModComp = 3;
		invect.open (argv[2], std::ifstream::in);
		if (!invect) {
			createB = true;
		}
		else
			createB = false;
	}


	if (ModComp) {
        std::cout<<"Computation is done over Z/("<<atoi(argv[ModComp])<<")"<<std::endl;
		typedef Givaro::Modular<double> Field;
		double q = atof(argv[ModComp]);
		typedef DenseVector<Field> DenseVector ;
		Field F(q);
		MatrixStream< Field > ms ( F, input );
		SparseMatrix<Field> A (ms);  // A.write(std::cout);
		std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
        if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cerr << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
		DenseVector X(F, A.coldim()),B(F, A.rowdim());
		if (createB) {
			std::cerr << "Creating a random {-1,1} vector U, B is AU (to have a consistent system)" << std::endl;
			DenseVector U(F, A.coldim() );
			for(DenseVector::iterator it=U.begin();
                it != U.end(); ++it)
				if (drand48() <0.5)
					F.assign(*it,F.mOne);
				else
					F.assign(*it,F.one);
			A.apply(B,U);
		}
		else {
			for(DenseVector::iterator it=B.begin();
                it != B.end(); ++it)
				F.read(invect,*it);
		}

		// A.write(std::cout << "A: ") << std::endl;

		std::cout << "B is " << B << std::endl;

		Timer chrono;

		// Sparse Elimination
		std::cout << "Sparse Elimination" << std::endl;
		chrono.clear();
		chrono.start();
		Method::SparseElimination M;
		solve (X, A, B, M);
		chrono.stop();

		std::cout << "(Sparse Gauss) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			F.write(std::cout, *it) << " ";
		std::cout << "]" << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<<std::endl;;

		// DenseElimination
		std::cout << "DenseElimination" << std::endl;
		chrono.start();
		solve (X, A, B, Method::DenseElimination());
		chrono.stop();

		std::cout << "(DenseElimination) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			F.write(std::cout, *it) << " ";
		std::cout << "]" << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<< std::endl;

		// Wiedemann
		std::cout << "Blackbox" << std::endl;
		chrono.clear();
		chrono.start();
		solve (X, A, B, Method::Blackbox());
		chrono.stop();

		std::cout << "(Wiedemann) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			F.write(std::cout, *it) << " ";
		std::cout << "]" << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<<std::endl;;
#if 0
		// Lanczos
		std::cout << "Lanczos" << std::endl;
		chrono.clear();
		chrono.start();
		solve (X, A, B, Method::Lanczos());
		chrono.stop();

		std::cout << "(Lanczos) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			F.write(std::cout, *it) << " ";
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
			F.write(std::cout, *it) << " ";
		std::cout << "]" << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<< std::endl;
#endif

	}
	else {
		std::cout<<"Computation is done over Q"<<std::endl;
		Givaro::ZRing<Integer> ZZ;
		typedef DenseVector<Givaro::ZRing<Integer> > DenseVector ;
		MatrixStream< Givaro::ZRing<Integer> > ms( ZZ, input );
		SparseMatrix<Givaro::ZRing<Integer> > A (ms);
		Givaro::ZRing<Integer>::Element d;
		std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
                if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cerr << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
		DenseVector X(ZZ, A.coldim()),B(ZZ, A.rowdim());

		if (createB) {
			std::cerr << "Creating a random {-1,1} vector U, B is AU" << std::endl;
			DenseVector U(ZZ, A.coldim() );
			for(DenseVector::iterator it=U.begin();
			    it != U.end(); ++it)
				if (drand48() <0.5)
					*it = -1;
				else
					*it = 1;
			A.apply(B,U);
		}
		else {
			for(DenseVector::iterator it=B.begin();
			    it != B.end(); ++it)
				invect >> *it;
		}

		std::cout << "B is " << B << std::endl;

		Timer chrono;
		// DenseElimination
        std::cout << "DenseElimination" << std::endl;
        chrono.start();
        solve (X, d, A, B, Method::DenseElimination());
        chrono.stop();

 		std::cout << "(DenseElimination) Solution is [";
        for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
 		ZZ.write(std::cout, *it) << " ";
        std::cout << "] / ";
        ZZ.write(std::cout, d)<< std::endl;
        std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;

		// Sparse Elimination
		std::cout << "Sparse Elimination" << std::endl;
		chrono.start();
		solve (X, d, A, B, Method::SparseElimination());
		chrono.stop();

		std::cout << "(SparseElimination) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			ZZ.write(std::cout, *it) << " ";
		std::cout << "] / ";
		ZZ.write(std::cout, d)<< std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;

                		// Wiedemann
		std::cout << "Wiedemann" << std::endl;
		chrono.start();
		solve (X, d, A, B, Method::Wiedemann());
		chrono.stop();

		std::cout << "(Wiedemann) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			ZZ.write(std::cout, *it) << " ";
		std::cout << "] / ";
		ZZ.write(std::cout, d) << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;



#if 0
		// Lanczos
		std::cout << "Lanczos" << std::endl;
		chrono.start();
		solve (X, d, A, B, Method::Lanczos());
		chrono.stop();

		std::cout << "(Lanczos) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			ZZ.write(std::cout, *it) << " ";
		std::cout << "] / ";
		ZZ.write(std::cout, d) << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;


		// Block Lanczos
		std::cout << "Block Lanczos" << std::endl;
		chrono.clear();
		chrono.start();
		solve (X, d, A, B, Method::BlockLanczos());
		chrono.stop();

		std::cout << "(Block Lanczos) Solution is [";
		for(DenseVector::const_iterator it=X.begin();it != X.end(); ++it)
			ZZ.write(std::cout, *it) << " ";
		std::cout << "] / ";
		ZZ.write(std::cout, d) << std::endl;
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;
#endif
	}

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
