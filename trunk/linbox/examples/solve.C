/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/**\file examples/solve.C
\brief Solving of sparse matrix over Z or Zp.
\ingroup examples
*/
//#include "linbox-config.h"
#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/field/PID-integer.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/methods.h"

using namespace LinBox;

int main (int argc, char **argv)
{

	commentator.setMaxDetailLevel (-1);
	commentator.setMaxDepth (-1);
	commentator.setReportStream (std::cerr);


	if (argc < 2 || argc > 4) {
		cerr << "Usage: solve <matrix-file-in-supported-format> [<dense-vector-file>] [<p>]" << endl;
		return -1;
	}

	std::ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc == 2) {
		PID_integer ZZ;
		MatrixStream< PID_integer > ms( ZZ, input );
		SparseMatrix<PID_integer> A (ms);
		std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

		std::vector<PID_integer::Element> X( A.coldim()),B(A.rowdim());
                for(std::vector<PID_integer::Element>::iterator it=B.begin();
                    it != B.end(); ++it)
                    if (rand() <0.5)
                        *it = -1;
                    else
                        *it = 1;
		std::cout << "B is [";
                for(std::vector<PID_integer::Element>::const_iterator it=B.begin();
                    it != B.end(); ++it)
                    ZZ.write(cout, *it) << " ";
                std::cout << "]" << std::endl;
                

		 Method::BlockLanczos BLz;
		 BLz.checkResult    (false) ;

		Timer chrono; chrono.start();
		solve (X, A, B, BLz);
		chrono.stop();

		std::cout << "Solution is [";
                for(std::vector<PID_integer::Element>::const_iterator it=X.begin();
                    it != X.end(); ++it)
                    ZZ.write(cout, *it) << " ";
                std::cout << "]" << std::endl;
		
		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;
	}
// 	if (argc == 3) { 

// 		typedef Modular<double> Field;
// 		double q = atof(argv[2]);
// 		Field F(q);
// 		MatrixStream< Field > ms ( F, input );
// 		SparseMatrix<Field> B (ms);
// 		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

// 		Field::Element det_B;
// 		det (det_B, B);

// 		cout << "Determinant is ";
// 		F.write(cout, det_B) << " mod " << q << endl;
// 	}

	return 0;
}
