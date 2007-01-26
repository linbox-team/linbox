/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/**\file examples/det.C examples/det.C
\brief Determinant of sparse matrix over Z or Zp.
\ingroup examples
*/

//#include "linbox-config.h"
#include <iostream>
#include <string>

#include "linbox/field/modular-double.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/det.h"
#include "linbox/util/matrix-stream.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
	if (argc > 3) {
		cerr << "Usage: det <matrix-file-in-supported-format> [<p>]" << endl;
		return -1;
	}
	//  CASE:  matrix file piped into this program
	//  NOTE:  parallel compuation not yet working for this method
	if (argc <= 1 || (argc == 2 && !strcmp(argv[1], "-p")) ) {
		// For a large integer matrix test, do "bigmat <n> | det", 
		// where <n> is a size parameter of your choice.

		int process = 0;
#ifdef __LINBOX_HAVE_MPI
   		Communicator C(&argc, &argv);
			process = C.rank();
#endif

			GMP_Integers ZZ;
			GMP_Integers::Element det_A;
			SparseMatrix<GMP_Integers> A(ZZ);
      if(argc <= 1 || process == 0){
         // std::cout << "# of processes = " << C.size() << std::endl;
			A.read(cin);
			cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
       }
			//GMP_Integers::Element det_A;
			
			if(argc == 2);
				//det(det_A, A, C);
			else
				det (det_A, A);

			cout << "Determinant is ";
			ZZ.write(cout, det_A) << endl;
		
	}
	//  CASE:  det called using a matrix file name on the cmd line
	//  ex:  ./det [matrix-file] [-p]
	else if (argc == 2 || (argc == 3 && !strcmp(argv[2], "-p")) ) {
	
		// For a small integer matrix test, do "det data/mat2.txt". 
		// It is a 2 by 2 matrix with determinant = -2.

		//  set up parallel code object
	   Communicator *Cptr = NULL;
#ifdef __LINBOX_HAVE_MPI
		if (argc == 3)
			Cptr = new Communicator(&argc, &argv);
#endif

		GMP_Integers ZZ;
		ifstream input (argv[1]);
		if (!input) 
		{ cerr << "Error opening matrix file " << argv[1] << endl; 
			return -1; 
		}

		SparseMatrix<GMP_Integers> A(ZZ);
		A.read(input);
#ifdef __LINBOX_HAVE_MPI
		if(argc == 2 || !Cptr->rank())
#endif
		   cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		GMP_Integers::Element det_A;

#ifdef __LINBOX_HAVE_MPI
		//  call parallel det with cra if applicable
		if(argc == 3)
		   cra_det (det_A, A, RingCategories::IntegerTag(), Method::Hybrid(*Cptr), Cptr);
		else
#endif
			det(det_A, A);

#ifdef __LINBOX_HAVE_MPI
		//  if not using parallel or if parent process in
		//  parallel, report the determinant
      if(argc == 2 || Cptr->rank() == 0){ 
#endif
			cout << "Determinant is ";
			ZZ.write(cout, det_A) << endl;
#ifdef __LINBOX_HAVE_MPI
		}
#endif

#ifdef __LINBOX_HAVE_MPI
		//  tie up parallel loose ends if necessary
		if(argc == 3)
			MPI_Finalize();
#endif
	}

	else if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		ifstream input (argv[1]);
		if (!input) 
		{ cerr << "Error opening matrix file " << argv[1] << endl; 
		  return -1; 
		}
		MatrixStream< Field > ms ( F, input );
		SparseMatrix<Field> B (ms);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

		Field::Element det_B;
		det (det_B, B);

		cout << "Determinant is ";
		F.write(cout, det_B) << " mod " << q << endl;
	}

	else{
		cout << "HELLO WORLD" << endl;
	}

	return 0;
}
