
/*
 * examples/mpidet.C
 *
 * Copyright (C) 2006, 2010 B Youse, D Saunders
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

/**\file examples/mpidet.C
 * @example examples/mpidet.C
  \brief Determinant of sparse matrix over Z or Zp.
  \ingroup examples
  */

#include <iostream>
#include <string>

#include <linbox/field/modular.h>
#include <linbox/matrix/sparse.h>
#include <linbox/solutions/det.h>
#include <linbox/util/matrix-stream.h>

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
#ifdef __LINBOX_HAVE_MPI
	if (argc < 2) {
		cerr << "Usage: det <matrix-file-in-supported-format>" << endl;
		return -1;
	}
	//  ex:  ./det [matrix-file]
	else{
		// For a small integer matrix test, do "make mpidet2 -f makefile.mpi"

		//  set up parallel code object
		Communicator *Cptr = NULL;
		Cptr = new Communicator(&argc, &argv);

		typedef PID_integer Integers;
		Integers ZZ;

		ifstream input (argv[1]);
		if (!input)
		{ cerr << "Error opening matrix file " << argv[1] << endl;
			return -1;
		}

		SparseMatrix<Integers>A(ZZ);
		A.read(input);
		if(!Cptr->rank()){
			cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
			cout << "Beginning parallel computation with " << Cptr->size()
			<< " processes." << endl;
		}

		Integers::Element det_A;

		//  call parallel det with cra
		cra_det(det_A, A, RingCategories::IntegerTag(), Method::Hybrid(*Cptr), Cptr);

		//  if parent process, report the determinant
		if(!Cptr->rank()){
			cout << "Determinant is ";
			ZZ.write(cout, det_A) << endl;
		}
		//  tie up parallel loose ends if necessary
		MPI_Finalize();
	}
	return 0;
#else
	cerr << "Compile with -D__LINBOX_HAVE_MPI" << endl;
	return -1 ;
#endif
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

