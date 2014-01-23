/*
 * examples/athadet.C
 *
 * Copyright (C) 2007, 2010 S. Guelton
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

/**\file examples/athadet.C
 * @example  examples/athadet.C
  \brief Determinant of sparse matrix over Z or Zp.
  \ingroup examples
  @author serge.guelton@imag.fr
  */


/*
 * to be  compiled with kaapi library **svn sources**, available at http://gforge.inria.fr
 * g++ -o atahdet athadet.C  $(pkg-config --cflags --libs kaapi) -g3 $( linbox-config --cflags --libs)
 * or
 * g++ -o atahdet athadet.C  $KAAPI_CPPFLAGS $KAAPI_LDFLAGS -g3 $( linbox-config --cflags --libs)
 */
#define __LINBOX_HAVE_KAAPI

#include <iostream>
#include <string>

#include "linbox/field/modular.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/matrix/sparse.h"
#include "linbox/solutions/det.h"
#include "linbox/util/matrix-stream.h"


using namespace LinBox;
using namespace std;

typedef PID_integer Integers;

struct  cra_det_task {


	void operator()(int argc, char **argv) {
		ifstream input (argv[1]);
		if (!input) {
			cerr << "Error opening matrix file '" << argv[1] << "'" << endl;
			return ;
		}


		Util::logfile() << "creating matrix" << std::endl;
		Integers ZZ;
		SparseMatrix2<Integers> sparseMatrix (ZZ);
		sparseMatrix.read(input);
		Util::logfile() << "matrix created" << std::endl;
		PID_integer::Element det_A;
		try {
			cra_det(det_A, sparseMatrix, RingCategories::IntegerTag(), Method::Hybrid() );
		}
		catch ( LinBox::LinboxError & err ) {
			std::cerr << err << std::endl;
		}
		cout << "Determinant is ";
		ZZ.write(cout, det_A) << endl;
	}
};

int main (int argc, char **argv)
{
	if(argc<2) {
		cerr << "not enough args : usage = athadet file" << endl;
		return 1;
	}

	a1::Community com = a1::System::join_community(argc,argv);

	Util::logfile() << "starting job" << std::endl;
	a1::ForkMain<cra_det_task>()(argc,argv);
	com.leave();
	a1::System::terminate();
	return 0;
}


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

