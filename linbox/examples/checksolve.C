/** 
 * examples/checksolve.C
 *
 * Copyright (C) 2007 C. Pernet
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */
/**\file examples/solve.C
\brief Solving of sparse matrix over Z or Zp.
\ingroup examples
*/
//#include "linbox-config.h"
#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/methods.h"

using namespace LinBox;
using namespace std;

template<class T, template <class ,class> class Container, template <class> class Alloc>
std::ostream& operator<< (std::ostream& o, const Container<T, Alloc<T> >& C) {
		for(typename Container<T, Alloc<T> >::const_iterator refs =  C.begin();
	    refs != C.end() ;
	    ++refs )
		o << (*refs) << " " ;
	return o;
}

int main (int argc, char **argv)
{

// 	commentator.setMaxDetailLevel (-1);
// 	commentator.setMaxDepth (-1);
// 	commentator.setReportStream (std::cerr);


    if (argc < 2 || argc > 4) {
        cerr << "Usage: checksolve <matrix-file-in-supported-format> <dense-vector-file> <p>" << endl;
        return 0;
    }

    
    std::ifstream input (argv[1]);
    if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

    std::ifstream invect(argv[2]);
    if (!input) { cerr << "Error opening vector file " << argv[2] << endl; return -1; }

    
    typedef Modular<double> Field;
    double q = atof(argv[3]);
    Field F(q);
    MatrixStream< Field > ms ( F, input );
    BlasBlackbox<Field> A (ms); // A.write(std::cout);
    cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
    
    std::vector<Field::Element> X( A.coldim()),B(A.rowdim());
	
    for(std::vector<Field::Element>::iterator it=B.begin();
	    it != B.end(); ++it)
	    invect >> *it;
    
    std::cout << "B is [ "<<B<< "]" << std::endl;
    
    solve (X, A, B, Method::BlasElimination());
    
    std::cout << "(BlasElimination) Solution is [ "<<X<< "]" << std::endl;		
    std::vector<Field::Element> r(A.rowdim());
    BlasMatrixDomain<Field> BMD(F);
    BMD.mul(r, static_cast<BlasMatrix<Field::Element>& >(A), X);
    //A.apply (r,X);
    VectorDomain<Field> VD(F);
    if (VD.areEqual (r,B))
	    std::cout<<"CHECK"<<endl;
    else{
	    std::cout<<"FAIL"<<endl;
	    cout<<"r = "<<r<<endl;
	    cout<<"B = "<<B<<endl;
    }

    return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
