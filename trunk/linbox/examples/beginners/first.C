/* File:	tests/Algorithms/wiedemann_linsolve1.h
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests on Wiedemann algorithm for solving nonhomogeneous linear
 * equations
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <utility>

#include "linbox/field/large-modular.h"
//#include "linbox/blackbox/sparse0.h"
#include "linbox/blackbox/sparse1.h" 

using namespace LinBox;

int main(int argc, char* argv[])
{

        char* in_file = argv[1];
        char* out_file = argv[2];

        typedef LargeModular  Field;
	typedef Field::element Element;
	typedef Field::RandIter RandIter;
	typedef std::list< pair<size_t, Element> > Row;
	typedef std::vector<Element> Vector;

        Field K(7);

        ofstream out_stream(out_file);
        ifstream in_stream(in_file);

//	SparseMatrix0<Field, Row, Vector>  A(K,4,4); 
//      A.read(in_stream);
//	A.write(out_stream);

        SparseBlackBoxDom< Field > A(K) ;

        A.read(in_stream);
        A.write(out_stream);

}
