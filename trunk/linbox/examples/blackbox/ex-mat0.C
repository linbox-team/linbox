/** @name examples/blackbox/ex-mat0.C
 * @author William J. Turner for the LinBox group
 *
 * @memo usage: ex-mat0 in-file out-file 
 *
 * @doc
 * Run tests on Wiedemann algorithm for solving nonhomogeneous linear
 * equations
 *
 * FIXME What does it do?  I think this may be a remnant, has evolved
 * into one of the other examples.  delete it?
 */
//@{

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <utility>

#include "linbox/field/modular.h"
//#include "linbox/blackbox/sparse0.h"
#include "linbox/blackbox/sparse1.h" 

using namespace LinBox;
using namespace std;

int main(int argc, char* argv[])
{

	if (argc != 3)
	{	cerr << "usage: " << argv[0] << " in_file out_file" << endl;
		return -1;
	}
        char* in_file = argv[1];
        char* out_file = argv[2];

        typedef Modular<uint32>  Field;
	typedef Field::Element Element;
	typedef Field::RandIter RandIter;
	typedef std::list< pair<size_t, Element> > Row;
	typedef std::vector<Element> Vector;

        Field K(7);

        ofstream out_stream(out_file);
        ifstream in_stream(in_file);

//	SparseMatrix<Field, Row, Vector>  A(K,4,4); 
//      A.read(in_stream);
//	A.write(out_stream);

        SparseBlackBoxDom< Field > A(K) ;

        A.read(in_stream);
        A.write(out_stream);

}
//@}
