/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/** \file examples/sparseelimrank.C
 * Time-stamp: <26 Mar 07 11:46:18 Jean-Guillaume.Dumas@imag.fr>
\brief Gaussian elimination Rank of sparse matrix over Z or Zp.
\ingroup examples
*/
//#include "linbox-config.h"

#include <iostream>
#include <vector>
#include <utility>

template<class T>
std::ostream& operator<< (std::ostream& o, const std::vector<std::pair<size_t, T> >& C) {
          for(typename std::vector<std::pair<size_t, T> >::const_iterator refs =  C.begin();
                                refs != C.end() ;
                                      ++refs )
		  o << '(' << refs->first << ';' << refs->second << ')';
            return o << std::endl;
}


#include "linbox/field/modular-double.h"
#include "linbox/field/gf2.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/solutions/rank.h"
#include "linbox/util/matrix-stream.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
    commentator.setMaxDetailLevel (-1);
    commentator.setMaxDepth (-1);
    commentator.setReportStream (std::cerr);

	if (argc < 2 || argc > 3) 
	{	cerr << "Usage: rank <matrix-file-in-supported-format> [<p>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	long unsigned int r;

	if (argc == 2) { // rank over the integers.

	   /* We could pick a random prime and work mod that prime, But the point here 
	   is that the rank function in solutions/ handles that issue.  Our matrix here 
	   is an integer matrix and our concept is that we are getting the rank of that 
	   matrix by some blackbox magic inside linbox.
	   */
		PID_integer ZZ;
		MatrixStream<PID_integer> ms( ZZ, input );
		SparseMatrix<PID_integer> A ( ms );
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		rank (r, A, Method::SparseElimination() );

		cout << "Rank is " << r << endl;
	}
	if (argc == 3) { 
		double q = atof(argv[2]);
                    typedef Modular<double> Field;
                    Field F(q);
		    MatrixStream<Field> ms( F, input );
                    SparseMatrix<Field, Vector<Field>::SparseSeq > B (ms);
                    cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
                    if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;


		    Method::SparseElimination SE;
		    SE.strategy(Specifier::PIVOT_NONE);
			// using Sparse Elimination
                    rank (r, B, SE);                    
                    if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;
                    cout << "Rank is " << r << endl;

		    SE.strategy(Specifier::PIVOT_LINEAR);
			// using Sparse Elimination
                    rank (r, B, SE);                    
                    if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;
                    cout << "Rank is " << r << endl;


	}

	return 0;
}
