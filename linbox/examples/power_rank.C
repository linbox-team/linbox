/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/** \file examples/power-rank.C
\brief Rank of sparse matrix over Z or Zp.
\ingroup examples
*/
#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/field/givaro-zpz.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/algorithms/smith-form-sparseelim-local.h"

using namespace LinBox;
using namespace std;



int main (int argc, char **argv)
{
    commentator.setMaxDetailLevel (-1);
    commentator.setMaxDepth (-1);
    commentator.setReportStream (std::cerr);

	if (argc < 4 || argc > 4) 
	{	cerr << "Usage: rank <matrix-file-in-supported-format> <prime> <prime-power>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	long unsigned int r;

	if (argc == 4) { 
            LinBox::int64 p = atoi(argv[2]);
            LinBox::int64 q = atoi(argv[3]);
                typedef GivaroZpz<Std64> Field;
                Field F(q);
                MatrixStream<Field> ms( F, input );
                SparseMatrix<Field, Vector<Field>::SparseSeq > B (ms);
                cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
                if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;

                    // using Sparse Elimination
                PowerGaussDomain< Field > PGD( F );
                std::vector<std::pair<size_t,size_t> > local;

		Timer tq; tq.clear(); tq.start();
                PGD(local, B, q, p);    
		tq.stop();

                
                std::cout << "Local Smith Form : ("; 
                for (std::vector<std::pair<size_t,size_t> >::const_iterator  p = local.begin(); 
                     p != local.end(); ++p) 
                    std::cout << p->first << " " << p->second << ", "; 
                cout << ")" << endl; 


		std::cerr << tq << std::endl;
 	}

	return 0;
}
