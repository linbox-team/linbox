/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
#include "linbox-config.h"

#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/field/gf2.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/rank.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
    commentator.setReportStream (std::cerr);

	if (argc < 2 || argc > 3) 
	{	cerr << "Usage: rank <matrix-file-in-SMS-format> [<p>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	long unsigned int r;

	if (argc == 2) { // rank over the integers.

	   /* We could pick a random prime and work mod that prime, But the point here 
	   is that the rank function in solutions/ handles that issue.  Our matrix here 
	   is an integer matrix and our concept is that we are getting the rank of that 
	   matrix by some blackbox magic inside linbox.
	   */
		GMP_Integers ZZ;
		SparseMatrix<GMP_Integers> A (ZZ);
		A.read (input);
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		rank (r, A);

		cout << "Rank is " << r << endl;
	}
	if (argc == 3) { 
		double q = atof(argv[2]);
		if (q == 2.0) {
			typedef LinBox::GF2 Field;
                        Field F;
                        SparseMatrix<Field, Vector<Field>::SparseSeq > B (F);
                        B.read (input);
                        cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
                        
                        GaussDomain<Field> GD ( F );
                        GD.rankinLinearPivoting (r, B, B.rowdim(), B.coldim());
                        
                        cout << "Rank mod 2 is " << r << endl;
                } else {
                    typedef Modular<double> Field;
                    Field F(q);
                    SparseMatrix<Field, Vector<Field>::SparseSeq > B (F);
                    B.read (input);
                    cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
                    
                    rankin (r, B);
                    
                    cout << "Rank is " << r << endl;
                }
                
	}

	return 0;
}
