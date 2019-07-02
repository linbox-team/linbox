/*
 * examples/power_rank.C
 *
 * Copyright (C) 2005, 2010 J-G Dumas
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

/** \file examples/power_rank.C
 * @example  examples/power_rank.C
 \brief Rank of sparse matrix over Z or Zp.
 \ingroup examples
*/
#include <linbox/linbox-config.h>

#include <iostream>

#include <givaro/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/algorithms/smith-form-sparseelim-local.h>

using namespace LinBox;
using namespace std;


template<typename Base>
int tmain (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::cerr);

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

    Base p; Givaro::Caster(p,Givaro::Integer(argv[2]));
    Base q; Givaro::Caster(q,Givaro::Integer(argv[3]));
    typedef Givaro::Modular<Base> Field;
    typedef SparseMatrix<Field, SparseMatrixFormat::SparseSeq > SparseMat;
    Field F(q);
    MatrixStream<Field> ms( F, input );
    SparseMat B (ms);
    cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
    if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout,Tag::FileFormat::Maple) << endl;

		// using Sparse Elimination
    PowerGaussDomain< Field > PGD( F );
    std::vector<std::pair<Base,size_t> > local;
    Permutation<Field> Q(F,B.coldim());

        // 1: StPr |= PRESERVE_UPPER_MATRIX
        // 2: StPr |= PRIVILEGIATE_REDUCING_FILLIN
        // 4: StPr |= PRIVILEGIATE_NO_COLUMN_PIVOTING
    size_t StPr( argc>5? atoi(argv[5]): 0);

    Givaro::Timer tq; tq.clear(); tq.start();
    if (StPr)
        PGD(local, B, Q, q, p, StPr);
    else
        PGD(local, B, Q, q, p);
    tq.stop();


    F.write(std::cout << "Local Smith Form ") << " : " << std::endl << '(';
	int num = B.rowdim();
    for (auto ip = local.begin(); ip != local.end(); ++ip) {
        std::cout << '[' << ip->first << ',' << ip->second << "] ";
		num -= ip->second;
	}
	if (num > 0) std::cout << '[' << F.zero << ',' << num << "] ";
	std::cout << ")" << std::endl;

        // Reposition Output with empty rows at the end
    auto newend = std::remove_if(
        B.rowBegin(), B.rowEnd(),
        [](typename SparseMat::ConstRow V)->bool { return V.size()==0; });
    B.refRep().erase(newend, B.rowEnd());
    B.refRep().resize(B.rowdim());

    if (B.rowdim() <= 20 && B.coldim() <= 20) {
        B.write(cerr,Tag::FileFormat::Maple) << endl;
        Q.write(cerr,Tag::FileFormat::Maple) << endl;
    }

    std::cerr << tq << std::endl;

	return 0;
}

int main(int argc, char ** argv) {
	if (argc < 4 || argc > 6) {
        cerr << "Usage: power_rank <matrix-file-in-supported-format> <prime> <prime-power> [<method>] [<flag>]" << endl;
        cerr << "       methods: \
						0=automatic, \
						1=int_64_t, \
						2=Integer, \
						6-11=ruint" << endl;
        return -1; }

    size_t method( argc>4? atoi(argv[4]) : 0);

    Givaro::Integer q(argv[3]);
    const size_t logq( (size_t)ceil(logtwo(q)) );

    if ( (method == 1) || ( (method==0) && (logq<63) ) ) {
        return tmain<int64_t>(argc,argv);
    } else {
        if ( (method == 2) ) {
            return tmain<Givaro::Integer>(argc,argv);
        } else {
            if (! method) method = (size_t)std::max(6.,  ceil(log(logq)/log(2.))  );
            switch (method) {
                case 6: return tmain<RecInt::ruint<6>>(argc,argv);
                case 7: return tmain<RecInt::ruint<7>>(argc,argv); 
                case 8: return tmain<RecInt::ruint<8>>(argc,argv);
                case 9: return tmain<RecInt::ruint<9>>(argc,argv);
                case 10: return tmain<RecInt::ruint<10>>(argc,argv);
                case 11: return tmain<RecInt::ruint<11>>(argc,argv);
                default: return tmain<Givaro::Integer>(argc,argv);
            }
        }
    }
    
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
