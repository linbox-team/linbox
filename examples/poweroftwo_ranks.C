/* examples/poweroftwo_ranks.C
 *
 * Copyright (C) 2012 LinBox
 * Written by J-G Dumas
 * Time-stamp: <21 Dec 18 10:08:04 Jean-Guillaume.Dumas@imag.fr>
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

/** \file examples/poweroftwo_ranks.C
 * @example  examples/poweroftwo_ranks.C
 \brief Ranks of sparse matrix modulo 2^k
 \ingroup examples
*/
#include <linbox/linbox-config.h>

#include <iostream>

#include <givaro/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/smith-form.h>
#include <linbox/algorithms/smith-form-sparseelim-poweroftwo.h>

using namespace LinBox;
using namespace std;

template<class Int_type, class Ring_type = Givaro::ZRing<Int_type> >
void runpoweroftworank(ifstream& input, const size_t exponent, size_t StPr) {
//     typedef std::vector<std::pair<size_t,Int_type> > Smith_t;
    typedef Ring_type Ring; // signed ?
    typedef LinBox::SparseMatrix<Ring, 
        LinBox::SparseMatrixFormat::SparseSeq > SparseMat;

    SmithList<Ring> local;
    Ring R;
    LinBox::MatrixStream<Ring> ms( R, input );
    SparseMat A (ms);

    input.close();
    LinBox::PowerGaussDomainPowerOfTwo< Int_type > PGD;
    LinBox::GF2 F2;
    Permutation<GF2> Q(F2,A.coldim());
            
    cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
    if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(cout,Tag::FileFormat::Maple) << endl;

    Givaro::Timer tim; 
    tim.clear(); tim.start();
    if (StPr)
        PGD(local, A, Q, exponent, StPr);
    else
        PGD(local, A, Q, exponent);
    tim.stop();

    R.write(std::cout << "Local Smith Form ") << " : " << std::endl << '(';
	int num = A.rowdim();
    for (auto  p = local.begin(); p != local.end(); ++p) {
        std::cout << '[' << p->first << ',' << p->second << "] ";
		num -= p->second;
	}
	if (num > 0) std::cout << '[' << F2.zero << ',' << num << "] ";
	std::cout << ')' << std::endl;


//         // Reposition Output with empty rows at the end
//     auto newend = std::remove_if(
//         A.rowBegin(), A.rowEnd(),
//         [](typename SparseMat::ConstRow V)->bool { return V.size()==0; });
//     A.refRep().erase(newend, A.rowEnd());
//     A.refRep().resize(A.rowdim());

    if (A.rowdim() <= 20 && A.coldim() <= 20) {
        A.write(cerr,Tag::FileFormat::Maple) << endl;
        Q.write(cerr,Tag::FileFormat::Maple) << endl;
    }

    std::cerr << tim << std::endl;
}

int main (int argc, char **argv) {
    commentator().setMaxDetailLevel (-1);
    commentator().setMaxDepth (-1);
    commentator().setReportStream (std::cerr);

    if (argc < 3 || argc > 5)
    {	cerr << "Usage: rank <matrix-file-in-supported-format> <power of two exponent> [<method>] [<flag>]" << endl; return -1; }

    ifstream input (argv[1]);
    if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }
    int method( argc>3? atoi(argv[3]): 0 );
    
    Givaro::Timer tim;
    int exponent = atoi(argv[2]);

        // 1: StPr |= PRESERVE_UPPER_MATRIX
        // 2: StPr |= PRIVILEGIATE_REDUCING_FILLIN
        // 4: StPr |= PRIVILEGIATE_NO_COLUMN_PIVOTING
    size_t StPr( argc>4? atoi(argv[4]): 0);

    if ((method == 2) || ((method == 0) && (exponent >= (1<<11))) ) {
        runpoweroftworank<Givaro::Integer>(input, exponent, StPr);
    } else {
        if ((method == 1) || ((method == 0) && (exponent < 64)) ) {
            runpoweroftworank<uint64_t>(input, exponent, StPr);
        } else {
            switch (method) {
                case 6: runpoweroftworank<RecInt::ruint<6>>(input, exponent, StPr); break;
                case 7: runpoweroftworank<RecInt::ruint<7>>(input, exponent, StPr); break;
                case 8: runpoweroftworank<RecInt::ruint<8>>(input, exponent, StPr); break;
                case 9: runpoweroftworank<RecInt::ruint<9>>(input, exponent, StPr); break;
                case 10: runpoweroftworank<RecInt::ruint<10>>(input, exponent, StPr); break;
                case 11: runpoweroftworank<RecInt::ruint<11>>(input, exponent, StPr); break;
                default: std::cerr << "Choose between ruint<6> ... ruint<11>" << std::endl;
                    break;
            }
        }
    }

    std::cerr << tim << std::endl;
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
