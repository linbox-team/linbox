/*
 * examples/smithvalence.C
 * Copyright (c) Linbox
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

/**\file examples/smithvalence.C
 * @example  examples/smithvalence.C
 \brief Valence of sparse matrix over Z or Zp.
 \ingroup examples
*/
#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif
#include <cassert>
#include <linbox/linbox-config.h>

#include <iostream>
#include <fflas-ffpack/paladin/parallel.h>
#include <linbox/algorithms/smith-form-valence.h>
#include <vector>

using namespace LinBox;


int main (int argc, char **argv)
{
        //     commentator().setMaxDetailLevel (-1);
        //     commentator().setMaxDepth (-1);
        //     commentator().setReportStream (std::cerr);


	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: smithvalence <matrix-file-in-supported-format> [-ata|-aat|valence] [coprime]" << std::endl;
        std::cerr << "       Optional parameters valence and coprime are integers." << std::endl;
        std::cerr << "       Prime factors of valence will be used for local computation." << std::endl;
        std::cerr << "       coprime will be used for overall rank computation." << std::endl;
		return -1;
	}

	std::ifstream input (argv[1]);
	if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }

	Givaro::ZRing<Integer> ZZ;
	MatrixStream< Givaro::ZRing<Integer> > ms( ZZ, input );
	typedef SparseMatrix<Givaro::ZRing<Integer>>  Blackbox;
	Blackbox A (ms);
	input.close();

	std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

	Givaro::ZRing<Integer>::Element val_A;

	LinBox::Timer chrono; chrono.start();
	if (argc >= 3) {
		Transpose<Blackbox> T(&A);
		if (strcmp(argv[2],"-ata") == 0) {
			Compose< Transpose<Blackbox>, Blackbox > C (&T, &A);
			std::cout << "A^T A is " << C.rowdim() << " by " << C.coldim() << std::endl;
			valence(val_A, C);
		}
		else if (strcmp(argv[2],"-aat") == 0) {
			Compose< Blackbox, Transpose<Blackbox> > C (&A, &T);
			std::cout << "A A^T is " << C.rowdim() << " by " << C.coldim() << std::endl;
			valence(val_A, C);
		}
		else {
			std::cout << "Suppose primes are contained in " << argv[2] << std::endl;
			val_A = LinBox::Integer(argv[2]);
		}
	}
	else {
		if (A.rowdim() != A.coldim()) {
			std::cerr << "Valence works only on square matrices, try either to change the dimension in the matrix file, or to compute the valence of A A^T or A^T A, via the -aat or -ata options."  << std::endl;
			exit(0);
		}
		else
			valence (val_A, A);
	}

	std::cout << "Valence is " << val_A << std::endl;

	std::vector<Givaro::Integer> Moduli;
	std::vector<size_t> exponents;
	Givaro::IntFactorDom<> FTD;

	std::vector< PairIntRk > smith;


    Givaro::Integer coprimeV=2;
	if (argc >= 4) {
		coprimeV =Givaro::Integer(argv[3]);
	}
	while ( gcd(val_A,coprimeV) > 1 ) {
		FTD.nextprimein(coprimeV);
	}

	if (argc >= 4) {
		std::cout << "Suppose " << argv[3] << " is coprime with Smith form" << std::endl;
	}

    std::cout << "Some factors (50000 factoring loop bound): ";
	
    FTD.set(Moduli, exponents, val_A, 50000);

	auto eit=exponents.begin();
	for(auto mit: Moduli) std::cout << mit << '^' << *eit << ' ';
	std::cout << std::endl;
    smith.resize(Moduli.size());

	std::cout << "num procs: " << omp_get_num_procs() << std::endl;
    PAR_BLOCK { 
        std::cout << "cur threads: " << NUM_THREADS << std::endl;
        std::cout << "max threads: " << MAX_THREADS << std::endl;
    }

        size_t coprimeR;

        std::vector<Givaro::Integer> SmithDiagonal(coprimeR,Givaro::Integer(1));
        std::vector<std::vector<size_t> > AllRanks(Moduli.size());


    PAR_BLOCK { 
       SYNCH_GROUP(
            for(size_t j=0; j<Moduli.size(); ++j) {
                smith[j].first = Moduli[j];
                { TASK(MODE(READ(smith[j].first) WRITE(smith[j].second) VALUE(j)),
                {
                    LRank(smith[j].second, argv[1], smith[j].first);
                })}
            }
        	{ TASK(MODE(READ(coprimeV) WRITE(coprimeR)),
            {
                LRank(coprimeR, argv[1], coprimeV);
            })}
        )


        SmithDiagonal.resize(coprimeR,Givaro::Integer(1));



        SYNCH_GROUP(
            for(size_t j=0; j<Moduli.size(); ++j) {
                if (smith[j].second != coprimeR) {
                    const PairIntRk& smithj(smith[j]);
                    const size_t& exponentsj(exponents[j]);
                    { TASK(MODE(CONSTREFERENCE(smithj, exponentsj, coprimeR) WRITE(AllRanks[j])),
                    {
                        AllPowersRanks(AllRanks[j], smithj, exponentsj, coprimeR, argv[1]);
                    })}
                }
            }
        )
    }
    
    for(size_t j=0; j<Moduli.size(); ++j) {
        if (smith[j].second != coprimeR) {
            populateSmithForm(SmithDiagonal, AllRanks[j], smith[j], coprimeR);
        }
    }    
    
	chrono.stop();

	std::cerr << "Integer Smith Form :" << std::endl;
    compressedSmith(std::cout, SmithDiagonal, A.rowdim(), A.coldim()) << std::endl;

	std::cerr << chrono << std::endl;
    

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
