/*
 * examples/smithvalence.C
 *
 * Copyright (C) 2010  J-G Dumas
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
//#include "linbox-config.h"
#include <iostream>

#include "smithvalence.h"

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

	PID_integer ZZ;
	MatrixStream< PID_integer > ms( ZZ, input );
	typedef SparseMatrix<PID_integer>  Blackbox;
	Blackbox A (ms);
	input.close();

	std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

	PID_integer::Element val_A;

 Givaro::Timer chrono; chrono.start();
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
			val_A = integer(argv[2]);
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

	std::vector<integer> Moduli;
	std::vector<size_t> exponents;
 Givaro::IntFactorDom<> FTD;

	typedef std::pair<integer,unsigned long> PairIntRk;
	std::vector< PairIntRk > smith;


	integer coprimeV=2;
	if (argc >= 4) {
		coprimeV = integer(argv[3]);
	}
	while ( gcd(val_A,coprimeV) > 1 ) {
		FTD.nextprimein(coprimeV);
	}

	if (argc >= 4) {
		std::cout << "Suppose " << argv[3] << " is coprime with Smith form" << std::endl;
	}

	std::cout << "integer rank: " << std::endl;

	unsigned long coprimeR; LRank(coprimeR, argv[1], coprimeV);
	smith.push_back(PairIntRk(coprimeV, coprimeR));
	//         std::cerr << "Rank mod " << coprimeV << " is " << coprimeR << std::endl;

	std::cout << "Some factors (5000 factoring loop bound): ";
	FTD.set(Moduli, exponents, val_A, 5000);
	std::vector<size_t>::const_iterator eit=exponents.begin();
	for(std::vector<integer>::const_iterator mit=Moduli.begin();
	    mit != Moduli.end(); ++mit,++eit)
		std::cout << *mit << '^' << *eit << ' ';
	std::cout << std::endl;

	std::vector<integer> SmithDiagonal(coprimeR,integer(1));


	for(std::vector<integer>::const_iterator mit=Moduli.begin();
	    mit != Moduli.end(); ++mit) {
		unsigned long r; LRank(r, argv[1], *mit);
		//             std::cerr << "Rank mod " << *mit << " is " << r << std::endl;
		smith.push_back(PairIntRk(*mit, r));
		for(size_t i=r; i < coprimeR; ++i)
			SmithDiagonal[i] *= *mit;
	}


	eit=exponents.begin();
	std::vector<PairIntRk>::const_iterator sit=smith.begin();
	for( ++sit; sit != smith.end(); ++sit, ++eit) {
		if (sit->second != coprimeR) {
			std::vector<size_t> ranks;
			ranks.push_back(sit->second);
            size_t effexp;
			if (*eit > 1) {
				PRank(ranks, effexp, argv[1], sit->first, *eit, coprimeR);
			}
			else {
				PRank(ranks, effexp, argv[1], sit->first, 2, coprimeR);
			}
			if (ranks.size() == 1) ranks.push_back(coprimeR);

            if (effexp < *eit) {
                for(size_t expo = effexp<<1; ranks.back() < coprimeR; expo<<=1) {
                    PRankInteger(ranks, argv[1], sit->first, expo, coprimeR);
                }
            } else {

                for(size_t expo = (*eit)<<1; ranks.back() < coprimeR; expo<<=1) {
                    PRank(ranks, effexp, argv[1], sit->first, expo, coprimeR);
                    if (ranks.size() < expo) {
                        std::cerr << "It seems we need a larger prime power, it will take longer ..." << std::endl;
                            // break;
                        PRankInteger(ranks, argv[1], sit->first, expo, coprimeR);
                    }
                }
            }

			std::vector<size_t>::const_iterator rit=ranks.begin();
			// unsigned long modrank = *rit;
			for(++rit; rit!= ranks.end(); ++rit) {
				if ((*rit)>= coprimeR) break;
				for(size_t i=(*rit); i < coprimeR; ++i)
					SmithDiagonal[i] *= sit->first;
				// modrank = *rit;
			}
		}
	}

	integer si=1;
	size_t num=0;
	for( std::vector<integer>::const_iterator dit=SmithDiagonal.begin();
	     dit != SmithDiagonal.end(); ++dit) {
		if (*dit == si) ++num;
		else {
			if (num > 0)
				std::cerr << '[' << si << ',' << num << "] ";
			num=1;
			si = *dit;
		}
	}
	std::cerr << '[' << si << ',' << num << "] " << std::endl;
	chrono.stop();
	std::cerr << chrono << std::endl;


	return 0;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

