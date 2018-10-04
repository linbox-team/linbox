/*
 * examples/omp_smithvalence.C
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

/**\file examples/omp_smithvalence.C
 * @example  examples/omp_smithvalence.C
 \brief Valence of sparse matrix over Z or Zp.
 \ingroup examples
*/
#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif

#include <iostream>
#include <omp.h>
#include "smithvalence.h"
using namespace LinBox;

int main (int argc, char **argv)
{
        //     commentator().setMaxDetailLevel (-1);
        //     commentator().setMaxDepth (-1);
        //     commentator().setReportStream (std::cerr);


	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: omp_smithvalence <matrix-file-in-supported-format> [-ata|-aat|valence] [coprime]" << std::endl;
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

	typedef std::pair<Givaro::Integer,unsigned long> PairIntRk;
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
	std::vector<size_t>::const_iterator eit=exponents.begin();
	for(std::vector<Givaro::Integer>::const_iterator mit=Moduli.begin();
	    mit != Moduli.end(); ++mit,++eit)
		std::cout << *mit << '^' << *eit << ' ';
	std::cout << std::endl;
    smith.resize(Moduli.size());

    unsigned long coprimeR; 
 
	std::cout << "num procs: " << omp_get_num_procs() << std::endl;
	std::cout << "max threads: " << omp_get_max_threads() << std::endl;

// #pragma omp parallel for shared(smith, Moduli) 
// 	for(size_t j=0; j<(Moduli.size()+1); ++j) {
//         if (j >= Moduli.size()) {
//             LRank(coprimeR, argv[1], coprimeV);
//             std::cerr << "Integer Rank (mod " << coprimeV << ") is " << coprimeR << " on thread: " << omp_get_thread_num() << std::endl;
//         } else {            
//             unsigned long r; LRank(r, argv[1], Moduli[j]);
//             std::cerr << "Rank mod " << Moduli[j] << " is " << r << " on thread: " << omp_get_thread_num() << std::endl;
//             smith[j] = PairIntRk( Moduli[j], r);
//         }
// 	}
#pragma omp parallel 
    {
#pragma omp single
        {
            for(size_t j=0; j<Moduli.size(); ++j) {
                smith[j].first = Moduli[j];
#pragma omp task shared(smith,argv) firstprivate(j)
                LRank(smith[j].second, argv[1], smith[j].first);
            }
            LRank(coprimeR, argv[1], coprimeV);
        }
	}

 
	std::vector<Givaro::Integer> SmithDiagonal(coprimeR,Givaro::Integer(1));

        /*
          for(std::vector<Givaro::Integer>::const_iterator mit=Moduli.begin();
          mit != Moduli.end(); ++mit) {
          unsigned long r; LRank(r, argv[1], *mit);
          std::cerr << "Rank mod " << *mit << " is " << r << std::endl;
          smith.push_back(PairIntRk(*mit, r));
          for(size_t i=r; i < coprimeR; ++i)
          SmithDiagonal[i] *= *mit;
          }
        */
#pragma omp parallel for shared(SmithDiagonal, smith, exponents)
	for(size_t j=0; j<Moduli.size(); ++j) {
        
        if (smith[j].second != coprimeR) {
            std::vector<size_t> ranks;
            ranks.push_back(smith[j].second);
            size_t effexp;
            if (exponents[j] > 1) {
                if (smith[j].first == 2)
                    PRankPowerOfTwo(ranks, effexp, argv[1], exponents[j], coprimeR);
                else
                    PRank(ranks, effexp, argv[1], smith[j].first, exponents[j], coprimeR);
            }
            else {
                if (smith[j].first == 2)
                    PRankPowerOfTwo(ranks, effexp, argv[1], 2, coprimeR);
                else
                    PRank(ranks, effexp, argv[1], smith[j].first, 2, coprimeR);
            }
            if (ranks.size() == 1) ranks.push_back(coprimeR);

            if (effexp < exponents[j]) {
                for(size_t expo = effexp<<1; ranks.back() < coprimeR; expo<<=1) {
                    if (smith[j].first == 2)
                        PRankIntegerPowerOfTwo(ranks, argv[1], expo, coprimeR);
                    else
                        PRankInteger(ranks, argv[1], smith[j].first, expo, coprimeR);
                }
            } else {

                for(size_t expo = (exponents[j])<<1; ranks.back() < coprimeR; expo<<=1) {
                    if (smith[j].first == 2)
                        PRankPowerOfTwo(ranks, effexp, argv[1], expo, coprimeR);
                    else
                        PRank(ranks, effexp, argv[1], smith[j].first, expo, coprimeR);
                    if (ranks.size() < expo) {
                        std::cerr << "It seems we need a larger prime power, it will take longer ..." << std::endl;
                            // break;
                        if (smith[j].first == 2)
                            PRankIntegerPowerOfTwo(ranks, argv[1], expo, coprimeR);
                        else
                            PRankInteger(ranks, argv[1], smith[j].first, expo, coprimeR);
                    }
                }
            }

#pragma omp critical
            {   
                for(size_t i=smith[j].second; i < coprimeR; ++i) {
                    SmithDiagonal[i] *= smith[j].first;
                }
                
                std::vector<size_t>::const_iterator rit=ranks.begin();
// 			unsigned long modrank = *rit;
                for(++rit; rit!= ranks.end(); ++rit) {
                    if ((*rit)>= coprimeR) break;
                    for(size_t i=(*rit); i < coprimeR; ++i)
                        SmithDiagonal[i] *= smith[j].first;
// 				modrank = *rit;
                }
            }
        }
	}
	chrono.stop();

    Givaro::Integer si=1;
	size_t num=0;
	std::cerr << "Integer Smith Form :" << std::endl;
    std::cout << '(';
	for( std::vector<Givaro::Integer>::const_iterator dit=SmithDiagonal.begin();
	     dit != SmithDiagonal.end(); ++dit) {
		if (*dit == si) ++num;
		else {
			std::cerr << '[' << si << ',' << num << "] ";
			num=1;
			si = *dit;
		}
	}
	std::cerr << '[' << si << ',' << num << "] " << std::endl;
	num = std::min(A.rowdim(),A.coldim()) - SmithDiagonal.size();
	si = ZZ.zero;
	if (num > 0) std::cout << '[' << si << ',' << num << ']';
	std::cout << ')' << std::endl;
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
