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

#include <linbox/linbox-config.h>

#include <iostream>
#include <givaro/givintfactor.h>
#include <fflas-ffpack/paladin/parallel.h>
#include <linbox/algorithms/smith-form-valence.h>
#include <vector>

using namespace LinBox;

template<class Blackbox>
std::vector<Givaro::Integer>& smithValence(std::vector<Givaro::Integer>& SmithDiagonal, 
                                           Givaro::Integer& valence,
                                           const Blackbox& A,
                                           const std::string& filename,
                                           Givaro::Integer& coprimeV, 
                                           size_t method=0) {
        // method for valence squarization: 0 for automatic, 1 for aat, 2 for ata
        // Blackbox provides the Integer matrix rereadable from filename
        // if valence != 0, then the valence is not computed and the parameter is used
        // if coprimeV != 0, then this value is supposed to be coprime with the valence

    if (valence == 0) {
        squarizeValence(valence, A, method);
    }

	std::clog << "Valence is " << valence << std::endl;

	std::vector<Givaro::Integer> Moduli;
	std::vector<size_t> exponents;

    std::clog << "Some factors (50000 factoring loop bound): ";
	
    FTD.set(Moduli, exponents, valence, 50000);

	auto eit=exponents.begin();
	for(auto mit: Moduli) std::clog << mit << '^' << *eit << ' ';
	std::clog << std::endl;
    smith.resize(Moduli.size());

	Givaro::IntFactorDom<> FTD;
    if (coprimeV == 1) {
        coprimeV=2;
        while ( gcd(valence,coprimeV) > 1 ) {
            FTD.nextprimein(coprimeV);
        }
    }

	std::clog << "num procs: " << omp_get_num_procs() << std::endl;
    PAR_BLOCK { 
        std::clog << "cur threads: " << NUM_THREADS << std::endl;
        std::clog << "max threads: " << MAX_THREADS << std::endl;
    }

	std::vector< PairIntRk > smith;
    std::vector<std::vector<size_t> > AllRanks(Moduli.size());
    size_t coprimeR;


    PAR_BLOCK { 
       SYNCH_GROUP(
            for(size_t j=0; j<Moduli.size(); ++j) {
                smith[j].first = Moduli[j];
                { TASK(MODE(READ(smith[j].first) WRITE(smith[j].second) VALUE(j)),
                {
                    LRank(smith[j].second, filename.c_str(), smith[j].first);
                })}
            }
        	{ TASK(MODE(READ(coprimeV) WRITE(coprimeR)),
            {
                LRank(coprimeR, filename.c_str(), coprimeV);
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
                        AllPowersRanks(AllRanks[j], smithj, exponentsj, coprimeR, filename.c_str());
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


    return SmithDiagonal;
}

    
    




int main (int argc, char **argv)
{

	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: smithvalence <matrix-file-in-supported-format> [-ata|-aat|valence] [coprime]" << std::endl;
        std::cerr << "       Optional parameters valence and coprime are integers." << std::endl;
        std::cerr << "       Prime factors of valence will be used for local computation." << std::endl;
        std::cerr << "       coprime will be used for overall Z rank computation." << std::endl;
		return -1;
	}

    const std::string filename(argv[1]);

	std::ifstream input (filename);
	if (!input) { std::cerr << "Error opening matrix file " << filename << std::endl; return -1; }

	Givaro::ZRing<Integer> ZZ;
	MatrixStream< Givaro::ZRing<Integer> > ms( ZZ, input );
	typedef SparseMatrix<Givaro::ZRing<Integer>>  Blackbox;
	Blackbox A (ms);
	input.close();

	std::clog << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

    size_t method=0;
	Givaro::Integer val_A(0);

	if (argc >= 3) {
		if (strcmp(argv[2],"-ata") == 0) {
            method=2;
		}
		else if (strcmp(argv[2],"-aat") == 0) {
            method=1;
		}
		else {
			std::clog << "Suppose primes are contained in " << argv[2] << std::endl;
			val_A = LinBox::Integer(argv[2]);
		}
	}

    Givaro::Integer coprimeV=1;
	if (argc >= 4) {
		std::clog << "Suppose " << argv[3] << " is coprime with Smith form" << std::endl;
		coprimeV =Givaro::Integer(argv[3]);
	}

    std::vector<Givaro::Integer> SmithDiagonal;

	LinBox::Timer chrono; chrono.start();
    smithValence(SmithDiagonal, val_A, A, filename, coprimeV, method);
	chrono.stop();

	std::clog << "Integer Smith Form :" << std::endl;
    compressedSmith(std::cout, SmithDiagonal, A.rowdim(), A.coldim()) << std::endl;

	std::clog << chrono << std::endl;
    

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
