/** 
 * examples/smithvalence.C
 *
 * Copyright (C) 2010  J-G Dumas
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */

/**\file examples/smithvalence.C
\brief Valence of sparse matrix over Z or Zp.
\ingroup examples
*/
//#include "linbox-config.h"
#include <iostream>

#include "linbox/field/gf2.h"
#include "linbox/field/modular-double.h"
#include "linbox/field/givaro-zpz.h"
#include "linbox/field/field-traits.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/valence.h"
#include "linbox/algorithms/smith-form-sparseelim-local.h"
#include "linbox/util/matrix-stream.h"

#include <givaro/givintnumtheo.h>

template<class Field>
unsigned long& TempLRank(unsigned long& r, char * filename, const Field& F) {
    std::ifstream input(filename);
    LinBox::MatrixStream< Field > msf( F, input );
    LinBox::SparseMatrix<Field,typename LinBox::Vector<Field>::SparseSeq> FA(msf);
    input.close();
    LinBox::Timer tim; tim.start();
    LinBox::rankin(r, FA);
    tim.stop();
    F.write(std::cout << "Rank over ") << " is " << r << ' ' << tim << std::endl;
    return r;
}
    
unsigned long& TempLRank(unsigned long& r, char * filename, const LinBox::GF2& F2) {
    std::ifstream input(filename);
    LinBox::ZeroOne<LinBox::GF2> A;
    A.read(input);
    input.close();
    LinBox::Timer tim; tim.start();
    LinBox::rankin(r, A, LinBox::Method::SparseElimination() );
    tim.stop();
    F2.write(std::cout << "Rank over ") << " is " << r << ' ' << tim << std::endl;
    return r;
}
    


unsigned long& LRank(unsigned long& r, char * filename, Integer p) 
{

    Integer maxmod16; LinBox::FieldTraits<LinBox::GivaroZpz<Std16> >::maxModulus(maxmod16);
    Integer maxmod32; LinBox::FieldTraits<LinBox::GivaroZpz<Std32> >::maxModulus(maxmod32);
    Integer maxmod53; LinBox::FieldTraits<LinBox::Modular<double> >::maxModulus(maxmod53);
    Integer maxmod64; LinBox::FieldTraits<LinBox::GivaroZpz<Std64> >::maxModulus(maxmod64);
    if (p == 2) {
        LinBox::GF2 F2;
        return TempLRank(r, filename, F2);
    } else if (p <= maxmod16) {
        typedef LinBox::GivaroZpz<Std16> Field;
        Field F(p);
        return TempLRank(r, filename, F);
    } else if (p <= maxmod32) {
        typedef LinBox::GivaroZpz<Std32> Field;
        Field F(p);
        return TempLRank(r, filename, F);
    } else if (p <= maxmod53) {
        typedef LinBox::Modular<double> Field;
        Field F(p);
        return TempLRank(r, filename, F);
    } else if (p <= maxmod64) {
        typedef LinBox::GivaroZpz<Std64> Field;
        Field F(p);
        return TempLRank(r, filename, F);
    } else {
        typedef LinBox::GivaroZpz<Integer> Field;
        Field F(p);
        return TempLRank(r, filename, F);
    }
    return r;
}

std::vector<size_t>& PRank(std::vector<size_t>& ranks, char * filename, Integer p, size_t e, size_t intr) 
{
    Integer maxmod;
    LinBox::FieldTraits<LinBox::GivaroZpz<Std64> >::maxModulus(maxmod);
    if (p <= maxmod) {
        typedef LinBox::GivaroZpz<Std64> Ring;
        int64 lp(p);
        Integer q = pow(p,e); int64 lq(q);
        if (q > Integer(lq)) {
            std::cerr << "Sorry power rank mod large composite not yet implemented" << std::endl;
            q = p;
            do {
                q *= p; lq = (int64)q;
            } while (q == Integer(lq));
            q/=p; lq = (int64)q;
            std::cerr << "Trying: " << lq << std::endl;
        }
        Ring F(lq);
        std::ifstream input(filename);
        LinBox::MatrixStream<Ring> ms( F, input );
        LinBox::SparseMatrix<Ring, LinBox::Vector<Ring>::SparseSeq > A (ms);
        input.close();
        LinBox::PowerGaussDomain< Ring > PGD( F );
        
        PGD.prime_power_rankin( lq, lp, ranks, A, A.rowdim(), A.coldim(), std::vector<size_t>());
        F.write(std::cout << "Ranks over ") << " are " ;
        for(std::vector<size_t>::const_iterator rit=ranks.begin(); rit != ranks.end(); ++rit)
            std::cout << *rit << ' ';
        std::cout << std::endl;
    } else {
        std::cerr << "Sorry power rank mod large composite not yet implemented" << std::endl;
        std::cerr << "Assuming integer rank" << std::endl;
        ranks.resize(0); ranks.push_back(intr);
    }               
    return ranks;
}




using namespace LinBox;

int main (int argc, char **argv)
{
//     commentator.setMaxDetailLevel (-1);
//     commentator.setMaxDepth (-1);
//     commentator.setReportStream (std::cerr);


	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: valence <matrix-file-in-supported-format> [-ata|-aat|valence] [coprime]" << std::endl;
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

        Timer chrono; chrono.start();
        if (argc >= 3) {
            Transpose<Blackbox> T(&A);
            if (strcmp(argv[2],"-ata") == 0) {
                Compose< Transpose<Blackbox>, Blackbox > C (&T, &A);
                std::cout << "A^T A is " << C.rowdim() << " by " << C.coldim() << std::endl;
                valence(val_A, C);
            } else if (strcmp(argv[2],"-aat") == 0) {
                Compose< Blackbox, Transpose<Blackbox> > C (&A, &T);
                std::cout << "A A^T is " << C.rowdim() << " by " << C.coldim() << std::endl;
                valence(val_A, C);                
            } else {
                std::cout << "Suppose primes are contained in " << argv[2] << std::endl;
                val_A = Integer(argv[2]);
            }
        } else {
            if (A.rowdim() != A.coldim()) {
                std::cerr << "Valence works only on square matrices, try either to change the dimension in the matrix file, or to compute the valence of A A^T or A^T A, via the -aat or -ata options."  << std::endl;
                exit(0);
            } else
                valence (val_A, A);
        }
        
        std::cout << "Valence is " << val_A << std::endl;

	std::vector<Integer> Moduli;
	std::vector<size_t> exponents;
	IntFactorDom<> FTD;

        typedef std::pair<Integer,unsigned long> PairIntRk;
        std::vector< PairIntRk > smith;
        

	Integer coprimeV=2;
        if (argc >= 4) {
            coprimeV = Integer(argv[3]);
        }
	while ( gcd(val_A,coprimeV) > 1 ) {
		FTD.nextprimein(coprimeV);
	}

        if (argc >= 4) {
            std::cout << "Suppose " << argv[3] << " is coprime with Smith form" << std::endl;
        }
        
	std::cout << "Integer rank: " << std::endl;
	
	unsigned long coprimeR; LRank(coprimeR, argv[1], coprimeV);
        smith.push_back(PairIntRk(coprimeV, coprimeR));
//         std::cerr << "Rank mod " << coprimeV << " is " << coprimeR << std::endl;

	std::cout << "Some factors (5000 factoring loop bound): ";
	FTD.set(Moduli, exponents, val_A, 5000);
        std::vector<size_t>::const_iterator eit=exponents.begin();
	for(std::vector<Integer>::const_iterator mit=Moduli.begin();
	    mit != Moduli.end(); ++mit,++eit)
            std::cout << *mit << '^' << *eit << ' ';
	std::cout << std::endl;
	
        std::vector<Integer> SmithDiagonal(coprimeR,Integer(1));


	for(std::vector<Integer>::const_iterator mit=Moduli.begin();
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
                if (*eit > 1) {
                    PRank(ranks, argv[1], sit->first, *eit, coprimeR);
                } else {
                    PRank(ranks, argv[1], sit->first, 2, coprimeR);
                }
                if (ranks.size() == 1) ranks.push_back(coprimeR);
                for(size_t expo = (*eit)<<1; ranks.back() < coprimeR; expo<<=1) {
                    PRank(ranks, argv[1], sit->first, expo, coprimeR);
                    if (ranks.size() < expo) {
                        std::cerr << "Larger prime power not yet implemented" << std::endl;
                        break;
                    }
                }
                std::vector<size_t>::const_iterator rit=ranks.begin();
                unsigned long modrank = *rit;
                for(++rit; rit!= ranks.end(); ++rit) {
                    if ((*rit)>= coprimeR) break;
                    for(size_t i=(*rit); i < coprimeR; ++i)
                        SmithDiagonal[i] *= sit->first;
                    modrank = *rit;
                }
            }
	}
	
        Integer si=1;
        size_t num=0;
        for( std::vector<Integer>::const_iterator dit=SmithDiagonal.begin();
             dit != SmithDiagonal.end(); ++dit) {
            if (*dit == si) ++num;
            else {
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
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
