/*
 * algorithms/smith-form-valence.h
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

/**\file algorithms/smith-form-valence.h
 * @example  algorithms/smith-form-valence.h
 \brief Valence of sparse matrix over Z or Zp.
 \ingroup examples
*/

#include <givaro/modular.h>
#include <givaro/givintnumtheo.h>
#include <givaro/gf2.h>
#include <linbox/field/field-traits.h>
#include <linbox/blackbox/transpose.h>
#include <linbox/blackbox/compose.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/rank.h>
#include <linbox/solutions/valence.h>
#include <linbox/solutions/smith-form.h>
#include <linbox/solutions/valence.h>
#include <linbox/algorithms/smith-form-sparseelim-local.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/util/timer.h>
#include <linbox/util/error.h>

#include <string>

#ifndef __VALENCE_FACTOR_LOOPS__
#define __VALENCE_FACTOR_LOOPS__ 50000
#endif

#ifndef __VALENCE_REPORTING__
# ifdef _LB_DEBUG
#  define __VALENCE_REPORTING__ 1
# else
#  define __VALENCE_REPORTING__ 0
# endif
#endif

#ifdef __LINBOX_USE_OPENMP
#include <omp.h>
#define THREAD_NUM omp_get_thread_num()
#else
#define THREAD_NUM 1
#endif

#if __VALENCE_REPORTING__
#include <sstream>
#include <iostream>
#endif

namespace LinBox {

template<class Field>
size_t& TempLRank(size_t& r, const char * filename, const Field& F)
{
	std::ifstream input(filename);
	MatrixStream< Field > msf( F, input );
	SparseMatrix<Field,SparseMatrixFormat::SparseSeq> FA(msf);
	input.close();
	Timer tim; tim.start();
	rankInPlace(r, FA);
	tim.stop();
	if (__VALENCE_REPORTING__) {
        std::ostringstream report;
        F.write(report << "Rank over ") << " is " << r << ' ' << tim <<  " on T" << THREAD_NUM << std::endl;
        std::clog << report.str();
    }
	return r;
}

size_t& TempLRank(size_t& r, const char * filename, const GF2& F2)
{
	std::ifstream input(filename);
	ZeroOne<GF2> A;
	A.read(input);
	input.close();

	Timer tim; tim.start();
	rankInPlace(r, A, Method::SparseElimination() );
	tim.stop();
	if (__VALENCE_REPORTING__) {
        std::ostringstream report;
        F2.write(report << "Rank over ") << " is " << r << ' ' << tim <<  " on T" << THREAD_NUM << std::endl;
        std::clog << report.str();
    }
	return r;
}

size_t& LRank(size_t& r, const char * filename,Givaro::Integer p)
{

	Givaro::Integer maxmod16; FieldTraits<Givaro::Modular<int16_t> >::maxModulus(maxmod16);
	Givaro::Integer maxmod32; FieldTraits<Givaro::Modular<int32_t> >::maxModulus(maxmod32);
	Givaro::Integer maxmod53; FieldTraits<Givaro::Modular<double> >::maxModulus(maxmod53);
	Givaro::Integer maxmod64; FieldTraits<Givaro::Modular<int64_t> >::maxModulus(maxmod64);
	if (p == 2) {
		GF2 F2;
		return TempLRank(r, filename, F2);
	}
	else if (p <= maxmod16) {
		typedef Givaro::Modular<int16_t> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else if (p <= maxmod32) {
		typedef Givaro::Modular<int32_t> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else if (p <= maxmod53) {
		typedef Givaro::Modular<double> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else if (p <= maxmod64) {
		typedef Givaro::Modular<int64_t> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else {
		typedef Givaro::Modular<Givaro::Integer> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	return r;
}

std::vector<size_t>& PRank(std::vector<size_t>& ranks, size_t& effective_exponent, const char * filename,Givaro::Integer p, size_t e, size_t intr)
{
#if __VALENCE_REPORTING__
    std::ostringstream logreport;
#endif
	effective_exponent = e;
	Givaro::Integer maxmod;
	FieldTraits<Givaro::Modular<int64_t> >::maxModulus(maxmod);
	if (p <= maxmod) {
		typedef Givaro::Modular<int64_t> Ring;
		int64_t lp(p);
		Givaro::Integer q = pow(p,uint64_t(e)); int64_t lq(q);
		if (q > Ring::maxCardinality()) {
#if __VALENCE_REPORTING__
                logreport << "Power rank might need extra large composite (" << p << '^' << e << ")." << std::endl;
#endif
			q = p;
			for(effective_exponent=1; q <= Ring::maxCardinality(); ++effective_exponent) {
				q *= p;
			}
			q/=p; --effective_exponent;
			lq = (int64_t)q;
            if (effective_exponent <= 1) {
                    // Not able to use Modular<int64_t>,
                    // modulus is ok, but is already too large when squared
                    // Return that no prime power was performed
#if __VALENCE_REPORTING__
                {
                    logreport << "Exceeding int64_t ... power rank useless, nothing done." << std::endl;
                    std::clog << logreport.str();
                }
#endif
                effective_exponent=1;
                return ranks;
            }
#if __VALENCE_REPORTING__
                logreport << "First trying: " << lq << " (=" << p << '^' << effective_exponent << ", without further warning this will be sufficient)." << std::endl;
#endif
		}
		Ring F(lq);
		std::ifstream input(filename);
		MatrixStream<Ring> ms( F, input );
		SparseMatrix<Ring,SparseMatrixFormat::SparseSeq > A (ms);
		input.close();

		PowerGaussDomain< Ring > PGD( F );
        Permutation<Ring> Q(F,A.coldim());

		Timer tim; tim.clear(); tim.start();
		PGD.prime_power_rankin( lq, lp, ranks, A, Q, A.rowdim(), A.coldim(), std::vector<size_t>());
		tim.stop();
#if __VALENCE_REPORTING__
		{
			F.write(logreport << "Ranks over ") << " are " ;
			for(auto const & rit:ranks) logreport << rit << ' ';
			logreport << ' ' << tim <<  " on T" << THREAD_NUM << std::endl;
		}
#endif
	} else {
            // Not able to use Modular<int64_t>, modulus too large
            // Return that no prime power was performed
#if __VALENCE_REPORTING__
            logreport << "Exceeding int64_t ... even for the prime, nothing done." << std::endl;
#endif
        effective_exponent=1;
	}

#if __VALENCE_REPORTING__
        std::clog << logreport.str();
#endif
	return ranks;
}
}

#include <linbox/field/gf2.h>
#include <linbox/algorithms/smith-form-sparseelim-poweroftwo.h>

namespace LinBox {

std::vector<size_t>& PRankPowerOfTwo(std::vector<size_t>& ranks, size_t& effective_exponent, const char * filename, size_t e, size_t intr)
{
#if __VALENCE_REPORTING__
    std::ostringstream logreport;
#endif
	effective_exponent = e;
	if (e > 63) {
#if __VALENCE_REPORTING__
		{
            logreport << "Power rank power of two might need extra large composite (2^" << e << ")." << std::endl;
            logreport << "First trying: 63, without further warning this will be sufficient)." << std::endl;
        }
#endif
		effective_exponent = 63;
	}


	typedef Givaro::ZRing<int64_t> Ring;
	Ring F;
	std::ifstream input(filename);
	MatrixStream<Ring> ms( F, input );
	SparseMatrix<Ring,SparseMatrixFormat::SparseSeq > A (ms);
	input.close();
	PowerGaussDomainPowerOfTwo< uint64_t > PGD;
    GF2 F2;
    Permutation<GF2> Q(F2,A.coldim());

	Timer tim; tim.clear(); tim.start();
	PGD.prime_power_rankin( effective_exponent, ranks, A, Q, A.rowdim(), A.coldim(), std::vector<size_t>());
	tim.stop();
#if __VALENCE_REPORTING__
	{
		F.write(logreport << "Ranks over ") << " modulo 2^" << effective_exponent << " are " ;
		for(auto const& rit: ranks) logreport << rit << ' ';
		logreport<< ' ' << tim <<  " on T" << THREAD_NUM << std::endl;
        std::clog << logreport.str();
	}
#endif
	return ranks;
}

std::vector<size_t>& PRankInteger(std::vector<size_t>& ranks, const char * filename,Givaro::Integer p, size_t e, size_t intr)
{
	typedef Givaro::Modular<Givaro::Integer> Ring;
	Givaro::Integer q = pow(p,uint64_t(e));
	Ring F(q);
	std::ifstream input(filename);
	MatrixStream<Ring> ms( F, input );
	SparseMatrix<Ring,SparseMatrixFormat::SparseSeq > A (ms);
	input.close();
	PowerGaussDomain< Ring > PGD( F );
    Permutation<Ring> Q(F,A.coldim());

	Timer tim; tim.clear(); tim.start();
	PGD.prime_power_rankin( q, p, ranks, A, Q, A.rowdim(), A.coldim(), std::vector<size_t>());
	tim.stop();
	if (__VALENCE_REPORTING__) {
        std::ostringstream logreport;
		F.write(logreport << "Ranks over ") << " are " ;
		for(auto const & rit: ranks) logreport << rit << ' ';
		logreport << ' ' << tim << std::endl;
        std::clog << logreport.str();
	}
	return ranks;
}

std::vector<size_t>& PRankIntegerPowerOfTwo(std::vector<size_t>& ranks, const char * filename, size_t e, size_t intr)
{
	typedef Givaro::ZRing<Givaro::Integer> Ring;
	Ring ZZ;
	std::ifstream input(filename);
	MatrixStream<Ring> ms( ZZ, input );
	SparseMatrix<Ring,SparseMatrixFormat::SparseSeq > A (ms);
	input.close();
	PowerGaussDomainPowerOfTwo< Givaro::Integer > PGD;
    Permutation<Ring> Q(ZZ, A.coldim());

	Timer tim; tim.clear(); tim.start();
	PGD.prime_power_rankin( e, ranks, A, Q, A.rowdim(), A.coldim(), std::vector<size_t>());
	tim.stop();
	if (__VALENCE_REPORTING__) {
        std::ostringstream logreport;
		ZZ.write(logreport << "Ranks over ") << " modulo 2^" << e << " are " ;
		for(auto const& rit: ranks) logreport << rit << ' ';
		logreport << ' ' << tim << std::endl;
	}
	return ranks;
}


typedef std::pair<Givaro::Integer,size_t> PairIntRk;


std::vector<size_t>& AllPowersRanks(
    std::vector<size_t>& ranks,
    const Givaro::Integer& squarefreePrime,// smith[j].first
    const size_t& squarefreeRank,// smith[j].second
    const size_t& exponentBound,	// exponents[j]
    const size_t& coprimeRank,		// coprimeR
    const char * filename) {		// argv[1]

    if (squarefreeRank != coprimeRank) {

        ranks.push_back(squarefreeRank);
        size_t effexp;
        if (exponentBound > 1) {
                // See if a not too small, not too large exponent would work
                // Usually, closest to word size
            if (squarefreePrime == 2)
                PRankPowerOfTwo(ranks, effexp, filename, exponentBound, coprimeRank);
            else
                PRank(ranks, effexp, filename, squarefreePrime, exponentBound, coprimeRank);
        } else {
                // Square does not divide valence
                // Try first with the smallest possible exponent: 2
            if (squarefreePrime == 2)
                PRankPowerOfTwo(ranks, effexp, filename, 2, coprimeRank);
            else
                PRank(ranks, effexp, filename, squarefreePrime, 2, coprimeRank);
        }

        if (effexp < exponentBound) {
                // Above report shows that more powers are needed,
                // try successive doublings Over abitrary precision
            for(size_t expo = effexp<<1; ranks.back() < coprimeRank; expo<<=1) {
                if (squarefreePrime == 2)
                    PRankIntegerPowerOfTwo(ranks, filename, expo, coprimeRank);
                else
                    PRankInteger(ranks, filename, squarefreePrime, expo, coprimeRank);
            }
        } else {
                // Larger exponents are needed
                // Try first small precision, then arbitrary
            for(size_t expo = (exponentBound)<<1; ranks.back() < coprimeRank; expo<<=1) {
                if (squarefreePrime == 2)
                    PRankPowerOfTwo(ranks, effexp, filename, expo, coprimeRank);
                else
                    PRank(ranks, effexp, filename, squarefreePrime, expo, coprimeRank);
                if (ranks.size() < expo) {
                    if (__VALENCE_REPORTING__)
                        std::clog << "It seems we need a larger prime power, it will take longer ...\n" << std::flush;
                        // break;
                    if (squarefreePrime == 2)
                        PRankIntegerPowerOfTwo(ranks, filename, expo, coprimeRank);
                    else
                        PRankInteger(ranks, filename, squarefreePrime, expo, coprimeRank);
                }
            }
        }
    }

    return ranks;
}

std::vector<Givaro::Integer>& populateSmithForm(
    std::vector<Givaro::Integer>& SmithDiagonal,
    const std::vector<size_t>& ranks,
    const Givaro::Integer& squarefreePrime,// smith[j].first
    const size_t& squarefreeRank,// smith[j].second
    const size_t& coprimeRank) {	// coprimeR


    SmithDiagonal.resize(coprimeRank, Givaro::Integer(1) );

    for(size_t i=squarefreeRank; i < coprimeRank; ++i) {
        SmithDiagonal[i] *= squarefreePrime;
    }
    auto rit=ranks.begin(); for(++rit; rit!= ranks.end(); ++rit) {
        if ((*rit)>= coprimeRank) break;
        for(size_t i=(*rit); i < coprimeRank; ++i)
            SmithDiagonal[i] *= squarefreePrime;
    }

    return SmithDiagonal;
}

template<class Blackbox>
std::vector<Givaro::Integer>& smithValence(std::vector<Givaro::Integer>& SmithDiagonal,
                                           Givaro::Integer& valence,
                                           const Blackbox& A,
                                           const std::string& filename,
                                           Givaro::Integer& coprimeV,
                                           size_t method=0) {
        // method for valence squarization:
		//	0 for automatic, 1 for aat, 2 for ata
        // Blackbox provides the Integer matrix rereadable from filename
        // if valence != 0:
		//	then the valence is not computed and the parameter is used
        // if coprimeV != 1:
		//  then this value is supposed to be coprime with the valence

    if (__VALENCE_REPORTING__)
        std::clog << "sV threads: " << NUM_THREADS << std::endl;

    if (valence == 0) {
        squarizeValence(valence, A, method);
    }

    if (__VALENCE_REPORTING__)
        std::clog << "Valence is " << valence << std::endl;

    std::vector<Givaro::Integer> Moduli;
	std::vector<size_t> exponents;

    if (__VALENCE_REPORTING__)
        std::clog << "Some factors (" << __VALENCE_FACTOR_LOOPS__ << " factoring loop bound): ";

    Givaro::IntFactorDom<> FTD;
    FTD.set(Moduli, exponents, valence, __VALENCE_FACTOR_LOOPS__);

    if (__VALENCE_REPORTING__) {
        auto eit=exponents.begin();
        for(auto const &mit: Moduli) std::clog << mit << '^' << *eit++ << ' ';
        std::clog << std::endl;
    }

	std::vector< size_t > smith(Moduli.size());

    if (coprimeV == 1) {
        coprimeV=2;
        while ( gcd(valence,coprimeV) > 1 ) {
            FTD.nextprimein(coprimeV);
        }
    }

    size_t coprimeR;
    std::vector<std::vector<size_t> > AllRanks(Moduli.size());

    for(size_t j=0; j<Moduli.size(); ++j) {
        { TASK(MODE(CONSTREFERENCE(Moduli,smith,filename) WRITE(smith[j]) ),
        {
            LRank(smith[j], filename.c_str(), Moduli[j]);
        })}
    }

//     { TASK(MODE(CONSTREFERENCE(coprimeV,filename) WRITE(coprimeR) ),
//     {
        LRank(coprimeR, filename.c_str(), coprimeV);
//     })}

    WAIT;

    SYNCH_GROUP(
        for(size_t j=0; j<Moduli.size(); ++j) {
            { TASK(MODE(CONSTREFERENCE(smith,Moduli,AllRanks,filename,coprimeR,exponents)
                        WRITE(AllRanks[j])),
            {
                AllPowersRanks(AllRanks[j], Moduli[j], smith[j], exponents[j],
                               coprimeR, filename.c_str());
            })}
        }
    )

    for(size_t j=0; j<Moduli.size(); ++j) {
        if (smith[j] != coprimeR) {
            populateSmithForm(SmithDiagonal, AllRanks[j], Moduli[j], smith[j], coprimeR);
        }
    }

    return SmithDiagonal;
}

template<class Blackbox>
std::vector<Givaro::Integer>& smithValence(
    std::vector<Givaro::Integer>& SmithDiagonal,
    const Blackbox& A,
    const std::string& filename,
    size_t method=0)
{
    Givaro::Integer valence(0);
    Givaro::Integer coprimeV(1);
    return smithValence(SmithDiagonal, valence, A, filename, coprimeV, method);
}


template<class PIR>
std::ostream& writeCompressedSmith(
    std::ostream& out,
    const std::vector<typename PIR::Element>& SmithDiagonal,
    const PIR& ZZ,
    const size_t m, const size_t n) {
        // Output is a list of pairs (integral value, number of repetitions)

    SmithList<PIR> tempSL;
    compressedSmith(tempSL, SmithDiagonal, ZZ, m, n);
    out << '(';
	for( auto sit : tempSL )
        out << '[' << sit.first << ',' << sit.second << "] ";
	return out << ')';
}

}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
