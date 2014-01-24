/*
 * examples/smithvalence.h
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

/**\file examples/smithvalence.h
 * @example  examples/smithvalence.h
 \brief Valence of sparse matrix over Z or Zp.
 \ingroup examples
 */

#include <givaro/givintnumtheo.h>
#include "linbox/field/gf2.h"
#include "linbox/field/modular.h"
#include "linbox/field/givaro.h"
#include "linbox/field/field-traits.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/compose.h"
#include "linbox/matrix/sparse.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/valence.h"
#include "linbox/algorithms/smith-form-sparseelim-local.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/util/timer.h"


template<class Field>
unsigned long& TempLRank(unsigned long& r, char * filename, const Field& F)
{
	std::ifstream input(filename);
	LinBox::MatrixStream< Field > msf( F, input );
	LinBox::SparseMatrix<Field,LinBox::SparseMatrixFormat::SparseSeq> FA(msf);
	input.close();
	LinBox::Timer tim; tim.start();
	LinBox::rankin(r, FA);
	tim.stop();
	F.write(std::cerr << "Rank over ") << " is " << r << ' ' << tim << std::endl;
	return r;
}

unsigned long& TempLRank(unsigned long& r, char * filename, const LinBox::GF2& F2)
{
	std::ifstream input(filename);
	LinBox::ZeroOne<LinBox::GF2> A;
	A.read(input);
	input.close();
	LinBox::Timer tim; tim.start();
	LinBox::rankin(r, A, LinBox::Method::SparseElimination() );
	tim.stop();
	F2.write(std::cerr << "Rank over ") << " is " << r << ' ' << tim << std::endl;
	return r;
}

unsigned long& LRank(unsigned long& r, char * filename,Givaro::Integer p)
{

	Givaro::Integer maxmod16; LinBox::FieldTraits<LinBox::GivaroZpz<Givaro::Std16> >::maxModulus(maxmod16);
	Givaro::Integer maxmod32; LinBox::FieldTraits<LinBox::GivaroZpz<Givaro::Std32> >::maxModulus(maxmod32);
	Givaro::Integer maxmod53; LinBox::FieldTraits<LinBox::Modular<double> >::maxModulus(maxmod53);
	Givaro::Integer maxmod64; LinBox::FieldTraits<LinBox::GivaroZpz<Givaro::Std64> >::maxModulus(maxmod64);
	if (p == 2) {
		LinBox::GF2 F2;
		return TempLRank(r, filename, F2);
	}
	else if (p <= maxmod16) {
		typedef LinBox::GivaroZpz<Givaro::Std16> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else if (p <= maxmod32) {
		typedef LinBox::GivaroZpz<Givaro::Std32> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else if (p <= maxmod53) {
		typedef LinBox::Modular<double> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else if (p <= maxmod64) {
		typedef LinBox::GivaroZpz<Givaro::Std64> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	else {
		typedef LinBox::GivaroZpz<Givaro::Integer> Field;
		Field F(p);
		return TempLRank(r, filename, F);
	}
	return r;
}

std::vector<size_t>& PRank(std::vector<size_t>& ranks, size_t& effective_exponent, char * filename,Givaro::Integer p, size_t e, size_t intr)
{
	effective_exponent = e;
	Givaro::Integer maxmod;
	LinBox::FieldTraits<LinBox::GivaroZpz<Givaro::Std64> >::maxModulus(maxmod);
	if (p <= maxmod) {
		typedef LinBox::GivaroZpz<Givaro::Std64> Ring;
		int64_t lp(p);
		Givaro::Integer q = pow(p,e); int64_t lq(q);
		if (q >Givaro::Integer(lq)) {
			std::cerr << "Power rank might need extra large composite (" << p << '^' << e << ")." << std::endl;
			q = p;
			for(effective_exponent=1; q <= Ring::getMaxModulus(); ++effective_exponent) {
				q *= p;
			}
			q/=p; --effective_exponent;
			lq = (int64_t)q;
			std::cerr << "First trying: " << lq << " (=" << p << '^' << effective_exponent << ", without further warning this will be sufficient)." << std::endl;
		}
		Ring F(lq);
		std::ifstream input(filename);
		LinBox::MatrixStream<Ring> ms( F, input );
		LinBox::SparseMatrix<Ring,LinBox::SparseMatrixFormat::SparseSeq > A (ms);
		input.close();
		LinBox::PowerGaussDomain< Ring > PGD( F );

		LinBox::Timer tim; tim.clear(); tim.start();
		PGD.prime_power_rankin( lq, lp, ranks, A, A.rowdim(), A.coldim(), std::vector<size_t>());
		tim.stop();
		F.write(std::cerr << "Ranks over ") << " are " ;
		for(std::vector<size_t>::const_iterator rit=ranks.begin(); rit != ranks.end(); ++rit)
			std::cerr << *rit << ' ';
		std::cerr << ' ' << tim << std::endl;
	}
	else {
		std::cerr << "*** WARNING *** Sorry power rank mod large composite not yet implemented" << std::endl;
		std::cerr << "*** WARNING *** Assuming integer rank, extra factors in the Smith form could be missing" << std::endl;
		ranks.resize(0); ranks.push_back(intr);
	}
	return ranks;
}

#include "linbox/algorithms/smith-form-sparseelim-poweroftwo.h"


std::vector<size_t>& PRankPowerOfTwo(std::vector<size_t>& ranks, size_t& effective_exponent, char * filename, size_t e, size_t intr)
{
	effective_exponent = e;
	if (e > 63) {
		std::cerr << "Power rank power of two might need extra large composite (2^" << e << ")." << std::endl;
		std::cerr << "First trying: 63, without further warning this will be sufficient)." << std::endl;
		effective_exponent = 63;
	}

	std::ifstream input(filename);
	typedef LinBox::UnparametricField<int64_t> Ring;
	Ring F;
	LinBox::MatrixStream<Ring> ms( F, input );
	LinBox::SparseMatrix<Ring,LinBox::SparseMatrixFormat::SparseSeq > A (ms);
	input.close();
	LinBox::PowerGaussDomainPowerOfTwo< uint64_t > PGD;

	LinBox::Timer tim; tim.clear(); tim.start();
	PGD.prime_power_rankin( effective_exponent, ranks, A, A.rowdim(), A.coldim(), std::vector<size_t>());
	tim.stop();
	std::cerr << "Ranks over 2^" << effective_exponent << " are " ;
	for(std::vector<size_t>::const_iterator rit=ranks.begin(); rit != ranks.end(); ++rit)
		std::cerr << *rit << ' ';
	std::cerr << ' ' << tim << std::endl;
	return ranks;
}

std::vector<size_t>& PRankInteger(std::vector<size_t>& ranks, char * filename,Givaro::Integer p, size_t e, size_t intr)
{
	typedef LinBox::GivaroZpz<Givaro::Integer> Ring;
	Givaro::Integer q = pow(p,e);
	Ring F(q);
	std::ifstream input(filename);
	LinBox::MatrixStream<Ring> ms( F, input );
	LinBox::SparseMatrix<Ring,LinBox::SparseMatrixFormat::SparseSeq > A (ms);
	input.close();
	LinBox::PowerGaussDomain< Ring > PGD( F );

	LinBox::Timer tim; tim.clear(); tim.start();
	PGD.prime_power_rankin( q, p, ranks, A, A.rowdim(), A.coldim(), std::vector<size_t>());
	tim.stop();
	F.write(std::cerr << "Ranks over ") << " are " ;
	for(std::vector<size_t>::const_iterator rit=ranks.begin(); rit != ranks.end(); ++rit)
		std::cerr << *rit << ' ';
	std::cerr << ' ' << tim << std::endl;
	return ranks;
}

std::vector<size_t>& PRankIntegerPowerOfTwo(std::vector<size_t>& ranks, char * filename, size_t e, size_t intr)
{
	typedef LinBox::PID_integer Ring;
	Ring ZZ;
	std::ifstream input(filename);
	LinBox::MatrixStream<Ring> ms( ZZ, input );
	LinBox::SparseMatrix<Ring,LinBox::SparseMatrixFormat::SparseSeq > A (ms);
	input.close();
	LinBox::PowerGaussDomainPowerOfTwo< Givaro::Integer > PGD;

	LinBox::Timer tim; tim.clear(); tim.start();
	PGD.prime_power_rankin( e, ranks, A, A.rowdim(), A.coldim(), std::vector<size_t>());
	tim.stop();
	std::cerr << "Ranks over 2^" << e << " are " ;
	for(std::vector<size_t>::const_iterator rit=ranks.begin(); rit != ranks.end(); ++rit)
		std::cerr << *rit << ' ';
	std::cerr << ' ' << tim << std::endl;
	return ranks;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
