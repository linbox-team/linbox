/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2010 LinBox
 * Written by <brice.boyer@imag.fr>
 * Adapted from nullspaceMP in Storjohann and Chen's IML.
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*! @file tests/test-nullspaceZ.C
 * @ingroup tests
 * @brief We test various implementation of dense integer nullsapce.
 * In details :
 *  - ?
 *  - ?
 *  - ?
 */


#include "linbox/linbox-config.h"
#include "linbox/field/modular.h"
#include "linbox/integer.h"
#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/randiter/generic.h"
// #include "linbox/algorithms/NullspaceZ.h"
#include "linbox/util/timer.h"

#include "test-matrix-utils.h"

#ifdef __LINBOX_HAVE_IML
#include "linbox/util/iml_wrapper.h"
#endif

#include "linbox/algorithms/iml.h"

using namespace LinBox;

template <class Field>
bool testIMLrank(const Field &F, size_t m, size_t n, size_t rank, int iterations)
{
	n += (size_t)Integer::random_lessthan(Integer(10)) ;
	m += (size_t)Integer::random_lessthan(Integer(10)) ;
	rank += (size_t)Integer::random_lessthan(Integer(7)) ;

	typedef typename Field::Element			Element;
	commentator().start ("Testing Rank from IML","testIMLrank",(unsigned int)iterations);
	bool ret = true;
	{
		size_t min = std::min(m,n);
		if (rank > min ) {
			rank = min; // rank <= min(m,n)...
		}
		if (min > 0)
			--min;
	}

	for (int k=0; k<iterations; ++k) {

		commentator().progress(k);
		BlasMatrix<Field> A(F,m,n);
		RandomMatrixWithRank(F,A.getWritePointer(),m,n,rank);
		assert(CheckRank(F,A.getWritePointer(),m,n,n,rank));

#ifdef __LINBOX_HAVE_IML
		BlasMatrix<Field> B(F,m,n);
		B.copy(A);
		long jrank = IML::mRank(F.characteristic(),B.getWritePointer(),A.rowdim(),A.coldim());
		if (jrank != (long)rank){
			// std::cout << jrank <<"!=" <<rank << std::endl;
			ret = false;
			break;
		}
#endif
		size_t irank = iml::mRank(A);
		if (irank != rank) {
			// std::cout << irank <<"!=" <<rank << std::endl;
			ret = false;
			break;
		}
		}

	commentator().stop (MSG_STATUS(ret),"testIMLrank");

	return ret;
}

#if 0
template <class Field>
bool testIMLnullspace(const Field &F, size_t m, size_t n, size_t rank, int iterations)
{
	n += (size_t)Integer::random_lessthan(Integer(10)) ;
	m += (size_t)Integer::random_lessthan(Integer(10)) ;
	rank += (size_t)Integer::random_lessthan(Integer(7)) ;

	typedef typename Field::Element			Element;
	commentator().start ("Testing Rank from IML","testNullSpace",(unsigned int)iterations);
	bool ret = true;
	{
		size_t min = std::min(m,n);
		if (rank > min ) {
			rank = min; // rank <= min(m,n)...
		}
		if (min > 0)
			--min;
	}

	for (int k=0; k<iterations; ++k) {

		commentator().progress(k);
		BlasMatrix<Field> A(F,m,n);
		RandomMatrixWithRank(F,A.getWritePointer(),m,n,rank);
		assert(CheckRank(F,A.getWritePointer(),m,n,n,rank));

#ifdef __LINBOX_HAVE_IML
		BlasMatrix<Field> B(F,m,n);
		B.copy(A);
		long jrank = IML::mRank(F.characteristic(),B.getWritePointer(),A.rowdim(),A.coldim());
		if (jrank != (long)rank){
			// std::cout << jrank <<"!=" <<rank << std::endl;
			ret = false;
			break;
		}
#endif
		size_t irank = iml::mRank(A);
		if (irank != rank) {
			// std::cout << irank <<"!=" <<rank << std::endl;
			ret = false;
			break;
		}
		}


	return ret;
}
#endif

template <class Field>
bool testIMLinverse(const Field &F, size_t m, int iterations)
{
	m += (size_t)Integer::random_lessthan(Integer(10)) ;
	typedef typename Field::Element			Element;
	commentator().start ("Testing det/inv from IML","testIMLinverse",(unsigned int)iterations);
	bool ret = true;

	typename Field::RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);

	BlasMatrixDomain<Field> BMD(F);

	for (int k=0; k<iterations; ++k) {

		commentator().progress(k);
		BlasMatrix<Field> A(F,m,m);
		Element det ;
	       	Gn.random(det);
		RandomMatrixWithDet(F,A.getWritePointer(),m,m,det);
		assert(CheckDet(F,A.getWritePointer(),m,m,det));

		BlasMatrix<Field> Ac(F,m,m);
		BlasMatrix<Field> Ad(F,m,m);
		Ac.copy(A);
		Ad.copy(A);

		Element idet = iml::mDeterminant(Ac);
		if (idet != det) {
			std::cout << idet <<"!=" <<det << std::endl;
			ret = false;
			break;
		}

		iml::mInverseIn(Ad);
		BlasMatrix<Field> Id(F,m,m);
		BMD.mul(Id,A,Ad);
		if (!BMD.isIdentity(Id)) {
			// A.write(std::cout<<"A:=",true)<<";" <<std::endl;
			// Ad.write(std::cout<<"B:=",true)<<";" <<std::endl;
			// Id.write(std::cout<<"Id:=",true)<<";" <<std::endl;
			ret = false ;
			break;
		}

#ifdef __LINBOX_HAVE_IML
		BlasMatrix<Field> Bc(F,m,m);
		Bc.copy(A);
		long jdet = IML::mDeterminant(F.characteristic(),Bc.getWritePointer(),A.rowdim());
		if (jdet != (long)det){
			// std::cout << jdet <<"!=" <<det << std::endl;
			ret = false;
			break;
		}
		BlasMatrix<Field> Bd(F,m,m);
		Bd.copy(A);
		IML::mInverse(F.characteristic(),Bd.getWritePointer(),Bd.rowdim());
		if (! BMD.areEqual(Bd,Ad)) {
			ret = false;
		}
#endif

	}

	commentator().stop (MSG_STATUS(ret),"testIMLinverse");
	return ret;
}

int main(int argc, char ** argv)
{
	typedef Modular<double> Field;
	//-----------------------------------------------------------------------
	// Choice of the finite field representation
	//typedef GivaroZpz<Std32> Field;
	typedef Modular<double> Field;
	//typedef Modular<float> Field;
	//typedef Modular<uint32_t> Field;
	//------------------------------------------------------------------------

	bool pass = true;

	static size_t n = 10;
	static size_t m = 10;
	static size_t r = 3;
	static integer q = 101;
	static int iterations =4;

	static Argument args[] = {
		{ 'n', "-n N", "Set width of test matrices.",			TYPE_INT,     &n },
		{ 'm', "-m M", "Set hight of test matrices.",			TYPE_INT,     &m },
		{ 'r', "-r R", "Set rank of test matrices.",			TYPE_INT,     &r },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].",		TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",           TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	commentator().start("IML/iml test suite", "iml");
	std::ostream& report = commentator().report();

	Field F (q);


	TESTE("IML rank");
	if (!testIMLrank (F, m,n,r, iterations))
		pass=false;
	RAPPORT("IML rank");

	TESTE("IML det");
	if (!testIMLinverse (F, m, iterations))
		pass=false;
	RAPPORT("IML det");


	commentator().stop(MSG_STATUS (pass),"IML/iml test suite");
	return (pass ? 0 : -1);

}


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
