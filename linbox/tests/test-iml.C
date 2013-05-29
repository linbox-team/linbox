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
	n += (size_t) Integer::random_lessthan(Integer(10)) ;
	m += (size_t)Integer::random_lessthan(Integer(10)) ;
	rank += (size_t)Integer::random_lessthan(Integer(7)) ;

	typedef typename Field::Element			Element;
	commentator().start ("Testing Rank from IML","testIMLrank",(unsigned int)iterations);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);

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

		commentator().startIteration((unsigned int)k);
		BlasMatrix<Field> A(F,m,n);
		RandomMatrixWithRank(F,A.getWritePointer(),m,n,rank);
		assert(CheckRank(F,A.getWritePointer(),m,n,n,rank));

#ifdef __LINBOX_HAVE_IML
		BlasMatrix<Field> B(F,m,n);
		B.copy(A);
		size_t jrank = (size_t)IML::mRank((IML::FiniteField)F.characteristic(),B.getWritePointer(),
						(long)A.rowdim(),(long)A.coldim());
		if (jrank != rank){
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
		commentator().stop("done");
		commentator().progress();
		}

	commentator().stop (MSG_STATUS(ret),"testIMLrank");

	return ret;
}


#ifdef __LINBOX_HAVE_IML
// WARNING this test checks equality with IML, not that IML results are correct
template <class Field>
bool testIMLstuff(const Field &F, size_t m, size_t n, size_t rank, int iterations)
{
	n += (size_t)Integer::random_lessthan(Integer(10)) ;
	m += (size_t)Integer::random_lessthan(Integer(10)) ;
	rank += (size_t)Integer::random_lessthan(Integer(7)) ;

	typedef typename Field::Element			Element;
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);

	commentator().start ("Testing RET from IML","testIMLrowechl",(unsigned int)iterations);
	bool ret = true;
	{
		size_t min = std::min(m,n);
		if (rank > min ) {
			rank = min; // rank <= min(m,n)...
		}
		if (min > 0)
			--min;
	}

	BlasMatrixDomain<Field> BMD(F);

	for (int k=0; k<iterations; ++k) {

		commentator().startIteration((unsigned int)k);
		{ // mAdjoint
			commentator().report()<< "mAdjoint" << std::endl;
			BlasMatrix<Field> A(F,m,m);
			RandomMatrixWithRank(F,A.getWritePointer(),m,m,rank);
			assert(CheckRank(F,A.getWritePointer(),m,m,m,rank));

			BlasMatrix<Field> B(F,m,m);
			BlasMatrix<Field> C(F,m,m);
			B.copy(A);
			C.copy(A);
			BlasMatrix<Field> C1(F,m,m);
			double * B1 = IML::mAdjoint((IML::FiniteField)F.characteristic(),B.getWritePointer(),
						    (long)A.rowdim());
			iml::mAdjoint(C1,C);
			BlasMatrix<Field> B2(F,B1,m,m);
			if (! BMD.areEqual(B2,C1)){
				ret = false;
				IML_XFREE(B1);
				break;
			}
			IML_XFREE(B1);
		}

		{ // mBasis 1,1
			commentator().report()<< "mBasis 11" << std::endl;
			BlasMatrix<Field> A(F,m,n);
			RandomMatrixWithRank(F,A.getWritePointer(),m,n,rank);
			assert(CheckRank(F,A.getWritePointer(),m,n,n,rank));

			BlasMatrix<Field> B(F,m,n);
			B.copy(A);
			double ** B1 = IML_XMALLOC(IML::Double*,1), ** B2=IML_XMALLOC(IML::Double*,1) ;
			size_t r = (size_t)IML::mBasis((IML::FiniteField)F.characteristic(),B.getWritePointer(),
					       (long)B.rowdim(),(long)B.coldim(),
					     1,1,B1,B2);
			BlasMatrix<Field> C1(F,B1[0],r,B.coldim());
			BlasMatrix<Field> C2(F,B2[0],B.rowdim()-r,B.rowdim());
			BlasMatrix<Field> A1(F),A2(F);
			BlasMatrix<Field> C(F,m,n);
			C.copy(A);
			iml::mBasis(C,1,1,A1,A2);
			if (!BMD.areEqual(A1,C1)||!BMD.areEqual(A2,C2)){
				ret = false;
				IML_XFREE(*B1);
			       	IML_XFREE(*B2);
				IML_XFREE(B1);
				IML_XFREE(B2);
				commentator().report()<< ((ret==true)?"ok":"ko") << std::endl;
				break;
			}
			else {
			IML_XFREE(*B1);
			IML_XFREE(*B2);
			IML_XFREE(B1);
			IML_XFREE(B2);
			}
		}

		{ // mBasis 1,0
			commentator().report()<< "mBasis 10" << std::endl;
			BlasMatrix<Field> A(F,m,n);
			RandomMatrixWithRank(F,A.getWritePointer(),m,n,rank);
			assert(CheckRank(F,A.getWritePointer(),m,n,n,rank));

			BlasMatrix<Field> B(F,m,n);
			B.copy(A);
			double ** B1 = IML_XMALLOC(IML::Double*,1), ** B2=IML_XMALLOC(IML::Double*,1) ;
			size_t r = (size_t)IML::mBasis((IML::FiniteField)F.characteristic(),B.getWritePointer(),
					       (long)B.rowdim(),(long)B.coldim(),
					     1,0,B1,B2);
			BlasMatrix<Field> C1(F,B1[0],r,B.coldim());
			BlasMatrix<Field> A1(F),A2(F);
			BlasMatrix<Field> C(F,m,n);
			C.copy(A);
			iml::mBasis(C,1,0,A1,A2);
			if (!BMD.areEqual(A1,C1)){
				ret = false;
				IML_XFREE(*B1);
				IML_XFREE(B1);
				IML_XFREE(B2);
				commentator().report()<< ((ret==true)?("ok"):("ko")) << std::endl;
				break;
			}
			else {
			IML_XFREE(*B1);
			IML_XFREE(B1);
			IML_XFREE(B2);
			}
		}

		{ // mBasis 0,1
			commentator().report()<< "mBasis 01" << std::endl;
			BlasMatrix<Field> A(F,m,n);
			RandomMatrixWithRank(F,A.getWritePointer(),m,n,rank);
			assert(CheckRank(F,A.getWritePointer(),m,n,n,rank));

			BlasMatrix<Field> B(F,m,n);
			B.copy(A);
			double ** B1 = IML_XMALLOC(IML::Double*,1), ** B2=IML_XMALLOC(IML::Double*,1) ;
			size_t r = (size_t)IML::mBasis((IML::FiniteField)F.characteristic(),B.getWritePointer(),
					      (long) B.rowdim(),(long)B.coldim(),
					     0,1,B1,B2);
			BlasMatrix<Field> C2(F,B2[0],B.rowdim()-r,B.rowdim());
			BlasMatrix<Field> A1(F),A2(F);
			BlasMatrix<Field> C(F,m,n);
			C.copy(A);
			iml::mBasis(C,1,1,A1,A2);
			if (!BMD.areEqual(A2,C2)){
				ret = false;
				IML_XFREE(*B2);
				IML_XFREE(B1);
				IML_XFREE(B2);
				commentator().report()<< ((ret==true)?("ok"):("ko")) << std::endl;
				break;
			}
			else {
				IML_XFREE(*B2);
				IML_XFREE(B1);
			       	IML_XFREE(B2);
			}
		}

		{ // mRankProfile
			commentator().report()<< "mRankProfile" << std::endl;
			BlasMatrix<Field> A(F,m,n);
			RandomMatrixWithRank(F,A.getWritePointer(),m,n,rank);
			assert(CheckRank(F,A.getWritePointer(),m,n,n,rank));

			BlasMatrix<Field> B(F,m,n);
			B.copy(A);
			long * rp = IML::mRankProfile((IML::FiniteField)F.characteristic(),
						      B.getWritePointer(),(long)B.rowdim(),(long)B.coldim());
			BlasMatrix<Field> C(F,m,n);
			C.copy(A);
			std::vector<size_t> Rp = iml::mRankProfile(C);
			for (size_t i = 0 ; i < n+1 ; ++i){
				if (rp[i] != (long)Rp[i]) {
					ret = false;
					break;
				}
			}
			commentator().report()<< ((ret==true)?("ok"):("ko")) << std::endl;
			IML_XFREE(rp);
			if (ret == false)
				break;
		}

		commentator().stop("done");
		commentator().progress();

		}

	commentator().stop (MSG_STATUS(ret),"testIMLrank");

	return ret;
}
#endif


template <class Field>
bool testIMLinverse(const Field &F, size_t m, int iterations)
{
	m += (size_t)Integer::random_lessthan(Integer(10)) ;
	typedef typename Field::Element			Element;
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);

	commentator().start ("Testing det/inv from IML","testIMLinverse",(unsigned int)iterations);
	bool ret = true;

	typename Field::RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);

	BlasMatrixDomain<Field> BMD(F);

	for (int k=0; k<iterations; ++k) {

		commentator().startIteration((unsigned int)k);
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
			// std::cout << idet <<"!=" <<det << std::endl;
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
		size_t jdet = (size_t)IML::mDeterminant((IML::FiniteField) F.characteristic(),Bc.getWritePointer(),
						(long)A.rowdim());
		if (jdet != det){
			// std::cout << jdet <<"!=" <<det << std::endl;
			ret = false;
			break;
		}
		BlasMatrix<Field> Bd(F,m,m);
		Bd.copy(A);
		IML::mInverse((IML::FiniteField)F.characteristic(),Bd.getWritePointer(),
			      (long)Bd.rowdim());
		if (! BMD.areEqual(Bd,Ad)) {
			ret = false;
		}
#endif

	commentator().stop("done");
		commentator().progress();

	}

	commentator().stop (MSG_STATUS(ret),"testIMLinverse");
	return ret;
}

#ifdef __LINBOX_HAVE_IML
template<class Field>
bool testIMLcra(const Field &F, int iterations)
{
	bool ret = true ;
	for(int i = 1 ; i < iterations ; ++i) {
	Integer A = Integer::random_lessthan<false>(200); // <2^200
	// std::cout << A << std::endl;
	size_t p =F.characteristic();
	size_t len;

	// number to reconstuct
	mpz_t B ; mpz_init_set(B,(mpz_srcptr)(Givaro::SpyInteger::get_rep(A)));

	long basislen =1;

	//XXX compute stuff
	IML::FiniteField qh = 2<<22;
	mpz_t mp_maxInter ; mpz_init(mp_maxInter);
	mpz_ui_pow_ui(mp_maxInter,2,200);

	IML::FiniteField ** basiscmb ;
	basiscmb = IML::findRNS(qh, mp_maxInter, &basislen);
	IML::FiniteField * basis = basiscmb[0];
	IML::FiniteField * cmbasis = basiscmb[1];

	// bdcoeff
	IML::FiniteField * bdcoeff = IML::repBound(basislen, basis, cmbasis);

	// compute mp_basisprod
	mpz_t mp_basisprod; mpz_init(mp_basisprod);
	IML::basisProd(basislen, basis, mp_basisprod);

	// compute dtemp (rep of B in RNS basis)
	IML::Double *dtemp;
	dtemp = IML_XMALLOC(IML::Double, (size_t)basislen);

	mpz_t mymod;
	mpz_init(mymod);
	for (size_t l = 0 ; l< (size_t)basislen ; ++l) {
		mpz_mod_ui(mymod,B,basis[l]);
		dtemp[l] = (IML::Double)mpz_get_ui(mymod);
	}

	IML::ChineseRemainder(basislen,mp_basisprod,basis,cmbasis,
			      bdcoeff, dtemp,B);
	// std::cout << B << std::endl;
	mpz_clear(mymod);
	mpz_clear(mp_basisprod);
	mpz_clear(mp_maxInter);
	IML_XFREE(basis);
	IML_XFREE(cmbasis);
	IML_XFREE(basiscmb);
	IML_XFREE(bdcoeff);
	if  (mpz_cmp((mpz_srcptr)(Givaro::SpyInteger::get_rep(A)),B)!=0){
		std::cout << "false" << std::endl;
		ret =false ;
	mpz_clear(B);
		break;
	}
	mpz_clear(B);
	}

	return ret;
}
#endif


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
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start("IML/iml test suite", "iml");

	Field F (q);


	TESTE("IML rank");
	if (!testIMLrank (F, m,n,r, iterations))
		pass=false;
	RAPPORT("IML rank");

	TESTE("IML det");
	if (!testIMLinverse (F, m, iterations))
		pass=false;
	RAPPORT("IML det");

#ifdef __LINBOX_HAVE_IML
	TESTE("IML stuff");
	if (!testIMLstuff(F, m,n,r, iterations))
		pass=false;
	RAPPORT("IML stuff");

	TESTE("IML rns");
	if (!testIMLcra(F, iterations))
		pass=false;
	RAPPORT("IML rns");
#endif


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
