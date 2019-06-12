/* linbox/algorithms/hybrid-det.h
 * Copyright (C) 2005 Anna Urbanska
 *
 * Written by Anna Urbanska <aniau@astronet.pl>
 * Modified by JGD <Jean-Guillaume.Dumas@imag.fr>
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_hybrid_det_H
#define __LINBOX_hybrid_det_H

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/matrix/sparse-matrix.h"
#include "givaro/zring.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/ring/modular.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/solutions/det.h"

// #define _LB_H_DET_TIMING

namespace LinBox
{

	template <class Blackbox, class MyMethod>
	struct IntegerModularDetReduced {
	private:

		const                        Blackbox &A;
		const                        MyMethod &M;
		/* contains the factor by which to divide det */
		integer                             beta;

		size_t                            factor;
		Givaro::ZRing<Integer>                           ZZ; //! @todo is it \c A.field()?
		size_t                        iter_count;
		size_t                       iter_count2;
		typedef  BlasVector<Givaro::ZRing<Integer> >  IVect ;
		IVect                            moduli ;

	public:

		IVect                            primes ;

		size_t iterations()
		{
			return iter_count;
		}

		size_t iterations2()
		{
			return iter_count2;
		}


		IntegerModularDetReduced(const Blackbox& b, const MyMethod& n, const integer& divisor, const size_t& fs) :
			A(b), M(n)
			, beta(divisor)
			, factor(fs)
			,ZZ(Givaro::ZRing<Integer>())
			,moduli(ZZ,fs),primes(ZZ,fs)
		{
			// moduli.resize(factor);
			// primes.resize(factor);
			iter_count = 0;
			iter_count2 = 0;

		}

		template<typename Field>
		IterationResult operator()(typename Field::Element& d, const Field& F)
		{

			if (beta > 1) {
				if (iter_count2 < factor) {
					Field D(primes[iter_count2]);
					typename Field::Element pbeta;
					typename Field::Element kbeta;
					typename Field::Element current_moduli;
					D.init(pbeta, beta);
					D.init(current_moduli, moduli[iter_count2]);
					D.div(kbeta, current_moduli,pbeta);
					++this->iter_count2;
					d=kbeta;
					return IterationResult::CONTINUE;
				}
			}

			typedef typename Blackbox::template rebind<Field>::other FBlackbox;
			FBlackbox Ap(A,F);
			detInPlace( d, Ap, M);

			if (beta > 1) {
				typename Field::Element y;
				F.init(y,beta);
				F.div(d,d,y);
			}

			if (iter_count < factor) {
				moduli[iter_count] = d;
			}
			++this->iter_count;

			return IterationResult::CONTINUE;
		}

		void Beta(Integer& b) { beta = b; iter_count2=0;}

	};

	/** \brief Compute the determinant of A over the integers
	 *
	 * The determinant of a linear operator A, represented as a
	 * black box, is computed over the integers.
	 *
	 * This variant is a hybrid between Last Invariant Factor and Chinese Remaindering
	 * It performs several determinants mod primes
	 * Then switches to the LIF method, producing a factor of det.
	 * It then comes back to the CRA if necessary to compute
	 * the remaining (usually small) factor of the determinant.
	 *
	 * @param d Field element into which to store the result
	 * @param A Black box of which to compute the determinant
	 * @param tag explicit over the integers
	 * @param M may be a Method::DenseElimination (default) or a Method::Wiedemann.
	 \ingroup solutions
	 */
	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element & lif_cra_det (typename Blackbox::Field::Element         &d,
							 const Blackbox                            &A,
							 const RingCategories::IntegerTag          &tag,
							 const MyMethod                            &M)
	{
		//commentator().setReportStream(std::cout);
		typedef Givaro::ModularBalanced<double> mymodular;
		typedef typename Blackbox::Field Integers;
		typedef typename Integers::Element Integer_t;

		commentator().start ("Integer Determinant - hybrid version ", "det");
		size_t myfactor=5;
		size_t early_counter=0;

		//double a = 0.0;//0.28;//0.0013//if time(lif)/time(10lu) < a * log(lif) then calculate bonus

		Integer_t lif = 1;
		Integer_t bonus = 1;
		Integer_t beta = 1;
		d=1;

		PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<mymodular>::bestBitSize(A.coldim()));
		//cout << "prime size: " << p_size << "\n";
		CRABuilderEarlySingle<mymodular> cra(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);
		IntegerModularDetReduced<Blackbox,MyMethod> iteration(A, M, beta,myfactor);

#if 0
		if (A.rowdim() < 200 ) {
			cra(d,iteration,genprime);
			commentator().stop ( "first step", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations() << "\n";
		}
		else {}
#endif
		Integer_t res;

#ifdef _LB_H_DET_TIMING
		Timer BT;
		double time1, time2;
		BT.start();
#endif
		{
			++genprime;
			mymodular D(*genprime);
			iteration.primes[early_counter] = *genprime;
			mymodular::Element r;
			D.assign(r,D.zero);
			iteration(r, D);
			cra.initialize( D, r);
			++early_counter;
		}

		while ( early_counter < myfactor && !cra.terminated() ) {
			++genprime;
			while (cra.noncoprime(*genprime)) ++genprime;
			mymodular D(*genprime);
			iteration.primes[early_counter] = *genprime;
			// prime(p, early_counter);
			mymodular::Element r;
			D.assign(r,D.zero);
			iteration(r,D);
			cra.progress( D, r);
			++early_counter;
		}

#ifdef _LB_H_DET_TIMING
		BT.stop();
		time1 = BT.usertime()/myfactor;
		if (time1 < 0) time1 = 0;
#endif
		cra.result(res);

		if (early_counter < myfactor) {
			/* determinant found */
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << myfactor << "\n";
			commentator().stop ( "first step", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations()<< "\n";
			return d=res;
		}
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << "no very early termination \n";
		/* turn to LU when matrix size small - not to be used at the moment */
#if 0
		if (A.rowdim() < 50 ) {
			while ( !cra.terminated() ) {
				genprime.randomPrime(p);
				mymodular D(p);
				iteration.primes[early_counter] = p;
				//		prime(p, early_counter);
				mymodular::Element r;
				D.assign(r,D.zero);
				iteration(r, D);
				cra.progress( D, r);
				++early_counter;
			}
			commentator().stop ( "zero step", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations() << "\n";
			return d;
		}

		PrimeIterator<IteratorCategories::HeuristicTag> genprime1(FieldTraits<mymodular>::bestBitSize(A.coldim()));

		Integers ZZ;
		DixonSolver < Integers , mymodular, PrimeIterator<IteratorCategories::HeuristicTag>, Method::DenseElimination > RSolver(A. field(), genprime);
#endif
		DixonSolver < Integers , mymodular, PrimeIterator<IteratorCategories::HeuristicTag>, Method::DenseElimination > RSolver;

		BlasVector<Integers> r_num1 (A.field(),A. coldim());

		LastInvariantFactor < Integers ,DixonSolver < Integers, mymodular, PrimeIterator<IteratorCategories::HeuristicTag>, Method::DenseElimination > >  LIF(RSolver);
#ifdef _LB_H_DET_TIMING
		BT.start();
#endif
		if (LIF.lastInvariantFactor1(lif, r_num1, A)==0) {
			//if (lif==0)
			d = 0;
			commentator().stop ("is 0", NULL, "det");
			return d;
		}
#ifdef _LB_H_DET_TIMING
		BT.stop();
		time2 = BT.usertime();
		if (time2 < 0) time2 =0;
#endif

		//commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		//            << "5 LU time: " << time1 << " LIF time: " << time2 << ")\n";
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << "lif calculated\n";
#if 0
		LIF.lastInvariantFactor_Bonus(lif, bonus, A);

		if (lif==0) {
			d = 0;
			commentator().stop ("done", NULL, "det");
			return d;
		}

		if (bonus == 1) {
			d = lif;
			commentator().stop ("done", NULL, "det");
			return d;
		}
#endif

		beta = lif*bonus;
		iteration.Beta(beta);

		CRABuilderEarlySingle<mymodular> cra2(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);
		Integer_t k = 1;

		early_counter = 0;
		while ( early_counter < myfactor && !cra2.terminated() ) {
			mymodular D(iteration.primes[early_counter]);
			mymodular::Element r;
			D.assign(r,D.zero);
			iteration(r,D);
			cra2.progress( D, r );
			++early_counter;
		}

		if (early_counter < myfactor) {
			/* determinant found */
			k = cra2.result(res);
			commentator().stop ("second step ", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";

		}
		else  if (0/* time2 < a*log2(lif)*time1/p_size*/) {
			BlasVector<Integers> r_num2 (A.field(),A. coldim());
			Integer_t lif2=1;
			LIF.lastInvariantFactor1(lif2,r_num2,A);
			LIF.bonus(bonus,lif, lif2, r_num1,r_num2);

			//lif2 = lif;
			if ((bonus > 1) || (lif2 != lif)) {
				//cout << "lif: " <<lif <<",\n     "<< lif2<< "przed\n";
				lif = lcm(lif, lif2);
				//cout << "lif: " <<lif << "po";
				beta=lif*bonus;
				iteration.Beta(beta);
				//iteration.Restart2() // included in Beta();
				//iteration.Moduli(moduli);
				//iteration.Primes(primes);
				k=1;
				CRABuilderEarlySingle<mymodular> cra3(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);

				early_counter = 0;
				while ( (early_counter < myfactor) && (!cra3.terminated() )) {
					mymodular D(iteration.primes[early_counter]);
					mymodular::Element r;
					D.assign(r,D.zero);
					iteration(r,D);
					cra3.progress( D, r);
					++early_counter;
					//iteration.Inc();
				}

				if (early_counter < myfactor) {
					/* determinant found based on the initial LU */
					cra3.result(k);
					commentator().stop ("third step - recalc", NULL, "det");
					commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
					<< "Iterations done " << iteration.iterations()<<"(" << iteration.iterations2() << ")\n";
					//<< "bonus size " << log2(bonus) << "\n";

				}
				else {
					/* enter the cra loop */
					//cra3(k,iteration, genprime);
					while (!cra3.terminated()) {
						++genprime;
						while (cra3.noncoprime(*genprime)) ++genprime;
						mymodular D(*genprime);
						mymodular::Element r;
						D.assign(r,D.zero);
						iteration(r,D);
						cra3.progress( D, r);
					}
					cra3.result(k);
					commentator().stop ("third step, bonus > 1", NULL, "det");
					commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
					<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
					//<< "bonus size " << log2(bonus) << "\n";
				}
			}
			else {
				//cra2(k,iteration, genprime);
				while (!cra2.terminated()) {
					++genprime;
					while (cra2.noncoprime(*genprime)) ++genprime;
					mymodular D(*genprime);
					mymodular::Element r;
					D.assign(r,D.zero);
					iteration(r,D);
					cra2.progress( D, r);
				}
				cra2.result(k);
				commentator().stop ("third step, bonus = 1", NULL, "det");
				commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
			}
		}
		else {
			while (!cra2.terminated()) {
				++genprime;
				while (cra2.noncoprime(*genprime)) ++genprime;
				mymodular D(*genprime);
				mymodular::Element r;
				D.assign(r,D.zero);
				iteration(r,D);
				cra2.progress( D, r);
			}
			cra2.result(k);
			commentator().stop ("second step+", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";

			/* enter the cra loop */
			//cra2(k,iteration, genprime);
		}

		//commentator().stop ("second step", NULL, "det");
		//commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		//        << "Iterations done " << iteration.iterations() << "("
		//        << iteration.iterations2() << " )\n";
		d = k*beta;

		Integer_t tmp;

#ifdef _LB_H_DET_TIMING
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "1 LU time: " << time1 << " LIF time: " << time2 << ")\n";
#endif
		commentator().report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "det/lif " << k<< "\n";
		//commentator().report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		//	<< "determinant size " << Integers::log2(tmp,abs(d))<< " lif size "<< Integers::log2(tmp,beta) << "\n";

		return d ;

	}

#if 0
	template <class Integers, class MyMethod>
	typename Integers::Element & lif_cra_det (typename Integers::Element                &d,
						  const SparseMatrix<Integers>              &A,
						  const RingCategories::IntegerTag          &tag,
						  const MyMethod                            &M)
	{

		//commentator().setReportStream(std::cout);
		typedef Givaro::ModularBalanced<double> mymodular;
		//typedef Givaro::ZRing<Integer> Integers;
		typedef typename Integers::Element Integer;

		commentator().start ("Integer Determinant - hybrid version for sparse matrices", "det");
		size_t myfactor=5;
		size_t early_counter=0;

		//double a = 0.0;//0.28;//0.0013i			//if time(lif)/time(10lu) < a * log(lif) then calculate bonus

		Integer lif = 1;
		Integer bonus = 1;
		Integer beta = 1;
		d=1;

                PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<mymodular>::bestBitSize(A.coldim()));

		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << "prime size: " << p_size << "\n";
		ChineseRemainder< mymodular > cra(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);
		IntegerModularDetReduced<SparseMatrix<Integers >,MyMethod> iteration(A, M, beta,myfactor);
#if 0
		if (A.rowdim() < 200 ) {
			cra(d,iteration,genprime);
			commentator().stop ( "first step", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations() << "\n";
		}
		else
#endif
			Integer p;
		Integer res;

		Timer BT;
		double time1, time2;
		BT.start();

                genprime.randomPrime(p);
                mymodular D(p);
                iteration.primes[early_counter] = p;
                mymodular::Element r;
                D.assign(r,D.zero);
		iteration(r,D);
                cra.initialize( D, r);
                ++early_counter;

		while ( early_counter < myfactor && !cra.terminated() ) {
			genprime.randomPrime(p);
			while (cra.noncoprime(p)) genprime.randomPrime(p);
			mymodular D(p);
			iteration.primes[early_counter] = p;
			//          prime(p, early_counter);
			mymodular::Element r;
			D.assign(r,D.zero);
			iteration(r,D);
			cra.progress( D,r);
			++early_counter;
		}

		BT.stop();
		time1 = BT.usertime()/myfactor;
		if (time1 < 0) time1 = 0;
		cra.result(res);

		if (early_counter < myfactor) {
			//determinant found
			commentator().stop ( "first step", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations()<< "\n";
			return d=res;
		}
#if 0
		cout << "no very early termination \n";
		turn to LU when matrix size small - not to be used at the moment
		if (A.rowdim() < 50 ) {
			while ( !cra.terminated() ) {
				genprime.randomPrime(p);
				mymodular D(p);
				iteration.primes[early_counter] = p;
				//          prime(p, early_counter);
				mymodular::Element r;
				D.assign(r,D.zero);
				iteration(r,D);
				cra.progress( D, r);
				++early_counter;
			}
			commentator().stop ( "zero step", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations() << "\n";
			return d;
		}

                PrimeIterator<IteratorCategories::HeuristicTag> genprime1(FieldTraits<mymodular>::bestBitSize(A.coldim()));
                Integers ZZ;
#endif
		typedef DixonSolver < Integers , mymodular, PrimeIterator<IteratorCategories::HeuristicTag>, Method::BlockHankel > Solver;
		Solver RSolver(A. field(), genprime);

		typename Vector<Integers>:: Dense r_num1 (A. coldim());

		LastInvariantFactor < Integers ,Solver >  LIF(RSolver);

		BT.start();
		if (LIF.lastInvariantFactor1(lif, r_num1, A)==0) {
			//if (lif==0)
			d = 0;
			commentator().stop ("is 0", NULL, "det");
			return d;
		}
		BT.stop();
		time2 = BT.usertime();
		if (time2 < 0) time2 =0;
#if 0
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "5 LU time: " << time1 << " LIF time: " << time2 << ")\n";
		cout << "lif calculated\n";
		LIF.lastInvariantFactor_Bonus(lif, bonus, A);

		if (lif==0) {
			d = 0;
			commentator().stop ("done", NULL, "det");
			return d;
		}

		if (bonus == 1) {
			d = lif;
			commentator().stop ("done", NULL, "det");
			return d;
		}
#endif
		beta = lif*bonus;
		iteration.Beta(beta);

		ChineseRemainder< Givaro::Modular<double> > cra2(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);
		Integer k = 1;

		early_counter = 0;
		while ( early_counter < myfactor && !cra2.terminated() ) {
			mymodular D(iteration.primes[early_counter]);
			mymodular::Element r;
			D.assign(r,D.zero);
			iteration(r,D);
			cra2.progress( D, r);
			++early_counter;
		}

		if (early_counter < myfactor) {
			// determinant found
			k = cra2.result(res);
			commentator().stop ("second step ", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";

		}
		else  if (0) { //time2 < a*log2(lif)*time1/p_size)
			typename Vector<Integers>:: Dense r_num2 (A. coldim());
			Integer lif2=1;
			LIF.lastInvariantFactor1(lif2,r_num2,A);
			LIF.bonus(bonus,lif, lif2, r_num1,r_num2);

			//lif2 = lif;
			if ((bonus > 1) || (lif2 != lif)) {
				//cout << "lif: " <<lif <<",\n     "<< lif2<< "przed\n";
				lif = lcm(lif, lif2);
				//cout << "lif: " <<lif << "po";
				beta=lif*bonus;
				iteration.Beta(beta);
				//iteration.Restart2() // included in Beta();
				//iteration.Moduli(moduli);
				//iteration.Primes(primes);
				k=1;
				ChineseRemainder< Givaro::Modular<double> > cra3(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);

				early_counter = 0;
				while ( (early_counter < myfactor) && (!cra3.terminated() )) {
					mymodular D(iteration.primes[early_counter]);
					mymodular::Element r;
					D.assign(r,D.zero);
					iteration(r,D);
					cra3.progress( D, r);
					++early_counter;
					//iteration.Inc();
				}

				if (early_counter < myfactor) {
					// determinant found based on the initial LU
					cra3.result(k);
					commentator().stop ("third step - recalc", NULL, "det");
					commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
					<< "Iterations done " << iteration.iterations()<<"(" << iteration.iterations2() << ")\n";
					//<< "bonus size " << log2(bonus) << "\n";

				}
				else {
					// enter the cra loop
					//cra3(k,iteration, genprime);
					while (!cra3.terminated()) {
						genprime.randomPrime(p);
						while (cra3.noncoprime(p)) genprime.randomPrime(p);
						mymodular D(p);
						mymodular::Element r;
						D.assign(r,D.zero);
						iteration(r,D);
						cra3.progress( D, r);
					}
					cra3.result(k);
					commentator().stop ("third step, bonus > 1", NULL, "det");
					commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
					<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
					//<< "bonus size " << log2(bonus) << "\n";
				}
			}
			else {
				//cra2(k,iteration, genprime);
				while (!cra2.terminated()) {
					genprime.randomPrime(p);
					while (cra2.noncoprime(p)) genprime.randomPrime(p);
					mymodular D(p);
					mymodular::Element r;
					D.assign(r,D.zero);
					iteration(r,D);
					cra2.progress( D, r);
				}
				cra2.result(k);
				commentator().stop ("third step, bonus = 1", NULL, "det");
				commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";
			}
		}
		else {
			while (!cra2.terminated()) {
				genprime.randomPrime(p);
				while (cra2.noncoprime(p)) genprime.randomPrime(p);
				mymodular D(p);
				mymodular::Element r;
				D.assign(r,D.zero);
				iteration(r,D);
				cra2.progress( D, r);
			}
			cra2.result(k);
			commentator().stop ("second step+", NULL, "det");
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Iterations done " << iteration.iterations()<< "(" << iteration.iterations2() << ")\n";

			// enter the cra loop
			//cra2(k,iteration, genprime);
		}

		//commentator().stop ("second step", NULL, "det");
		//commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		//        << "Iterations done " << iteration.iterations() << "("
		//        << iteration.iterations2() << " )\n";
		d = k*beta;

		Integer tmp;

		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "1 LU time: " << time1 << " LIF time: " << time2 << ")\n";
		commentator().report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "det/lif " << k<< "\n";
		//commentator().report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		//        << "determinant size " << Integers::log2(tmp,abs(d))<< " lif size "<< Integers::log2(tmp,beta) << "\n";

		return d ;

	}
#endif

} // end of LinBox namespace

#undef _LB_H_DET_TIMING

#endif // __LINBOX_hybrid_det_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
