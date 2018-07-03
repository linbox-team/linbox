/* linbox/blackbox/rational-reconstruction-base.h
 * Copyright (C) 2009 Anna Marszalek
 *
 * Written by Anna Marszalek <aniau@astronet.pl>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_varprec_cra_early_single_H
#define __LINBOX_varprec_cra_early_single_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/field/gmp-rational.h"
#include "linbox/solutions/methods.h"
#include <utility>
#include "linbox/algorithms/cra-single.h"
#include "linbox/algorithms/cra-full-multip.h"
#include "linbox/algorithms/lazy-product.h"


namespace LinBox
{

	/*
	 * Implements early terminated CRA with preconditioning of the result
	 * factor | result should be given
	 * does not detect bad factors but can run until full termination [for result] in this case
	 * factor may be changed by changeFactor function, residues stored in FullMultipCRA are recomputed in this case,
	 * FullMulpitCRA should consist of vectors of size 0, but errors happen.
	 */
	struct VarPrecEarlySingleCRA: public EarlySingleCRA, FullMultipCRA {

		//typedef Givaro::ZRing<Integer> Integers;
		//typedef Integers::Element Integer;
		typedef GMPRationalField Rationals;
		typedef Rationals::Element Quotient;

		typedef VarPrecEarlySingleCRA Self_t;
        using SingleParent = EarlySingleCRA;
        using MultiParent = FullMultipCRA;

	public:
		Integer factor_;
		Integer multip_;

		VarPrecEarlySingleCRA(const unsigned long EARLY = DEFAULT_EARLY_TERM_THRESHOLD
				      , const double b=0.0
				      , const Integer& f=Integer(1)
				      , const Integer& m=Integer(1)) :
			SingleParent(EARLY)
			, MultiParent(b)
			, factor_(f)
			, multip_(m)
		{
			if (factor_ == 0 )
			       	factor_ = 1;
		}

		VarPrecEarlySingleCRA(const VarPrecEarlySingleCRA& other) :
			SingleParent(other.EARLY_TERM_THRESHOLD)
			, MultiParent(other.LOGARITHMIC_UPPER_BOUND)
			, factor_(other.factor_), multip_(other.multip_)
		{
			factor_ = 1;
		}

		int getThreshold(int& t)
		{
			return t = (int)SingleParent::EARLY_TERM_THRESHOLD;
		}

		Integer& getModulus(Integer& m)
		{
			SingleParent::getModulus(m);
			return m;
		}

		/*
		 * Residue is Result * M / F - this is the value that we are reconstructing
		 */
		Integer& getResidue(Integer& m)
		{
			SingleParent::getResidue(m);
			return m;
		}

		void initialize (const Integer& D, const Integer& e)
		{
			Integer z;
			inv(z,factor_,D);
			z*=e;
			z*=multip_;
			z%=D;

			SingleParent::initialize(D, z);
			Givaro::ZRing<Integer> ZZ ;
			BlasVector<Givaro::ZRing<Integer> > v(ZZ);
			v.push_back(e);
			MultiParent::initialize(D, v);
		}

        template <class Domain>
		void initialize (const Domain& D, const typename Domain::Element& e)
		{
			typename Domain::Element z;
			D.init(z,factor_);
			D.invin(z);
			D.mulin(z,e);
			typename Domain::Element m; D.init(m, multip_);
			D.mulin(z,m);

			SingleParent::initialize(D, z);
			BlasVector<Domain> v(D);
			v.push_back(e);
			MultiParent::initialize(D, v);
		}

		void progress (const Integer& D,const  Integer& e)
		{

			Integer z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.

			inv(z,factor_,D);
			z*=e;
			z*=multip_;
			z%=D;
			//z = e / factor mod D;
			SingleParent::progress(D, z);

			Givaro::ZRing<Integer> ZZ ;
			BlasVector<Givaro::ZRing<Integer> > v(ZZ);
			v.push_back(e);
			MultiParent::progress(D, v);
		}

        template <class Domain>
		void progress (const Domain& D, const typename Domain::Element& e)
		{

			typename Domain::Element z;
			D.init(z,factor_);
			D.invin(z);
			D.mulin(z,e);
			typename Domain::Element m; D.init(m, multip_);
			D.mulin(z,m);
			//z = (e/ factor mod D)

			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.

			SingleParent::progress(D, z);

			BlasVector<Domain> v(D);
			v.push_back(e);

			MultiParent::progress(D, v);
		}

		bool terminated()
		{
			bool ET = SingleParent::terminated();
			if (MultiParent::LOGARITHMIC_UPPER_BOUND> 1.0) ET = ET || MultiParent::terminated();
			return ET;
		}

		bool noncoprime(const Integer& i) const
		{
			return SingleParent::noncoprime(i);
		}

		Integer& getFactor(Integer& f)
		{
			return f = factor_;
		}

		Integer& getMultip(Integer& m)
		{
			return m=multip_;
		}

		Integer& getPreconditioner(Integer& f, Integer& m)
		{
			getMultip(m);
			getFactor(f);
			return f;
		}

		Quotient& getPreconditioner(Quotient& q)
		{
			Integer m,f;
			getPreconditioner(f,m);
			return q=Quotient(m,f);
		}

		Integer& result(Integer& r)
		{
			if ((MultiParent::LOGARITHMIC_UPPER_BOUND> 1.0) && ( MultiParent::terminated() )) {
				Givaro::ZRing<Integer> ZZ ;
				BlasVector<Givaro::ZRing<Integer> > v(ZZ);
				MultiParent::result(v);
				return r = v.front();
			}
			else {
				Integer z;
				SingleParent::result(z);
				return (r=factor_*z/multip_);
			}
		}

		Quotient& result(Quotient& q)
		{
			if ((MultiParent::LOGARITHMIC_UPPER_BOUND> 1.0) && ( MultiParent::terminated() )) {
				Givaro::ZRing<Integer> ZZ ;
				BlasVector<Givaro::ZRing<Integer> > v(ZZ);
				MultiParent::result(v);
				return q = v.front();
			}
			else {
				Integer z;
				SingleParent::result(z);//residue
				return (q=Quotient(factor_*z,multip_));

			}
		}

		Integer& result(Integer& num, Integer& den)
		{
			if ((MultiParent::LOGARITHMIC_UPPER_BOUND> 1.0) && ( MultiParent::terminated() )) {
				den =1;
				Givaro::ZRing<Integer> ZZ ;
				BlasVector<Givaro::ZRing<Integer> > v(ZZ);
				MultiParent::result(v);
				return num = v.front();
			}
			else {
				Integer z;
				SingleParent::result(z);
				den = multip_;
				num = factor_*z;
				return num;
			}
		}

		/*
		 * changes preconditioners and recomputes the residue
		 */
		bool changePreconditioner(const Integer& f, const Integer& m)
		{
			if ((factor_ == f) && (multip_==m)) return SingleParent::terminated();

			factor_ = f;
			multip_ = m;
			if (factor_==0) {factor_ = 1;}//no factor if factor==0

			Integer e=0;

			//clear CRAEarlySingle;
			SingleParent::occurency_ = 0;
			SingleParent::nextM_ = 1;
			SingleParent::primeProd_ = 1;
			SingleParent::residue_ = 0;

			//Computation of residue_
            for (auto it = MultiParent::shelves_begin();
                 it != MultiParent::shelves_end();
                 ++it)
            {
                if (it->occupied) {
					Integer D = it->mod();
					Integer e_i;

					inv(e_i,factor_,D);
					//!@todo use faster mul/mod here !
					Integer::mulin(e_i, (it->residue.front()));
					Integer::mulin(e_i, multip_);
					e_i %=D ;


					Integer prev_residue_ = SingleParent::residue_;
					SingleParent::progress(D,e_i);

					if (prev_residue_ == SingleParent::residue_ ) {
						SingleParent::occurency_ += it->count;
					}



					if ( SingleParent::terminated() ) {
						return true;
					}
                }
            }

			return false;
		}

	};

} //namespace LinBox

#endif //__LINBOX_varprec_cra_early_single_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
