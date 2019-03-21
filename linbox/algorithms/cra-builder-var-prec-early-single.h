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
#include "linbox/algorithms/cra-builder-single.h"
#include "linbox/algorithms/cra-builder-full-multip.h"
#include "linbox/algorithms/lazy-product.h"


namespace LinBox
{

	/*
	 * Implements early terminated CRA with preconditioning of the result
	 * factor | result should be given
	 * does not detect bad factors but can run until full termination [for result] in this case
	 * factor may be changed by changeFactor function, residues stored in CRABuilderFullMultip are recomputed in this case,
	 * FullMulpitCRA should consist of vectors of size 0, but errors happen.
	 */
	template<class Domain_Type>
	struct CRABuilderVarPrecEarlySingle: public CRABuilderEarlySingle<Domain_Type>, CRABuilderFullMultip<Domain_Type> {

		//typedef Givaro::ZRing<Integer> Integers;
		//typedef Integers::Element Integer;
		typedef GMPRationalField Rationals;
		typedef Rationals::Element Quotient;

		typedef Domain_Type                     Domain;
		typedef typename Domain::Element DomainElement;
		typedef CRABuilderVarPrecEarlySingle<Domain> Self_t;

	public:
		Integer factor_;
		Integer multip_;

		CRABuilderVarPrecEarlySingle(const size_t EARLY = LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD
				      , const double b=0.0
				      , const Integer& f=Integer(1)
				      , const Integer& m=Integer(1)) :
			CRABuilderEarlySingle<Domain>(EARLY)
			, CRABuilderFullMultip<Domain>(b)
			, factor_(f)
			, multip_(m)
		{
			if (factor_ == 0 )
			       	factor_ = 1;
		}

		CRABuilderVarPrecEarlySingle(const CRABuilderVarPrecEarlySingle& other) :
			CRABuilderEarlySingle<Domain>(other.EARLY_TERM_THRESHOLD)
			, CRABuilderFullMultip<Domain>(other.LOGARITHMIC_UPPER_BOUND)
			, factor_(other.factor_), multip_(other.multip_)
		{
			factor_ = 1;
		}

		int getThreshold(int& t)
		{
			return t = (int)CRABuilderEarlySingle<Domain>::EARLY_TERM_THRESHOLD;
		}

		Integer& getModulus(Integer& m)
		{
			CRABuilderEarlySingle<Domain>::getModulus(m);
			return m;
		}

		/*
		 * Residue is Result * M / F - this is the value that we are reconstructing
		 */
		Integer& getResidue(Integer& m)
		{
			CRABuilderEarlySingle<Domain>::getResidue(m);
			return m;
		}

		void initialize (const Integer& D, const Integer& e)
		{
			Integer z;
			inv(z,factor_,D);
			z*=e;
			z*=multip_;
			z%=D;

			CRABuilderEarlySingle<Domain>::initialize(D, z);
			Givaro::ZRing<Integer> ZZ ;
			BlasVector<Givaro::ZRing<Integer> > v(ZZ);
			v.push_back(e);
			CRABuilderFullMultip<Domain>::initialize(D, v);
		}

		void initialize (const Domain& D, const DomainElement& e)
		{
			DomainElement z;
			D.init(z,factor_);
			D.invin(z);
			D.mulin(z,e);
			DomainElement m; D.init(m, multip_);
			D.mulin(z,m);

			CRABuilderEarlySingle<Domain>::initialize(D, z);
			BlasVector<Domain> v(D);
			v.push_back(e);
			CRABuilderFullMultip<Domain>::initialize(D, v);
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
			CRABuilderEarlySingle<Domain>::progress(D, z);

			Givaro::ZRing<Integer> ZZ ;
			BlasVector<Givaro::ZRing<Integer> > v(ZZ);
			v.push_back(e);
			CRABuilderFullMultip<Domain>::progress(D, v);
		}

		void progress (const Domain& D, const DomainElement& e)
		{

			DomainElement z;
			D.init(z,factor_);
			D.invin(z);
			D.mulin(z,e);
			DomainElement m; D.init(m, multip_);
			D.mulin(z,m);
			//z = (e/ factor mod D)

			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.

			CRABuilderEarlySingle<Domain>::progress(D, z);

			BlasVector<Domain> v(D);
			v.push_back(e);

			CRABuilderFullMultip<Domain>::progress(D, v);
		}

		bool terminated()
		{
			bool ET = CRABuilderEarlySingle<Domain>::terminated();
			if (CRABuilderFullMultip<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) ET = ET || CRABuilderFullMultip<Domain>::terminated();
			return ET;
		}

		bool noncoprime(const Integer& i) const
		{
			return CRABuilderEarlySingle<Domain>::noncoprime(i);
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
			if ((CRABuilderFullMultip<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( CRABuilderFullMultip<Domain>::terminated() )) {
				Givaro::ZRing<Integer> ZZ ;
				BlasVector<Givaro::ZRing<Integer> > v(ZZ);
				CRABuilderFullMultip<Domain>::result(v);
				return r = v.front();
			}
			else {
				Integer z;
				CRABuilderEarlySingle<Domain>::result(z);
				return (r=factor_*z/multip_);
			}
		}

		Quotient& result(Quotient& q)
		{
			if ((CRABuilderFullMultip<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( CRABuilderFullMultip<Domain>::terminated() )) {
				Givaro::ZRing<Integer> ZZ ;
				BlasVector<Givaro::ZRing<Integer> > v(ZZ);
				CRABuilderFullMultip<Domain>::result(v);
				return q = v.front();
			}
			else {
				Integer z;
				CRABuilderEarlySingle<Domain>::result(z);//residue
				return (q=Quotient(factor_*z,multip_));

			}
		}

		Integer& result(Integer& num, Integer& den)
		{
			if ((CRABuilderFullMultip<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( CRABuilderFullMultip<Domain>::terminated() )) {
				den =1;
				Givaro::ZRing<Integer> ZZ ;
				BlasVector<Givaro::ZRing<Integer> > v(ZZ);
				CRABuilderFullMultip<Domain>::result(v);
				return num = v.front();
			}
			else {
				Integer z;
				CRABuilderEarlySingle<Domain>::result(z);
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
			if ((factor_ == f) && (multip_==m)) return CRABuilderEarlySingle<Domain>::terminated();

			factor_ = f;
			multip_ = m;
			if (factor_==0) {factor_ = 1;}//no factor if factor==0

			Integer e=0;

			//clear CRAEarlySingle;
			CRABuilderEarlySingle<Domain>::occurency_ = 0;
			CRABuilderEarlySingle<Domain>::nextM_ = 1;
			CRABuilderEarlySingle<Domain>::primeProd_ = 1;
			CRABuilderEarlySingle<Domain>::residue_ = 0;

			//Computation of residue_
            for (auto it = CRABuilderFullMultip<Domain>::shelves_begin();
                 it != CRABuilderFullMultip<Domain>::shelves_end();
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


					Integer prev_residue_ = CRABuilderEarlySingle<Domain>::residue_;
					CRABuilderEarlySingle<Domain>::progress(D,e_i);

					if (prev_residue_ == CRABuilderEarlySingle<Domain>::residue_ ) {
						CRABuilderEarlySingle<Domain>::occurency_ += it->count;
					}



					if ( CRABuilderEarlySingle<Domain>::terminated() ) {
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
