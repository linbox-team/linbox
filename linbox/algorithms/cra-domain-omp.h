/* linbox/algorithms/cra-domain-omp.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Naive parallel chinese remaindering
 * Launch NN iterations in parallel, where NN=omp_get_max_threads()
 * Then synchronization and termintation test.
 * Time-stamp: <13 Mar 12 13:49:58 Jean-Guillaume.Dumas@imag.fr>
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

/*! @file algorithms/cra-domain-omp.h
 * @brief Parallel (OMP) version of \ref CRA
 * @ingroup CRA
 */

#ifndef __LINBOX_omp_cra_H
#define __LINBOX_omp_cra_H
// commentator is not thread safe
#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif
#include <omp.h>
#include <set>
#include "linbox/algorithms/cra-domain-sequential.h"

namespace LinBox
{

	template<class CRABase>
	struct ChineseRemainderOMP : public ChineseRemainderSequential<CRABase> {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
		typedef ChineseRemainderSequential<CRABase>    Father_t;

		template<class Param>
		ChineseRemainderOMP(const Param& b) :
			Father_t(b)
		{}

		ChineseRemainderOMP(const CRABase& b) :
			Father_t(b)
		{}

		template <class ResultType, class Function, class PrimeIterator>
		ResultType& operator() (ResultType& res, Function& Iteration, PrimeIterator& primeiter)
		{
			using ResidueType = typename CRAResidue<ResultType,Function>::template ResidueType<Domain>;
			size_t NN = omp_get_max_threads();
			//std::cerr << "Blocs: " << NN << " iterations." << std::endl;
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);

			std::vector<Domain> ROUNDdomains; ROUNDdomains.reserve(NN);
			std::vector<ResidueType> ROUNDresidues; ROUNDresidues.reserve(NN);
			std::vector<IterationResult> ROUNDresults(NN);
			std::set<Integer> coprimeset;

			while (! this->Builder_.terminated()) {
				ROUNDdomains.clear();
				ROUNDresidues.clear();
				coprimeset.clear();

				while (coprimeset.size() < NN) {
					coprimeset.emplace(this->get_coprime(primeiter));
					++primeiter;
				}

				for(auto coprimesetiter = coprimeset.cbegin(); coprimesetiter != coprimeset.cend(); ++coprimesetiter) {
					// std::cerr << "With prime: " << *coprimesetiter << std::endl;
					ROUNDdomains.emplace_back(*coprimesetiter);
					ROUNDresidues.emplace_back(CRAResidue<ResultType,Function>::create(ROUNDdomains.back()));
				}

#pragma omp parallel for
				for(size_t i=0;i<NN;++i) {
					ROUNDresults[i] = Iteration(ROUNDresidues[i], ROUNDdomains[i]);
				}
#pragma omp barrier

				// if any thread says RESTART, then all CONTINUEs become SKIPs
				bool anyrestart = false;
				for (auto res : ROUNDresults) {
					if (res == IterationResult::RESTART) anyrestart = true;
				}
				if (anyrestart) {
					this->nbad_ += this->ngood_;
					this->ngood_ = 0;
				}

				for (size_t i = 0; i < NN; ++i) {
					if (ROUNDresults[i] == IterationResult::SKIP) {
						this->doskip();
					}
					else if (anyrestart && ROUNDresults[i] == IterationResult::CONTINUE) {
						// commentator should indicate that this prime is bad
						++this->nbad_;
					}
					else if (this->ngood_ == 0) {
						this->ngood_ = 1;
						this->Builder_.initialize(ROUNDdomains[i], ROUNDresidues[i]);
					}
					else {
						++this->ngood_;
						this->Builder_.progress(ROUNDdomains[i], ROUNDresidues[i]);
					}
				}
				//std::cerr << "Computed: " << iterCount() << " primes." << std::endl;
			}

			// commentator().stop ("done", NULL, "mmcrait");
			//std::cerr << "Used: " << this->iterCount() << " primes." << std::endl;
			return this->Builder_.result(res);
		}
	};
}

#endif //__LINBOX_omp_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
