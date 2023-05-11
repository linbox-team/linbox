/* linbox/algorithms/cra-domain-parallel.h
 * Copyright (C) 1999-2022 The LinBox group
 *
 * Time-stamp: <28 Oct 22 18:35:25 Jean-Guillaume.Dumas@imag.fr>
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

/*! @file algorithms/cra-domain-parallel.h
 * @brief Parallel (PALADIN) version of \ref CRA
 * @brief Naive parallel chinese remaindering
 * @brief Launch by blocks of NN iterations
 * @brief in parallel, by default, NN is numver of available threads
 * @brief Then synchronization and termintation test.
 * @ingroup CRA
 */

#ifndef __LINBOX_parallel_cra_H
#define __LINBOX_parallel_cra_H

#  ifndef DISABLE_COMMENTATOR
#    warning "commentator is not thread safe"
#    define DISABLE_COMMENTATOR
#  endif

#include "linbox/algorithms/cra-domain-sequential.h"

#ifndef __LB_CRA_REPORTING__
# ifdef _LB_DEBUG
#  define __LB_CRA_REPORTING__ 1
# else
#  define __LB_CRA_REPORTING__ 0
# endif
#endif

namespace LinBox
{

	template<class CRABase>
	struct ChineseRemainderParallel
        : public ChineseRemainderSequential<CRABase> {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
		typedef ChineseRemainderSequential<CRABase>    Father_t;
        typedef ChineseRemainderParallel<CRABase>    Self_t;

        friend std::ostream& operator<< (std::ostream& out, const Self_t& cra) {
            std::ostringstream report;
            report << "Parallel Chinese Remaindering on " << NUM_THREADS << " threads out of " << MAX_THREADS << std::endl;
            return out << report.str();
        }

		template<class Param>
		ChineseRemainderParallel(const Param& b) :
			Father_t(b)
		{
#if __LB_CRA_REPORTING__
            std::clog << *this << std::endl;
#endif
        }

		ChineseRemainderParallel(const CRABase& b) :
			Father_t(b)
		{
#if __LB_CRA_REPORTING__
            std::clog << *this << std::endl;
#endif
        }

		template <class ResultType, class Function, class PrimeIterator>
		ResultType& operator() (ResultType& res, Function& Iteration, PrimeIterator& primeiter)
		{
            (*this)(-1, res, Iteration, primeiter);
            return res;
        }

		template <class ResultType, class Function, class PrimeIterator>
		bool operator() (int k, ResultType& res, Function& Iteration, PrimeIterator& primeiter, size_t NN = NUM_THREADS)
        {
//             std::clog << "Parallel Givaro::Modular iteration, blocks " << NN << " iterations." << std::endl;
			using ResidueType = typename CRAResidue<ResultType,Function>::template ResidueType<Domain>;
			if (NN == 1) return Father_t::operator()(k, res,Iteration,primeiter);

			std::vector<Domain> ROUNDdomains; ROUNDdomains.reserve(NN);
			std::vector<ResidueType> ROUNDresidues; ROUNDresidues.reserve(NN);
			std::vector<IterationResult> ROUNDresults(NN);
			std::set<Integer> coprimeset;

			while (k != 0 && ! this->Builder_.terminated()) {
                if ( (k>0) && (size_t(k)<NN) ) NN = k;
                k -= NN;
				ROUNDdomains.clear();
				ROUNDresidues.clear();
				coprimeset.clear();

				while (coprimeset.size() < NN) {
					coprimeset.emplace(this->get_coprime(primeiter));
					++primeiter;
				}

				for(auto coprimesetiter = coprimeset.cbegin(); coprimesetiter != coprimeset.cend(); ++coprimesetiter) {
					ROUNDdomains.emplace_back(*coprimesetiter);
					ROUNDresidues.emplace_back(CRAResidue<ResultType,Function>::create(ROUNDdomains.back()));
				}


                SYNCH_GROUP(
                for(size_t i=0;i<NN;++i) {
                { TASK(MODE(CONSTREFERENCE(ROUNDdomains,ROUNDresidues,ROUNDresults)
                            WRITE(ROUNDresults[i], ROUNDresidues[i]) ),
                {
#if __LB_CRA_REPORTING__
                    std::ostringstream report;
                    ROUNDdomains[i].write(report <<
                                          "Iteration lauch " << i <<
                                          " on T" << THREAD_INDEX <<
                                          " over ") << std::endl;
                    std::clog << report.str();
#endif

                    ROUNDresults[i] = Iteration(ROUNDresidues[i], ROUNDdomains[i]);

                })}
                }
                )


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
					else if (anyrestart && ROUNDresults[i] == IterationResult::SKIP) {
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

#if __LB_CRA_REPORTING__
                std::clog << "Current good/bad residues: "
                          << this->ngood_ << '/'
                          << this->nbad_ << std::endl;
#endif

			}

			this->Builder_.result(res);
			return this->Builder_.terminated();
		}
	};
}

#endif //__LINBOX_parallel_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
