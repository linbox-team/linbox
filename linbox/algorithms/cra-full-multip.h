/* linbox/algorithms/cra-full-multip.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Time-stamp: <15 Dec 10 15:54:00 Jean-Guillaume.Dumas@imag.fr>
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

/*!@file algorithms/cra-full-multip.h
 * @ingroup algorithms
 * @brief NO DOC
 */

#ifndef __LINBOX_cra_full_multip_H
#define __LINBOX_cra_full_multip_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include "linbox/vector/blas-vector.h"
#include <utility>

#include "linbox/algorithms/lazy-product.h"

namespace LinBox
{

	/*! NO DOC...
	 * @ingroup CRA
	 * @bib
	 * - Jean-Guillaume Dumas, Thierry Gautier et Jean-Louis Roch.  <i>Generic design
	 * of Chinese remaindering schemes</i>  PASCO 2010, pp 26-34, 21-23 juillet,
	 * Grenoble, France.
	 */
	template<class Domain_Type>
	struct FullMultipCRA {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef FullMultipCRA<Domain> 		Self_t;

    public:
        struct Shelf {
            bool occupied = false;
            std::vector<Integer> residue;
            LazyProduct mod;
            double logmod = 0.; // natural log of the modulus on this shelf
            int count = 0; // how many images combined to get this one

            Shelf(size_t dim=0) :residue(dim) { };
        };

	protected:
        std::vector<Shelf> shelves_;
		const double				LOGARITHMIC_UPPER_BOUND;
		double totalsize_ = 0.; // natural log of the current modulus
        size_t dimension_ = 0; // dimension of the vector being reconstructed
        bool collapsed_ = false;
        bool normalized_ = false;
        // INVARIANT: shelves_.empty() || shelves_.back().occupied
        // INVARIANT: forall (shelf : shelves_) { shelf.residue.size() == dimension_ }

	public:
        /** @brief Creates a new vector CRA object.
         * @param bnd  upper bound on the natural logarithm of the result
         * @param dim  dimension of the vector to be reconstructed
         */
		FullMultipCRA(const double bnd=0.0, size_t dim=0) :
			LOGARITHMIC_UPPER_BOUND(bnd), dimension_(dim)
		{}

		Integer& getModulus(Integer& m) const
		{
            if (shelves_.empty()) return m = 1;
            collapse();
            return m = shelves_.back().mod();
		}

        const Integer& getModulus() const {
            collapse();
            return shelves_.back().mod();
        }

        // alias for result
		template<class Vect>
		inline Vect& getResidue(Vect& r) const
		{
			return result(r);
		}

		//! init
		template<typename ModType, class Vect>
		inline void initialize (const ModType& D, const Vect& e)
		{
            initialize_iter(D, e.begin(), e.size());
		}

        template <typename ModType, class Iter>
        inline void initialize_iter (const ModType& D, Iter e_it, size_t e_size)
        {
            shelves_.clear();
            totalsize_ = 0;
            dimension_ = e_size;
            progress_iter(D, e_it, e_size);
        }

        // generic version
		template <typename ModType, class Vect>
		inline void progress (const ModType& D, const Vect& e)
		{
            // resize existing residues if necessary
            if (e.size() > dimension_) {
                dimension_ = e.size();
                for (auto& shelf : shelves_) {
                    shelf.residue.resize(dimension_);
                }
            }

            // call iterator version
            progress_iter(D, e.begin(), e.size());
        }

        template <typename ModType, class Iter>
        void progress_iter (const ModType& D, Iter e_it, size_t e_size) {
            // update collapsed_ and normalized_
            collapsed_ = shelves_.empty();
            normalized_ = false;

            // put new result into the proper shelf
            Integer Dval;
            getMod(Dval, D);
            double logD = Givaro::naturallog(Dval);
            auto cur = getShelf(logD);

            totalsize_ += logD;

            ensureShelf(cur, shelves_, dimension_);
            if (! shelves_[cur].occupied) {
                // shelf is empty, so just copy it there
                std::copy_n(e_it, e_size, shelves_[cur].residue.begin());
                shelves_[cur].mod.initialize(Dval);
                shelves_[cur].logmod = logD;
                shelves_[cur].count = 1;
                shelves_[cur].occupied = true;
                return;
            }

            Integer invprod;

            // shelf is nonempty, so we incorporate the new result there
            {
                // FIXME incorporate invP0 and fieldreconstruct when D is a domain
                precomputeInvProd(invprod, shelves_[cur].mod(), Dval);
                auto r_it = shelves_[cur].residue.begin();
                for (size_t i=0; i < e_size; ++i, ++e_it, ++r_it) {
                    smallbigreconstruct(*r_it, *e_it, invprod);
                }
                // in case e is shorter than dimension_, treat missing values as zeros
                for (; r_it != shelves_[cur].residue.end(); ++r_it) {
                    *r_it *= invprod;
                }
                shelves_[cur].mod.mulin(Dval);
                shelves_[cur].logmod += logD;
                shelves_[cur].count += 1;
            }

            // combine further shelves as necessary
            decltype(cur) next;
            while ((next = getShelf(shelves_[cur].logmod)) != cur) {
                ensureShelf(next, shelves_, dimension_);
                if (shelves_[next].occupied) {
                    // combine cur shelf with next shelf
                    combineShelves(shelves_[next], shelves_[cur], invprod);
                    shelves_[cur].occupied = false;
                } else {
                    // put cur shelf data in next shelf position
                    std::swap(shelves_[cur], shelves_[next]);
                }

                cur = next;
            }
		}

		//! result
		inline const std::vector<Integer>& result () const
		{
            normalize();
            return shelves_.back().residue;
        }

        template <class Vect>
        inline Vect& result(Vect& r) const
        {
            r.resize(dimension_);
            result_iter(r.begin());
            return r;
        }

        template <class Iter>
        void result_iter (Iter r_it) const {
            if (shelves_.empty()) {
                for (size_t i=0; i < dimension_; ++i)
                    *r_it = 0;
            }
            else {
                normalize();
                std::copy_n(shelves_.back().residue.begin(), dimension_, r_it);
            }
        }

		bool terminated() const
		{
			return totalsize_ > LOGARITHMIC_UPPER_BOUND;
		}

		bool noncoprime(const Integer& i) const
		{
            for (auto& shelf : shelves_) {
                if (shelf.occupied && shelf.mod.noncoprime(i)) return true;
            }
            return false;
		}

        // XXX iterator invalidated by many other method calls
        decltype(shelves_.crbegin()) shelves_begin() const {
            return shelves_.rbegin();
        }

        decltype(shelves_.crend()) shelves_end() const {
            return shelves_.rend();
        }

	protected:
        /** Returns the index where the shelf (with specified natural log of modulus) belongs.
         */
        static inline size_t getShelf(double logmod) {
            // note, the (-3) is because the smallest "reasonable" modulus is a 15-bit or 23-bit prime
            return std::max(std::ilogb(logmod), 3) - 3;
        }

        static inline Integer& getMod(Integer& m, const Integer& D) {
            return m = D;
        }

        template <class Domain>
        static inline Integer& getMod(Integer& m, const Domain& D) {
            return D.characteristic(m);
        }

        static inline void combineShelves(Shelf& dest, const Shelf& src, Integer& temp) {
            // assumption: dest is already occupied
            precomputeInvProd(temp, dest.mod(), src.mod());
            for (auto dest_it = dest.residue.begin(), src_it = src.residue.begin();
                 dest_it != dest.residue.end();
                 ++dest_it, ++src_it)
            {
                smallbigreconstruct(*dest_it, *src_it, temp);
            }
            dest.mod.mulin(src.mod());
            dest.logmod += src.logmod;
            dest.count += src.count;
        }

        static inline void ensureShelf(size_t index, std::vector<Shelf>& shelves, size_t dim) {
            while (index >= shelves.size()) {
                shelves.emplace_back(dim);
            }
        }

        // collapse all residues to a single shelf
        void collapse() const {
            if (collapsed_) return;
            auto& ncshelves = const_cast<std::vector<Shelf>&>(shelves_);
            if (ncshelves.empty()) {
                ncshelves.emplace_back(dimension_);
                ncshelves.back().occupied = true;
            }
            else {
                Integer temp;
                auto cur = ncshelves.begin();
                while(true) {
                    auto next = cur;
                    if (++next == ncshelves.end()) break;
                    while (! next->occupied) ++next;
                    combineShelves(*next, *cur, temp);
                    cur->occupied = false;
                    cur = next;
                }
                // move top shelf to higher index if necessary
                auto newtop = getShelf(ncshelves.back().logmod);
                if (newtop >= ncshelves.size()) {
                    auto curtop = ncshelves.size()-1;
                    ensureShelf(newtop, ncshelves, dimension_);
                    std::swap(ncshelves[curtop], ncshelves[newtop]);
                }
            }
            const_cast<bool&>(collapsed_) = true;
        }

        // collapse and then normalize the top shelf
        void normalize() const {
            if (normalized_) return;
            collapse();
            Integer halfm = shelves_.back().mod();
            --halfm;
            halfm >>= 1;
            for (auto& x : const_cast<std::vector<Integer>&>(shelves_.back().residue)) {
                Integer::modin(x, shelves_.back().mod());
                if (x > halfm) x -= shelves_.back().mod();
            }
            const_cast<bool&>(normalized_) = true;
        }

		static Integer& precomputeInvProd(Integer& res, const Integer& m1, const Integer& m0)
		{
			inv(res, m0, m1);
			return res *= m0; // res <-- (m0^{-1} mod m1)*m0
		}

		static DomainElement& precomputeInvP0(DomainElement& invP0, const Domain& D1, const Integer& P0)
		{
			return D1.invin( D1.init(invP0, P0) ); // res <-- (P0^{-1} mod m1)
		}


		static Integer& smallbigreconstruct(Integer& u1, const Integer& u0, const Integer& invprod)
		{
			u1 -= u0;	  // u1 <-- (u1-u0)
			u1 *= invprod;    // u1 <-- (u1-u0)( m0^{-1} mod m1 ) m0
			return u1 += u0;  // u1 <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0
		}


		static Integer& normalize(Integer& u1, Integer& tmp, const Integer& m1)
		{
			if (u1 < 0)
				tmp += m1;
			else
				tmp -= m1;
			return ((absCompare(u1,tmp) > 0)? u1 = tmp : u1 );
		}


		static Integer& fieldreconstruct(Integer& res, const Domain& D1, const DomainElement& u1, const Integer& r0, const DomainElement& invP0, const Integer& P0)
		{
                    	// u0 = r0 mod m1
			DomainElement u0; D1.init(u0, r0);
			if (D1.areEqual(u1, u0))
				return res=r0;
			else
				return fieldreconstruct(res, D1, u1, u0, r0, invP0, P0);
		}

		static Integer& fieldreconstruct(Integer& res, const Domain& D1, const DomainElement& u1, DomainElement& u0, const Integer& r0, const DomainElement& invP0, const Integer& P0)
		{
			// u0 and m0 are modified
			D1.negin(u0);   	// u0 <-- -u0
			D1.addin(u0,u1);   	// u0 <-- u1-u0
			D1.mulin(u0, invP0);    // u0 <-- (u1-u0)( m0^{-1} mod m1 )
			D1.convert(res, u0);    // res <-- (u1-u0)( m0^{-1} mod m1 )         and res <  m1
			res *= P0;      	// res <-- (u1-u0)( m0^{-1} mod m1 ) m0      and res <= (m0m1-m0)
			return res += r0;	// res <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0 and res <  m0m1
		}

	};

}


#endif //__LINBOX_cra_full_multip_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
