/* linbox/algorithms/cra-builder-full-multip.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Time-stamp: <11 Feb 25 13:39:19 Jean-Guillaume.Dumas@imag.fr>
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

/*!@file algorithms/cra-builder-full-multip.h
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

	/** @brief Chinese remaindering of a vector of elements without early termination.
	 * @ingroup CRA
	 * @bib
	 * - Jean-Guillaume Dumas, Thierry Gautier et Jean-Louis Roch.  <i>Generic design
	 * of Chinese remaindering schemes</i>  PASCO 2010, pp 26-34, 21-23 juillet,
	 * Grenoble, France.
     *
     * The idea is that each "shelf" contains a vector of residues with some modulus.
     * We want to combine shelves that have roughly the same size of modulus, for
     * efficiency. The method is that any submitted residue is assigned to a unique
     * shelf according to log2(log(modulus)), as computed by the getShelf() helper.
     * When two residues belong on the same shelf, they are combined and re-assigned
     * to another shelf, recursively.
	 */
	template<class Domain_Type>
	struct CRABuilderFullMultip {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef CRABuilderFullMultip<Domain>		Self_t;

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
		const double				LOGARITHMIC_UPPER_BOUND; // log2 of upper bound
		double totalsize_ = 0.; // log2 of the current modulus
        size_t dimension_ = 0; // dimension of the vector being reconstructed
        bool collapsed_ = false;
        bool normalized_ = false;
        // INVARIANT: shelves_.empty() || shelves_.back().occupied
        // INVARIANT: forall (shelf : shelves_) { shelf.residue.size() == dimension_ }

	public:
        friend std::ostream& operator<< (std::ostream& out, const Self_t& cra) {
            std::ostringstream report;
            report << "CRA Builder: "
                   << "[BoundedTermination] [MultipleReconstructions]";
            return out << report.str();
        }

        /** @brief Creates a new vector CRA object.
         * @param bnd  upper bound on the natural logarithm of the result
         * @param dim  dimension of the vector to be reconstructed
         */
		CRABuilderFullMultip(const double bnd=0.0, size_t dim=0) :
			LOGARITHMIC_UPPER_BOUND(bnd), dimension_(dim)
		{
#if __LB_CRA_REPORTING__
            std::clog << *this << std::endl;
#endif
        }
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
            const integer& Dval = mod_to_integer(D);
            double logD = Givaro::naturallog(Dval);
            auto cur = getShelf(logD);

            totalsize_ += Givaro::logtwo(Dval);

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

            // shelf is nonempty, so we incorporate the new result there
            {
                auto invprod = precompInv(shelves_[cur].mod(), D);
                auto r_it = shelves_[cur].residue.begin();
                for (size_t i=0; i < e_size; ++i, ++e_it, ++r_it) {
                    reconstruct(*r_it, shelves_[cur].mod(), *e_it, invprod, D);
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
                    combineShelves(shelves_[next], shelves_[cur]);
                    shelves_[cur].occupied = false;
                } else {
                    // put cur shelf data in next shelf position
                    std::swap(shelves_[cur], shelves_[next]);
                }

                cur = next;
            }
		}

		//! result
		inline const std::vector<Integer>& result (bool normalized=true) const
		{
            normalize();
            return shelves_.back().residue;
        }

        template <class Vect>
        inline Vect& result(Vect& r, bool normalized=true) const
        {
            r.resize(dimension_);
            result_iter(r.begin());
            return r;
        }

        template <class Iter>
        void result_iter (Iter r_it, bool normalized=true) const {
            if (shelves_.empty()) {
                for (size_t i=0; i < dimension_; ++i)
                    *r_it = 0;
            }
            else {
                normalize();
                std::copy_n(shelves_.back().residue.begin(), dimension_, r_it);
            }
        }

        // alias for result
		inline const std::vector<Integer>& getResidue() const
		{
			return result();
		}

        // alias for result
		template<class Vect>
		inline Vect& getResidue(Vect& r) const
		{
			return result(r);
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

        size_t getDimension() const
        { return dimension_; }

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

        /** @brief Returns a reference to D.
         * This is needed to automatically handle whether D is a Domain or an actual
         * integer.
         */
        static inline const integer& mod_to_integer(const Integer& D) {
            return D;
        }

        /** @brief Returns the characteristic of D.
         */
        template <class Domain>
        static inline integer mod_to_integer(const Domain& D) {
            integer m;
            D.characteristic(m);
            return m;
        }

        /** @brief Incorporates the residue in src into dest and updates the modulus.
         */
        static inline void combineShelves(Shelf& dest, const Shelf& src) {
            // assumption: dest is already occupied
            auto invprod = precompInv(dest.mod(), src.mod());
			auto src_it = src.residue.begin();
            for (auto dest_it = dest.residue.begin();
                 dest_it != dest.residue.end();
                 ++dest_it, ++src_it)
            {
                reconstruct(*dest_it, dest.mod(), *src_it, invprod, src.mod());
            }
            dest.mod.mulin(src.mod());
            dest.logmod += src.logmod;
            dest.count += src.count;
        }

        /** @brief Expands the shelves as necessary so that the given index
         * exists in the array.
         */
        static inline void ensureShelf(size_t index, std::vector<Shelf>& shelves, size_t dim) {
            while (index >= shelves.size()) {
                shelves.emplace_back(dim);
            }
        }

        /** @brief Collapses all shelves by combining residues.
         *
         * After this, there will be a single (top) shelf containing the current
         * full residue.
         */
        void collapse() const {
            if (collapsed_) return;
            auto& ncshelves = const_cast<std::vector<Shelf>&>(shelves_);
            if (ncshelves.empty()) {
                ncshelves.emplace_back(dimension_);
                ncshelves.back().occupied = true;
            }
            else {
                auto cur = ncshelves.begin();
                while (! cur->occupied) ++cur;
                auto next = cur;
                while(true) {
                    if (++next == ncshelves.end()) break;
                    while (! next->occupied) ++next;
                    combineShelves(*next, *cur);
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

        /** @brief Collapses (if necessary) the top shelf and normalizes the
         * result into the symmetric modulus range.
         */
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

        // return: (m0^-1 mod m1) * m0
        static Integer precompInv(const Integer& m1, const Integer& m0)
		{
            Integer res;
			inv(res, m0, m1); // res <- m0^{-1} mod m1
			res *= m0; // res <-- (m0^{-1} mod m1)*m0
            return res;
		}

        // return: m1^-1 mod m0
        template <class Domain>
        static typename Domain::Element precompInv(const Integer& m1, const Domain& D0)
        {
            typename Domain::Element res;
            D0.invin( D0.init(res, m1) ); // res <- m1^{-1} mod m0
            return res;
        }

        // precond: invprod = (m0^-1 mod m1)*m0
        // return: u1' such that u1' mod m0 = u0 and u1' mod m1 = u1
        template <typename IntegerLike>
        static Integer& reconstruct(Integer& u1, const Integer&, const IntegerLike& u0, const Integer& invprod, const Integer&)
        {
			u1 -= u0;	  // u1 <-- (u1-u0)
			u1 *= invprod;    // u1 <-- (u1-u0)( m0^{-1} mod m1 ) m0
			return u1 += u0;  // u1 <-- u0 + (u1-u0)( m0^{-1} mod m1 ) m0
        }

        // precond: invprod = m1^-1 mod m0
        // return: u1' such that u1' mod m0 = u0 and u1' mod m1 = u1
        template <class Domain>
        static Integer& reconstruct(Integer& u1, const Integer& m1, const typename Domain::Element& u0, const typename Domain::Element& invprod, const Domain& D0)
        {
            typename Domain::Element u1elt;
            D0.init(u1elt, u1);
            if (D0.areEqual(u0, u1elt))
                return u1;

            D0.negin(u1elt); // u1elt <-- -u1
            D0.addin(u1elt, u0); // u1elt <-- u0 - u1
            D0.mulin(u1elt, invprod); // u1elt <-- (u0-u1)/m1 mod m0
            integer temp;
            D0.convert(temp, u1elt); // temp <-- (u0-u1)/m1 mod m0
            temp *= m1; // temp <-- ((u0-u1)/m1 mod m0) * m1
            u1 += temp; // u1 <-- u1 + ((u0-u1)/m1 mod m0) * m1
            return u1;
        }

#ifdef __LB_CRA_TIMING__
    public:
        std::ostream& reportTimes(std::ostream& os) const
        {
            return os <<  "FullMultip CRA total size:" << totalsize_;
        }
#endif

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
