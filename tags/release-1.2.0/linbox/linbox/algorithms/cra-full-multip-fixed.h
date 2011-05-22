/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr>
 *
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

/*! @file algorithms/cra-full-multip-fixed.h
 * @ingroup algorithms
 * @brief CRA for multi-residues.
 *
 * An upper bound is given on the size of the data to reconstruct.
 */

#ifndef __LINBOX_cra_full_multip_fixed_H
#define __LINBOX_cra_full_multip_fixed_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

#include "linbox/algorithms/lazy-product.h"
#include "linbox/algorithms/cra-full-multip.h"


namespace LinBox
{



	/*! @ingroup CRA
	 * @brief Chinese Remaindering Algorithm for multiple residues.
	 * An upper bound is given on the size of the data to reconstruct.
	 */
	template<class Domain_Type>
	struct FullMultipFixedCRA : FullMultipCRA<Domain_Type> {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element 	DomainElement;
		typedef FullMultipFixedCRA<Domain> 	Self_t;

		typedef std::vector<double>::iterator        DoubleVect_Iterator ;
		typedef std::vector< bool >::iterator          BoolVect_Iterator ;
		typedef std::vector< LazyProduct >::iterator   LazyVect_Iterator ;
		typedef std::vector< Integer >                           IntVect ;
		typedef IntVect::iterator                       IntVect_Iterator ;
		typedef std::vector< IntVect >::iterator    IntVectVect_Iterator ;
		typedef IntVect::const_iterator            IntVect_ConstIterator ;

	protected:
		const size_t				size;

	private :
		/*! \internal
		 *  Intialize the Radix ladder.
		 */
		void _initialize ()
		{
			this->RadixSizes_.resize(1);
			this->RadixPrimeProd_.resize(1);
			this->RadixResidues_.resize(1);
			this->RadixOccupancy_.resize(1);
			this->RadixOccupancy_.front() = false;
		}

	public:
		/*! Constructor.
		 * @param p is a pair such that
		 * - \c p.first is the size of a residue (ie. it would be 1 for \"FullSingle\")
		 * - \c p.second is the theoretical upperbound (natural logarithm) on the size of the integers to reconstruct.
		 * .
		 */
		FullMultipFixedCRA(const std::pair<size_t,double>& p ) :
		       	FullMultipCRA<Domain>(p.second), size(p.first)
	       	{
		       	this->_initialize();
		}

		/*! Intialize to the first residue/prime.
		 * @param D domain
		 * @param e iterator on the first residue
		 * @pre any CRA should first call \c initialize before \c progress
		 */
		template<class Iterator>
		void initialize (const Domain& D, Iterator& e)
		{
			this->_initialize();
			this->progress(D, e);
		}

		/*! Add a new residue (ie take into account a new prime).
		 * @param D domain
		 * @param e iterator for the new residue, for instance, a
		 * <code>std::vector<T>::iterator</code>.
		 */
		template<class Iterator>
		void progress (const Domain& D, Iterator& e)
		{
			// Radix shelves
			DoubleVect_Iterator  _dsz_it = this->RadixSizes_.begin();
			LazyVect_Iterator    _mod_it = this->RadixPrimeProd_.begin();
			IntVectVect_Iterator _tab_it = this->RadixResidues_.begin();
			BoolVect_Iterator    _occ_it = this->RadixOccupancy_.begin();

			IntVect ri(this->size);
			LazyProduct mi;
			double di;
			if (*_occ_it) {
				// If lower shelf is occupied
				// Combine it with the new residue
				// The for loop will transmit the resulting
				// combination to the upper shelf
				Iterator e_it = e;
				IntVect_Iterator       ri_it = ri.begin();
				IntVect_ConstIterator  t0_it = _tab_it->begin();

				DomainElement invP0;
				this->precomputeInvP0(invP0, D, _mod_it->operator()() );

				for( ; ri_it != ri.end(); ++e_it, ++ri_it, ++ t0_it){
					this->fieldreconstruct(*ri_it, D, *e_it, *t0_it, invP0,
							       (*_mod_it).operator()() );
				}
				Integer tmp;
				D.characteristic(tmp);
				double ltp = ::Givaro::naturallog(tmp);
				di = *_dsz_it + ltp;
				this->totalsize += ltp;
				mi.mulin(tmp);
				mi.mulin(*_mod_it);
				*_occ_it = false;
			}
			else {
				// Lower shelf is free
				// Put the new residue here and exit
				Integer tmp;
				D.characteristic(tmp);
				double ltp = ::Givaro::naturallog(tmp);
				_mod_it->initialize(tmp);
				*_dsz_it = ltp;
				this->totalsize += ltp;
				Iterator e_it = e;
				_tab_it->resize(this->size);
				IntVect_Iterator t0_it= _tab_it->begin();
				for( ; t0_it != _tab_it->end(); ++e_it, ++ t0_it){
					D.convert(*t0_it, *e_it);
				}
				*_occ_it = true;
				return;
			}

			// We have a combination to put in the upper shelf
			for(++_dsz_it, ++_mod_it, ++_tab_it, ++_occ_it ;
			    _occ_it != this->RadixOccupancy_.end() ;
			    ++_dsz_it, ++_mod_it, ++_tab_it, ++_occ_it) {
				if (*_occ_it) {
					// This shelf is occupied
					// Combine it with the new combination
					// The loop will try to put it on the upper shelf
					IntVect_Iterator      ri_it = ri.begin();
					IntVect_ConstIterator t_it  = _tab_it->begin();

					Integer invprod;
					this->precomputeInvProd(invprod, mi(), _mod_it->operator()());

					for( ; ri_it != ri.end(); ++ri_it, ++ t_it)
						this->smallbigreconstruct(*ri_it, *t_it, invprod);

					// Product (lazy) computation
					mi.mulin(*_mod_it);

					// Moding out
					for(ri_it = ri.begin() ; ri_it != ri.end(); ++ri_it) {
						*ri_it %= mi();
					}

					di += *_dsz_it;
					*_occ_it = false;
				}
				else {
					// This shelf is free
					// Put the new combination here and exit
					*_dsz_it = di;
					*_mod_it = mi;
					*_tab_it = ri;
					*_occ_it = true;
					return;
				}
			}
			// All the shelfves were occupied
			// We create a new top shelf
			// And put the new combination there
			this->RadixSizes_.push_back( di );
			this->RadixResidues_.push_back( ri );
			this->RadixPrimeProd_.push_back( mi );
			this->RadixOccupancy_.push_back ( true );
		}

		/*! Compute the result.
		 * moves low occupied shelves up.
		 * @param[out] d an iterator for the result.
		 */
		template<class Iterator>
		Iterator& result (Iterator &d)
		{
			LazyVect_Iterator     _mod_it = this->RadixPrimeProd_.begin();
			IntVectVect_Iterator  _tab_it = this->RadixResidues_.begin();
			BoolVect_Iterator     _occ_it = this->RadixOccupancy_.begin();
			LazyProduct Product;
			// We have to find to lowest occupied shelf
			for( ; _occ_it != this->RadixOccupancy_.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
				if (*_occ_it) {
					// Found the lowest occupied shelf
					Product = *_mod_it;
					Iterator t0_it = d;
					IntVect_Iterator t_it = _tab_it->begin();
					if (++_occ_it == this->RadixOccupancy_.end()) {
						// It is the only shelf of the radix
						// We normalize the result and output it
						for( ; t_it != _tab_it->end(); ++t0_it, ++t_it)
							normalize(*t0_it = *t_it, *t_it,
								  _mod_it->operator()());
						this->RadixPrimeProd_.resize(1);
						return d;
					}
					else {
						// There are other shelves
						// The result is initialized with this shelf
						// The for loop will combine the other shelves m with the actual one
						for( ; t_it != _tab_it->end(); ++t0_it, ++t_it)
							*t0_it  = *t_it;
						++_mod_it; ++_tab_it;
						break;
					}
				}
			}
			for( ; _occ_it != this->RadixOccupancy_.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
				if (*_occ_it) {
					// This shelf is occupied
					// We need to combine it with the actual value of the result
					Iterator t0_it = d;
					IntVect_ConstIterator t_it = _tab_it->begin();

					Integer invprod;
					this->precomputeInvProd(invprod, Product(), _mod_it->operator()() );

					for( ; t_it != _tab_it->end(); ++t0_it, ++t_it)
						this->smallbigreconstruct(*t0_it, *t_it, invprod);

					// Overall product computation
					Product.mulin(*_mod_it);

					// Moding out and normalization
					t0_it = d;
					for(size_t i=0;i<this->size; ++i, ++t0_it) {
						*t0_it %= Product();
						Integer tmp(*t0_it);
						this->normalize(*t0_it, tmp, Product());
					}

				}
			}

			// We put it also the final prime product in the first shelf of products
			// JGD : should we also put the result
			//       in the first shelf of residues and resize it to 1
			//       and set to true the first occupancy and resize it to 1
			//       in case result is not the last call (more progress to go) ?
			this->RadixPrimeProd_.resize(1);
			this->RadixPrimeProd_.front() = Product;
			return d;
		}

	};

}

namespace LinBox
{
	/*! NO DOC..
	 * @ingroup CRA
	 * Version of LinBox::FullMultipCRA for matrices.
	 */
	template<class Domain_Type>
	struct FullMultipBlasMatCRA : FullMultipCRA<Domain_Type> {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element 	DomainElement;
		typedef FullMultipBlasMatCRA<Domain> 	Self_t;

		typedef std::vector<double>::iterator        DoubleVect_Iterator ;
		typedef std::vector< bool >::iterator          BoolVect_Iterator ;
		typedef std::vector< LazyProduct >::iterator   LazyVect_Iterator ;
		typedef std::vector< Integer >                           IntVect ;
		typedef IntVect::iterator                       IntVect_Iterator ;
		typedef std::vector< IntVect >::iterator    IntVectVect_Iterator ;
		typedef IntVect::const_iterator            IntVect_ConstIterator ;

	protected:
		const size_t				size;

	private :
		/*! \internal
		 *  Intialize the Radix ladder.
		 */
		void _initialize ()
		{
			this->RadixSizes_.resize(1);
			this->RadixPrimeProd_.resize(1);
			this->RadixResidues_.resize(1);
			this->RadixOccupancy_.resize(1);
			this->RadixOccupancy_.front() = false;
		}

	public:
		/*! Constructor.
		 * @param p is a pair such that
		 * - \c p.first is the size of a residue, it would be 1 for \"FullSingle\"
		 * - \c p.second is the theoretical upperbound (natural logarithm) on the size of the integers to reconstruct.
		 * .
		 */
		FullMultipBlasMatCRA(const std::pair<size_t,double>& p ) :
		       	FullMultipCRA<Domain>(p.second), size(p.first)
		{ }

		/*! Intialize to the first residue/prime.
		 * @param D domain
		 * @param e
		 * @pre any CRA should first call \c initialize before \c progress
		 */
		template<class Matrix>
		void initialize (const Domain& D, Matrix& e)
		{
			this->_initialize();
			this->progress(D, e);
		}

		/*! Add a new residue (ie take into account a new prime).
		 * @param D domain
		 * @param e
		 */
		template<class Matrix>
		void progress (const Domain& D, Matrix& e)
		{
			// Radix shelves
			DoubleVect_Iterator  _dsz_it = this->RadixSizes_.begin();
			LazyVect_Iterator    _mod_it = this->RadixPrimeProd_.begin();
			IntVectVect_Iterator _tab_it = this->RadixResidues_.begin();
			BoolVect_Iterator    _occ_it = this->RadixOccupancy_.begin();

			IntVect ri(this->size);
			LazyProduct mi;
			double di;
			if (*_occ_it) {
				// If lower shelf is occupied
				// Combine it with the new residue
				// The for loop will transmit the resulting
				// combination to the upper shelf
				typename Matrix::RawIterator e_it = e.rawBegin();
				IntVect_Iterator            ri_it = ri.begin();
				IntVect_ConstIterator       t0_it = _tab_it->begin();

				DomainElement invP0;
			       	this->precomputeInvP0(invP0, D, _mod_it->operator()() );

				for( ; ri_it != ri.end(); ++e_it, ++ri_it, ++ t0_it){
					this->fieldreconstruct(*ri_it, D, *e_it, *t0_it, invP0,
							       (*_mod_it).operator()() );
				}
				Integer tmp;
				D.characteristic(tmp);
				double ltp = ::Givaro::naturallog(tmp);
				di = *_dsz_it + ltp;
				this->totalsize += ltp;
				mi.mulin(tmp);
				mi.mulin(*_mod_it);
				*_occ_it = false;
			}
			else {
				// Lower shelf is free
				// Put the new residue here and exit
				Integer tmp;
				D.characteristic(tmp);
				double ltp = ::Givaro::naturallog(tmp);
				_mod_it->initialize(tmp);
				*_dsz_it = ltp;
				this->totalsize += ltp;
				typename Matrix::RawIterator e_it = e.rawBegin();
				_tab_it->resize(this->size);
				IntVect_Iterator t0_it= _tab_it->begin();
				for( ; t0_it != _tab_it->end(); ++e_it, ++ t0_it){
					D.convert(*t0_it, *e_it);
				}
				*_occ_it = true;
				return;
			}

			// We have a combination to put in the upper shelf
			for(++_dsz_it, ++_mod_it, ++_tab_it, ++_occ_it ;
			    _occ_it != this->RadixOccupancy_.end() ;
			    ++_dsz_it, ++_mod_it, ++_tab_it, ++_occ_it) {
				if (*_occ_it) {
					// This shelf is occupied
					// Combine it with the new combination
					// The loop will try to put it on the upper shelf
					IntVect_Iterator      ri_it = ri.begin();
					IntVect_ConstIterator t_it  = _tab_it->begin();

					Integer invprod;
					this->precomputeInvProd(invprod, mi(), _mod_it->operator()());

					for( ; ri_it != ri.end(); ++ri_it, ++ t_it)
						this->smallbigreconstruct(*ri_it, *t_it, invprod);

					// Product (lazy) computation
					mi.mulin(*_mod_it);

					// Moding out
					for(ri_it = ri.begin() ; ri_it != ri.end(); ++ri_it) {
						*ri_it %= mi();
					}

					di += *_dsz_it;
					*_occ_it = false;
				}
				else {
					// This shelf is free
					// Put the new combination here and exit
					*_dsz_it = di;
					*_mod_it = mi;
					*_tab_it = ri;
					*_occ_it = true;
					return;
				}
			}
			// All the shelfves were occupied
			// We create a new top shelf
			// And put the new combination there
			this->RadixSizes_.push_back( di );
			this->RadixResidues_.push_back( ri );
			this->RadixPrimeProd_.push_back( mi );
			this->RadixOccupancy_.push_back ( true );
		}

		/*! Compute the result.
		 * moves low occupied shelves up.
		 * @param[out] d
		 */
		template<class Matrix>
		Matrix& result (Matrix &d)
		{
			LazyVect_Iterator     _mod_it = this->RadixPrimeProd_.begin();
			IntVectVect_Iterator  _tab_it = this->RadixResidues_.begin();
			BoolVect_Iterator     _occ_it = this->RadixOccupancy_.begin();
			LazyProduct Product;
			size_t _j=0;
			// We have to find to lowest occupied shelf
			for( ; _occ_it != this->RadixOccupancy_.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
				++_j;
				if (*_occ_it) {
					// Found the lowest occupied shelf
					Product = *_mod_it;
					typename Matrix::RawIterator t0_it = d.rawBegin();
					IntVect_Iterator t_it = _tab_it->begin();
					if (++_occ_it == this->RadixOccupancy_.end()) {
						// It is the only shelf of the radix
						// We normalize the result and output it
						for( ; t_it != _tab_it->end(); ++t0_it, ++t_it)
							normalize(*t0_it = *t_it, *t_it,
								  _mod_it->operator()());
						this->RadixPrimeProd_.resize(1);
						return d;
					}
					else {
						// There are other shelves
						// The result is initialized with this shelf
						// The for loop will combine the other shelves m with the actual one
						for( ; t_it != _tab_it->end(); ++t0_it, ++t_it)
							*t0_it  = *t_it;
						++_mod_it; ++_tab_it;
						break;
					}
				}
			}
			for( ; _occ_it != this->RadixOccupancy_.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
				++_j;
				if (*_occ_it) {
					// This shelf is occupied
					// We need to combine it with the actual value of the result
					typename Matrix::RawIterator t0_it = d.rawBegin();
					IntVect_ConstIterator t_it = _tab_it->begin();

					Integer invprod;
					this->precomputeInvProd(invprod, Product(), _mod_it->operator()() );

					for( ; t_it != _tab_it->end(); ++t0_it, ++t_it)
						this->smallbigreconstruct(*t0_it, *t_it, invprod);

					// Overall product computation
					Product.mulin(*_mod_it);

					// Moding out and normalization
					t0_it = d.rawBegin();
					for(size_t i=0;i<this->size; ++i, ++t0_it) {
						*t0_it %= Product();
						Integer tmp(*t0_it);
						this->normalize(*t0_it, tmp, Product());
					}

				}
			}

			// We put it also the final prime product in the first shelf of products
			// JGD : should we also put the result
			//       in the first shelf of residues and resize it to 1
			//       and set to true the first occupancy and resize it to 1
			//       in case result is not the last call (more progress to go) ?
			this->RadixPrimeProd_.resize(1);
			this->RadixPrimeProd_.front() = Product;
			return d;
		}

	};

}

#endif //__LINBOX_cra_full_multip_fixed_H

