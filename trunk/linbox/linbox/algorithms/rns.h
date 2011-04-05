/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2011 LinBox
 * Written by BB <brice.boyer@imag.fr>
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

/*! @file algorithms/rns.h
 *  @ingroup algorithms
 * @brief <b>R</b>esidue <b>N</b>umber <b>S</b>ystem tools.
 * NO DOC
 */

#ifndef __LINBOX_algorithms_rns_H
#define __LINBOX_algorithms_rns_H

#include "linbox/integer.h"
#include <givaro/givrns.h> // Chinese Remainder of an array of elements

#include <givaro/givrnsfixed.h>    // Chinese Remainder with fixed primes


namespace LinBox
{

	/*! RNS.
	 * Creates a RNS than can recover any number between \c 0 and \c q-1 if
	 * \c Unsigned=true or \c -q+1 and \c q-1 otherwise (where \c q=2<up>\c
	 * _bits_</up>).
	 * @todo template by field and ring
	 */
	template<bool Unsigned>
	class RNS {
	private:
		typedef std::vector<unsigned long>  Fvect ;
		typedef std::vector<integer>        Ivect ;
		typedef Modular<double>             Field ;

		Fvect               _primes_; //!< vector of integers, pairwise coprime (or pairwise different primes)
		unsigned long         _size_; //!< number of primes
		unsigned long          _bit_; //!< max number of bit reachable. @todo why not max int reachable to save maybe a couple primes ? bof.
		integer             _maxint_; //!< max reachable integer.
		integer             _midint_; //!< in signed case, the mid one.
		unsigned long           _ps_; //!< prime size (minimum)


		typedef RNSsystem<Integer, Field >      CRTSystem;
		typedef typename CRTSystem::domains       Domains;
		typedef typename CRTSystem::array        Elements;
		// typedef typename CRTSystem::ring             Ring;

		CRTSystem     _CRT_ ;
		Domains _PrimeDoms_ ;

#ifdef __LINBOX_HAVE_IML
#endif
	public:
		/*! Create a RNS able to recover any integer of at most l bits.
		 * Builds \c ceil(l/ps) different primes for the RNS
		 * @param l  max recoverable bits
		 * @param ps bitsize of the primes (defaulting to 21 because...)
		 */
		RNS(unsigned long l, unsigned long ps=21) ;
		/*x Create a RNS with given primes.
		 * @param primes given basis of primes
		 * @param l      recoverable bits. If not given or 0, then it is computed. Giving it will reduce initialization time.
		 * @warning both \c l and \c primes 'integrity' are not checked. Use \c checkRNS() member to perform the check.
		 *
		 */
		// RNS(Ivect & primes, unsigned long l = 0) ;

		// void addPrime( double prime);

		// void addPrime( unsigned long newl);

		// bool checkRNS() ;

		// cra
		/*! Inits cra.
		 */
		void initCRA() ;
		/*! Computes \c result corresponding to the \c residues.
		 *
		 */
		void cra(integer & result, const std::vector<double> & residues);
		/*! Computes \c result corresponding to the \c residues.
		 *
		 */
		void cra(std::vector<integer> & result, const std::vector<std::vector<double> > & residues);

		/*! Computes \c result corresponding to the iteration.
		 *
		 */
		template<typename Iteration>
		void cra(Ivect & result, Iteration & iter) ;

		template<class Tinteger, class Tresidue>
		void cra(Tinteger & result, Tresidue & residues);

		template<class Tinteger, class Tresidue>
		void convert(Tinteger & result, Tresidue & residues) ;

		// mixed radix
	};


	template<bool Unsigned>
	class RNSfixed {
	private:
		typedef std::vector<unsigned long>  Fvect ;
		typedef std::vector<integer>        Ivect ;
		typedef Modular<double>             Field ;

		Fvect               _primes_; //!< vector of integers, pairwise coprime (or pairwise different primes)
		unsigned long         _size_; //!< number of primes
		unsigned long          _bit_; //!< max number of bit reachable. @todo why not max int reachable to save maybe a couple primes ? bof.
		integer             _maxint_; //!< max reachable integer.
		integer             _midint_; //!< in signed case, the mid one.
		unsigned long           _ps_; //!< prime size (minimum)


		typedef RNSsystemFixed<Integer>       CRTSystemFixed;
		typedef CRTSystemFixed::array         Prime_t;
		// typedef typename CRTSystem::domains       Domains;
		// typedef typename CRTSystem::array        Elements;
		// typedef typename CRTSystem::ring             Ring;

		CRTSystemFixed     _CRT_ ;
		Prime_t         _Primes_ ;

#ifdef __LINBOX_HAVE_IML
#endif
	public:
		/*! Create a RNSfixed able to recover any integer of at most l bits.
		 * Builds \c ceil(l/ps) different primes for the RNSfixed
		 * @param l  max recoverable bits
		 * @param ps bitsize of the primes (defaulting to 21 because...)
		 */
		RNSfixed(unsigned long l, unsigned long ps=21) ;
		/*x Create a RNSfixed with given primes.
		 * @param primes given basis of primes
		 * @param l      recoverable bits. If not given or 0, then it is computed. Giving it will reduce initialization time.
		 * @warning both \c l and \c primes 'integrity' are not checked. Use \c checkRNSfixed() member to perform the check.
		 *
		 */
		// RNSfixed(Ivect & primes, unsigned long l = 0) ;

		// void addPrime( double prime);

		// void addPrime( unsigned long newl);

		// bool checkRNSfixed() ;

		// cra
		/*! Inits cra.
		 */
		void initCRA() ;
		/*! Computes \c result corresponding to the \c residues.
		 *
		 */
		void cra(integer & result, const std::vector<double> & residues);
		/*! Computes \c result corresponding to the \c residues.
		 *
		 */
		void cra(std::vector<integer> & result, const std::vector<std::vector<double> > & residues);

		/*! Computes \c result corresponding to the iteration.
		 *
		 */
		template<class Iteration>
		void cra(Ivect & result, Iteration & iter) ;

		template<class Tinteger, class Tresidue>
		void cra(Tinteger & result, Tresidue & residues);

		template<class Tinteger, class Tresidue>
		void convert(Tinteger & result, Tresidue & residues) ;

		// mixed radix
	};

}

#include "rns.inl"

#endif // __LINBOX_algorithms_rns_H
