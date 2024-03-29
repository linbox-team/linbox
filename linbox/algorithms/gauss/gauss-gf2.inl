/* linbox/algorithms/gauss-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <03 Dec 21 09:06:58 Jean-Guillaume.Dumas@imag.fr>
 *
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
 *.
 */
#ifndef __LINBOX_gauss_gf2_INL
#define __LINBOX_gauss_gf2_INL
// SparseSeqMatrix is container< container< size_t > >

#include "linbox/algorithms/gauss.h"
#include "linbox/util/commentator.h"
#include <utility>

#ifdef __LINBOX_ALL__ //BB: ???
#ifndef __LINBOX_COUNT__
#define __LINBOX_COUNT__
#endif
#ifndef __LINBOX_OFTEN__
#define __LINBOX_OFTEN__  __LINBOX_ALL__
#endif
#ifndef __LINBOX_FILLIN__
#define __LINBOX_FILLIN__
#endif
#endif

namespace LinBox
{
	// Specialization over GF2
	template <class SparseSeqMatrix, class Perm>
	inline size_t&
	GaussDomain<GF2>::InPlaceLinearPivoting (size_t &Rank,
						 bool          &determinant,
						 SparseSeqMatrix        &LigneA,
						 Perm           &P,
						 size_t Ni,
						 size_t Nj) const
	{
		// Requirements : LigneA is an array of sparse rows
		// In place (LigneA is modified)
		// With reordering (D is a density type. Density is allocated here)
		//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
		commentator().start ("Gaussian elimination with reordering over GF2",
				   "IPLRGF2", Ni);
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Gaussian QLUP elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

#ifdef __LINBOX_COUNT__
		long long nbelem = 0;
#endif

		determinant = true;
		// allocation of the column density
		std::vector<size_t> col_density (Nj);


		// assignment of LigneA with the domain object
		for (size_t jj = 0; jj < Ni; ++jj)
			for (size_t k = 0; k < LigneA[jj].size(); ++k)
				++col_density[LigneA[jj][k]];

		long last = (long)Ni - 1;
		long c;
		Rank = 0;

#ifdef __LINBOX_OFTEN__
		long sstep = last/40;
		if (sstep > __LINBOX_OFTEN__) sstep = __LINBOX_OFTEN__;
#else
		long sstep = 1000;
#endif
		// Elimination steps with reordering

		typename SparseSeqMatrix::iterator LigneA_k = LigneA.begin();
		for (long k = 0; k < last; ++k, ++LigneA_k) {
			long p = k, s = 0;

#ifdef __LINBOX_FILLIN__
			if ( ! (k % 100) )
#else
				if ( ! (k % sstep) )
#endif
				{
					commentator().progress (k);
#ifdef __LINBOX_FILLIN__
					long sl(0);
					for (size_t l = 0; l < Ni; ++l)
						sl += LigneA[(size_t)l].size ();

					commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
					<< "Fillin (" << Rank << "/" << Ni << ") = "
					<< sl
					<< " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
					<< double(sl)/double(Ni-k) << " avg)"
					<< std::endl;
#endif
				}

			long l;
			for(l = k; l < static_cast<long>(Ni); ++l) {
				if ( (s = (long)LigneA[(size_t)l].size()) != 0 ) {
					p = l;
					break;
				}
			}

			if (s) {
				// Row permutation for the sparsest row
				for (; l < static_cast<long>(Ni); ++l) {
				long sl;
					if (((sl = (long)LigneA[(size_t)l].size ()) < s) && (sl)) {
						s = sl;
						p = l;
					}
					}

				if (p != k) {
					//                         std::cerr << "Permuting rows: " << k << " <--> " << p << std::endl;
					std::swap( *LigneA_k, LigneA[(size_t)p]);
				}


				SparseFindPivotBinary (*LigneA_k, Rank, c, col_density, determinant);

				if (c != -1) {
					long ll;
					if ( c != (static_cast<long>(Rank)-1) ) {
						P.permute(Rank-1,(size_t)c);
						for (ll=0      ; ll < k ; ++ll)
							permuteBinary( LigneA[(size_t)ll], Rank, c);
					}
					long npiv=(long)LigneA_k->size();
					for (ll = k+1; ll < static_cast<long>(Ni); ++ll) {
						bool elim=false;
						eliminateBinary (elim, LigneA[(size_t)ll], *LigneA_k, Rank, c, (size_t)npiv, col_density);
					}
				}

				// LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
				nbelem += LigneA_k->size ();
#endif
			}
			// LigneA.write(rep << "U:= ", Tag::FileFormat::Maple) << std::endl;
		}//for k

		SparseFindPivotBinary ( LigneA[(size_t)last], Rank, c, determinant);
		if (c != -1) {
			if ( c != (static_cast<long>(Rank)-1) ) {
				P.permute(Rank-1,(size_t)c);
				for (long ll=0      ; ll < last ; ++ll)
					permuteBinary( LigneA[(size_t)ll], Rank, c);
			}
		}

#ifdef __LINBOX_COUNT__
		nbelem += LigneA[(size_t)last].size ();
		commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
		<< "Left elements : " << nbelem << std::endl;
#endif

#ifdef __LINBOX_FILLIN__
		long sl(0);
		for (size_t l=0; l < Ni; ++l)
			sl += LigneA[(size_t)l].size ();

		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Fillin (" << Rank << "/" << Ni << ") = " << sl
		<< std::endl;
#endif

		if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
			determinant = false;

		integer card;
		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Determinant : " << determinant
		<< " over GF (2)" << std::endl;

		// LigneA.write(rep << "U:= ", Tag::FileFormat::Maple) << ':' << std::endl;

		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Rank : " << Rank
		<< " over GF (2)" << std::endl;
		commentator().stop ("done", 0, "IPLRGF2");



		return Rank;
	}

	// Specialization over GF2
	template <class SparseSeqMatrix, class Perm> inline size_t&
	GaussDomain<GF2>::QLUPin (size_t &Rank,
				  bool          &determinant,
				  Perm          &Q,
				  SparseSeqMatrix        &LigneL,
				  SparseSeqMatrix        &LigneA,
				  Perm          &P,
				  size_t Ni,
				  size_t Nj) const
	{
		linbox_check( Q.coldim() == Q.rowdim() );
		linbox_check( P.coldim() == P.rowdim() );
		linbox_check( Q.coldim() == LigneL.size() );

		// Requirements : LigneA is an array of sparse rows
		// In place (LigneA is modified)
		// With reordering (D is a density type. Density is allocated here)
		//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
		commentator().start ("Gaussian elimination with reordering over GF2",
				   "IPLRGF2", Ni);
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Gaussian QLUP elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

#ifdef __LINBOX_COUNT__
		long long nbelem = 0;
#endif

		determinant = true;
		// allocation of the column density
		std::vector<size_t> col_density (Nj);


		for(typename SparseSeqMatrix::iterator LigneL_it = LigneL.begin() ;
		    LigneL_it != LigneL.end(); ++LigneL_it)
			LigneL_it->reserve(16);

		std::deque<std::pair<size_t,size_t> > invQ;

		// assignment of LigneA with the domain object
		for (size_t jj = 0; jj < Ni; ++jj)
			for (size_t k = 0; k < LigneA[jj].size (); k++)
				++col_density[LigneA[jj][(size_t)k]];

		long last = (long)Ni - 1;
		long c;
		Rank = 0;

#ifdef __LINBOX_OFTEN__
		long sstep = last/40;
		if (sstep > __LINBOX_OFTEN__) sstep = __LINBOX_OFTEN__;
#else
		long sstep = 1000;
#endif
		// Elimination steps with reordering

		typename SparseSeqMatrix::iterator LigneA_k = LigneA.begin();
		for (long k = 0; k < last; ++k, ++LigneA_k) {
			long p = k, s = 0;

#ifdef __LINBOX_FILLIN__
			if ( ! (k % 100) )
#else
				if ( ! (k % sstep) )
#endif
				{
					commentator().progress (k);
#ifdef __LINBOX_FILLIN__
					long sl(0);
					for (size_t l = 0; l < Ni; ++l)
						sl += LigneA[(size_t)l].size ();

					commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
					<< "Fillin (" << Rank << "/" << Ni << ") = "
					<< sl
					<< " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
					<< double(sl)/double(Ni-k) << " avg)"
					<< std::endl;
#endif
				}

			long l;
			for(l = k; l < static_cast<long>(Ni); ++l) {
				if ( (s = (long)LigneA[(size_t)l].size()) ) {
					p = l;
					break;
				}
			}

			if (s) {
				// Row permutation for the sparsest row
				for (; l < static_cast<long>(Ni); ++l){
				long sl;
					if (((sl = (long) LigneA[(size_t)l].size ()) < s) && (sl)) {
						s = sl;
						p = l;
					}
					}

				if (p != k) {
					//                         std::cerr << "Permuting rows: " << k << " <--> " << p << std::endl;
					invQ.emplace_front((size_t)k, (size_t)p);
					std::swap( *LigneA_k, LigneA[(size_t)p]);
					std::swap( LigneL[(size_t)k], LigneL[(size_t)p]);
				}


				SparseFindPivotBinary (*LigneA_k, Rank, c, col_density, determinant);

				if (c != -1) {
					long ll;
					if ( c != (static_cast<long>(Rank)-1) ) {
						P.permute(Rank-1,(size_t)c);
						for (ll=0      ; ll < k ; ++ll)
							permuteBinary( LigneA[(size_t)ll], Rank, c);
					}
					long npiv=(long)LigneA_k->size();
					for (ll = k+1; ll < static_cast<long>(Ni); ++ll) {
						bool elim=false;
						eliminateBinary (elim, LigneA[(size_t)ll], *LigneA_k, Rank, c, (size_t)npiv, col_density);
						if(elim) LigneL[(size_t)ll].emplace_back(Rank-1);
					}
				}

				// LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
				nbelem += LigneA_k->size ();
#endif
			}
			LigneL[(size_t)k].emplace_back((size_t)k);
			//  LigneL.write(rep << "L:= ", Tag::FileFormat::Maple) << std::endl;
			//  LigneA.write(rep << "U:= ", Tag::FileFormat::Maple) << std::endl;
		}//for k

		SparseFindPivotBinary ( LigneA[(size_t)last], Rank, c, determinant);
		if (c != -1) {
			if ( c != (static_cast<long>(Rank)-1) ) {
				P.permute(Rank-1,(size_t)c);
				for (long ll=0      ; ll < last ; ++ll)
					permuteBinary( LigneA[(size_t)ll], Rank, c);
			}
		}

		LigneL[(size_t)last].emplace_back((size_t)last);

#ifdef __LINBOX_COUNT__
		nbelem += LigneA[(size_t)last].size ();
		commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
		<< "Left elements : " << nbelem << std::endl;
#endif

#ifdef __LINBOX_FILLIN__
		long sl(0);
		for (size_t l=0; l < Ni; ++l)
			sl += LigneA[(size_t)l].size ();

		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Fillin (" << Rank << "/" << Ni << ") = " << sl
		<< std::endl;
#endif

		if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
			determinant = false;

		integer card;
		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Determinant : " << determinant
		<< " over GF (2)" << std::endl;

		for(std::deque<std::pair<size_t,size_t> >::const_iterator it = invQ.begin(); it!=invQ.end();++it)
			Q.permute( it->first, it->second );

#if 0
		Q.write(rep << "Q:= ", Tag::FileFormat::Maple) << ':' << std::endl;
		LigneL.write(rep << "L:= ", Tag::FileFormat::Maple) << ':' << std::endl;
		LigneA.write(rep << "U:= ", Tag::FileFormat::Maple) << ':' << std::endl;
		P.write(rep << "P:= ", Tag::FileFormat::Maple) << ':' << std::endl;
#endif
		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Rank : " << Rank
		<< " over GF (2)" << std::endl;
		commentator().stop ("done", 0, "IPLRGF2");



		return Rank;
	}


} // namespace LinBox

#endif // __LINBOX_gauss_gf2_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
