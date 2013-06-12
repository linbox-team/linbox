/* linbox/algorithms/gauss.inl
 * Copyright (C) 1999 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <18 Jun 10 15:48:38 Jean-Guillaume.Dumas@imag.fr>
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
// ========================================================================= //
// (C) The Linbox Group 1999
// Calcul de rang par la méthode de Gauss pivot par ligne, sur matrice creuse
// ========================================================================= //

#ifndef __LINBOX_gauss_INL
#define __LINBOX_gauss_INL

#include "linbox/algorithms/gauss.h"
#include "linbox/util/commentator.h"
#include <utility>

#ifdef __LINBOX_ALL__
#define __LINBOX_COUNT__
#define __LINBOX_OFTEN__ __LINBOX_ALL__ // BB: ???
#define __LINBOX_FILLIN__
#endif

namespace LinBox
{
	template <class _Field>
	template <class Matrix, class Perm> inline unsigned long&
	GaussDomain<_Field>::QLUPin (unsigned long &Rank,
				     Element       &determinant,
				     Perm          &Q,
				     Matrix        &LigneL,
				     Matrix        &LigneA,
				     Perm          &P,
				     unsigned long Ni,
				     unsigned long Nj) const
	{
		linbox_check( Q.coldim() == Q.rowdim() );
		linbox_check( P.coldim() == P.rowdim() );
		linbox_check( LigneL.coldim() == LigneL.rowdim() );
		linbox_check( Q.coldim() == LigneL.rowdim() );
		linbox_check( LigneL.coldim() == LigneA.rowdim() );
		linbox_check( LigneA.coldim() == P.rowdim() );

		typedef typename Matrix::Row        Vector;
		typedef typename Vector::value_type E;

		// Requirements : LigneA is an array of sparse rows
		// In place (LigneA is modified)
		// With reordering (D is a density type. Density is allocated here)
		//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
		commentator().start ("Gaussian elimination with reordering",
				   "IPLR", Ni);
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Gaussian QLUP elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

#ifdef __LINBOX_COUNT__
		long long nbelem = 0;
#endif

		Element Eone; field().init(Eone,1UL);
		field().init(determinant,1UL);
		// allocation of the column density
		std::vector<size_t> col_density (Nj);


		for(typename Matrix::RowIterator LigneL_it = LigneL.rowBegin() ;
		    LigneL_it != LigneL.rowEnd(); ++LigneL_it)
			LigneL_it->reserve(16);

		std::deque<std::pair<size_t,size_t> > invQ;

		// assignment of LigneA with the domain object
		for (unsigned long jj = 0; jj < Ni; ++jj)
			for (unsigned long k = 0; k < LigneA[(size_t)jj].size (); k++)
				++col_density[LigneA[(size_t)jj][k].first];

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

		typename Matrix::RowIterator LigneA_k = LigneA.rowBegin(), LigneA_p;
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
					sl +=(long) LigneA[(size_t)l].size ();

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
				if ( (s = (long) LigneA[(size_t)l].size()) ) {
					p = l;
					break;
				}
			}

			if (s) {
				long sl;
				// Row permutation for the sparsest row
				for (; l < static_cast<long>(Ni); ++l)
					if (((sl =(long) LigneA[(size_t)l].size ()) < s) && (sl)) {
						s = sl;
						p = l;
					}

				if (p != k) {
					// std::cerr << "Permuting rows: " << k << " <--> " << p << std::endl;
					invQ.push_front( std::pair<size_t,size_t>((size_t)k,(size_t)p) );
					field().negin(determinant);
					std::swap( *LigneA_k, LigneA[(size_t)p]);
					std::swap( LigneL[(size_t)k], LigneL[(size_t)p]);
				}


				SparseFindPivot (*LigneA_k, Rank, c, col_density, determinant);

				if (c != -1) {
					long ll;
					if ( c != (static_cast<long>(Rank)-1) ) {
						P.permute(Rank-1,(size_t)c);
						for (ll=0      ; ll < k ; ++ll)
							permute( LigneA[(size_t)ll], Rank, c);
					}
					long npiv=(long)LigneA_k->size();
					for (ll = k+1; ll < static_cast<long>(Ni); ++ll) {
						E hc;
						hc.first=(unsigned)Rank-1;
						eliminate (hc.second, LigneA[(size_t)ll], *LigneA_k, Rank, c, (size_t)npiv, col_density);
						if(! field().isZero(hc.second)) LigneL[(size_t)ll].push_back(hc);
					}
				}

				//                     LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
				nbelem += LigneA_k->size ();
#endif
			}
			E one((unsigned)k,Eone);
			LigneL[(size_t)k].push_back(one);
			//                 LigneL.write(rep << "L:= ", FORMAT_MAPLE) << std::endl;
			//                 LigneA.write(rep << "U:= ", FORMAT_MAPLE) << std::endl;
		}//for k

		SparseFindPivot ( LigneA[(size_t)last], Rank, c, determinant);
		if (c != -1) {
			if ( c != (static_cast<long>(Rank)-1) ) {
				P.permute(Rank-1,(size_t)c);
				for (long ll=0      ; ll < last ; ++ll)
					permute( LigneA[(size_t)ll], Rank, c);
			}
		}

		E one((unsigned)last,Eone);
		LigneL[(size_t)last].push_back(one);

#ifdef __LINBOX_COUNT__
		nbelem += LigneA[(size_t)last].size ();
		commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
		<< "Left elements : " << nbelem << std::endl;
#endif

#ifdef __LINBOX_FILLIN__
		long sl(0);
		for (size_t l=0; l < Ni; ++l)
			sl += (long) LigneA[(size_t)l].size ();

		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Fillin (" << Rank << "/" << Ni << ") = " << sl
		<< std::endl;
#endif

		if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
			field().init(determinant,0UL);

		integer card;
		field().write(commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
			 << "Determinant : ", determinant)
		<< " over GF (" << field().cardinality (card) << ")" << std::endl;

		for(std::deque<std::pair<size_t,size_t> >::const_iterator it = invQ.begin(); it!=invQ.end();++it)
			Q.permute( it->first, it->second );


		//             std::ostream& rep = commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT);
		//             Q.write(rep << "Q:= ", FORMAT_MAPLE) << ':' << std::endl;
		//             LigneL.write(rep << "L:= ", FORMAT_MAPLE) << ':' << std::endl;
		//             LigneA.write(rep << "U:= ", FORMAT_MAPLE) << ':' << std::endl;
		//             P.write(rep << "P:= ", FORMAT_MAPLE) << ':' << std::endl;

		commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Rank : " << Rank
		<< " over GF (" << card << ")" << std::endl;
		commentator().stop ("done", 0, "IPLR");



		return Rank;
	}

	template <class _Field>
	template <class Matrix> inline unsigned long&
	GaussDomain<_Field>::InPlaceLinearPivoting (unsigned long &Rank,
						    Element        &determinant,
						    Matrix         &LigneA,
						    unsigned long   Ni,
						    unsigned long   Nj) const
	{
		typedef typename Matrix::Row        Vector;
		typedef typename Vector::value_type E;

		// Requirements : LigneA is an array of sparse rows
		// In place (LigneA is modified)
		// With reordering (D is a density type. Density is allocated here)
		//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
		commentator().start ("Gaussian elimination with reordering",
				   "IPLR", Ni);
		field().write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			  << "Gaussian elimination on " << Ni << " x " << Nj << " matrix, over: ") << std::endl;

#ifdef __LINBOX_COUNT__
		long long nbelem = 0;
#endif

		field().init(determinant,1UL);
		Vector Vzer (0);
		// allocation of the column density
		std::vector<size_t> col_density (Nj);

		// assignment of LigneA with the domain object
		for (unsigned long jj = 0; jj < Ni; ++jj)
			for (unsigned long k = 0; k < LigneA[(size_t)jj].size (); k++)
				++col_density[LigneA[(size_t)jj][k].first];

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
		for (long k = 0; k < last; ++k) {
			unsigned long l;
			long p = k, s = (long)LigneA[(size_t)k].size (), sl;

#ifdef __LINBOX_FILLIN__
			if ( ! (k % 100) ) {
#else
				if ( ! (k % sstep) ) {
#endif
					commentator().progress (k);
#ifdef __LINBOX_FILLIN__
					for (sl = 0, l = 0; l < Ni; ++l)
						sl += (long)LigneA[(size_t)l].size ();

					commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
					<< "Fillin (" << Rank << "/" << Ni << ") = "
					<< sl
					<< " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
					<< double(sl)/double(Ni-k) << " avg)"
					<< std::endl;
#endif
				}

				if (s) {
					// Row permutation for the sparsest row
					for (l = (unsigned long)k + 1; l < (unsigned long)Ni; ++l)
						if (((sl = (long)LigneA[(size_t)l].size ()) < s) && (sl)) {
							s = sl;
							p = (long)l;
						}

					if (p != k) {
						field().negin(determinant);
						Vector vtm = LigneA[(size_t)k];
						LigneA[(size_t)k] = LigneA[(size_t)p];
						LigneA[(size_t)p] = vtm;
					}

					//                     LigneA.write(std::cerr << "BEF, k:" << k << ", Rank:" << Rank << ", c:" << c)<<std::endl;

					SparseFindPivot (LigneA[(size_t)k], Rank, c, col_density, determinant);
					//                     LigneA.write(std::cerr << "PIV, k:" << k << ", Rank:" << Rank << ", c:" << c)<<std::endl;
					if (c != -1) {
						for (l = (unsigned long)k + 1; l < (unsigned long)Ni; ++l)
							eliminate (LigneA[(size_t)l], LigneA[(size_t)k], Rank, c, col_density);
					}

					//                     LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
					nbelem += LigneA[(size_t)k].size ();
#endif
					LigneA[(size_t)k] = Vzer;
				}

			}//for k

			SparseFindPivot (LigneA[(size_t)last], Rank, c, determinant);

#ifdef __LINBOX_COUNT__
			nbelem += LigneA[(size_t)last].size ();
			commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
			<< "Left elements : " << nbelem << std::endl;
#endif

#ifdef __LINBOX_FILLIN__
			long sl(0);
			for (size_t l=0; l < Ni; ++l)
				sl += (long)LigneA[(size_t)l].size ();

			commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
			<< "Fillin (" << Rank << "/" << Ni << ") = " << sl
			<< std::endl;
#endif

			integer card;

			if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
				field().init(determinant,0UL);

			field().write(commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
				 << "Determinant : ", determinant)
			<< " over GF (" << field().cardinality (card) << ")" << std::endl;

			commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
			<< "Rank : " << Rank
			<< " over GF (" << card << ")" << std::endl;
			commentator().stop ("done", 0, "IPLR");
			return Rank;
		}

		template <class _Field>
		template <class Matrix, class Perm> inline unsigned long&
		GaussDomain<_Field>::InPlaceLinearPivoting (unsigned long &Rank,
							    Element        &determinant,
							    Matrix         &LigneA,
							    Perm           &P,
							    unsigned long   Ni,
							    unsigned long   Nj) const
		{
			typedef typename Matrix::Row        Vector;
			typedef typename Vector::value_type E;

			// Requirements : LigneA is an array of sparse rows
			// In place (LigneA is modified)
			// With reordering (D is a density type. Density is allocated here)
			//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
			commentator().start ("Gaussian elimination with reordering",
					   "IPLR", Ni);
			field().write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				  << "Gaussian elimination on " << Ni << " x " << Nj << " matrix, over: ") << std::endl;

#ifdef __LINBOX_COUNT__
			long long nbelem = 0;
#endif

			field().init(determinant,1UL);
			// allocation of the column density
			std::vector<size_t> col_density (Nj);

			// assignment of LigneA with the domain object
			for (unsigned long jj = 0; jj < Ni; ++jj)
				for (unsigned long k = 0; k < LigneA[(size_t)jj].size (); k++)
					++col_density[LigneA[(size_t)jj][k].first];

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
			for (long k = 0; k < last; ++k) {
				unsigned long l;
				long p = k, s =(long) LigneA[(size_t)k].size (), sl;

#ifdef __LINBOX_FILLIN__
				if ( ! (k % 100) )
#else
				if ( ! (k % sstep) )
#endif
				{
					commentator().progress (k);
#ifdef __LINBOX_FILLIN__
					for (sl = 0, l = 0; l < Ni; ++l)
						sl +=(long) LigneA[(size_t)l].size ();

					commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
					<< "Fillin (" << Rank << "/" << Ni << ") = "
					<< sl
					<< " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
					<< double(sl)/double(Ni-k) << " avg)"
					<< std::endl;
#endif
				}

				if (s) {
					// Row permutation for the sparsest row
					for (l = (unsigned long)k + 1; l < (unsigned long)Ni; ++l)
						if (((sl =(long) LigneA[(size_t)l].size ()) < s) && (sl)) {
							s = sl;
							p = (long)l;
						}

					if (p != k) {
						field().negin(determinant);
						Vector vtm = LigneA[(size_t)k];
						LigneA[(size_t)k] = LigneA[(size_t)p];
						LigneA[(size_t)p] = vtm;
					}

					//                     LigneA.write(std::cerr << "BEF, k:" << k << ", Rank:" << Rank << ", c:" << c)<<std::endl;

					SparseFindPivot (LigneA[(size_t)k], Rank, c, col_density, determinant);
					//                     LigneA.write(std::cerr << "PIV, k:" << k << ", Rank:" << Rank << ", c:" << c)<<std::endl;
					if (c != -1) {
						if ( c != (static_cast<long>(Rank)-1) )
							P.permute(Rank-1,(size_t)c);
						for (long ll=0; ll < k ; ++ll)
							permute( LigneA[(size_t)ll], Rank, c);

						for (l = (unsigned long)k + 1; l < (unsigned long)Ni; ++l)
							eliminate (LigneA[(size_t)l], LigneA[(size_t)k], Rank, c, col_density);
					}

					//                     LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
					nbelem += LigneA[(size_t)k].size ();
#endif
				}

			}//for k

			SparseFindPivot (LigneA[(size_t)last], Rank, c, determinant);
			if ( (c != -1) && (c != (static_cast<long>(Rank)-1) ) ) {
				P.permute(Rank-1,(size_t)c);
				for (long ll=0; ll < last ; ++ll)
					permute( LigneA[(size_t)ll], Rank, c);
			}


#ifdef __LINBOX_COUNT__
			nbelem += LigneA[(size_t)last].size ();
			commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
			<< "Left elements : " << nbelem << std::endl;
#endif

#ifdef __LINBOX_FILLIN__
			long sl(0);
			for (size_t l=0; l < Ni; ++l)
				sl +=(long) LigneA[(size_t)l].size ();

			commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
			<< "Fillin (" << Rank << "/" << Ni << ") = " << sl
			<< std::endl;
#endif

			integer card;

			if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
				field().init(determinant,0UL);

			field().write(commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
				 << "Determinant : ", determinant)
			<< " over GF (" << field().cardinality (card) << ")" << std::endl;

			commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
			<< "Rank : " << Rank
			<< " over GF (" << card << ")" << std::endl;
			commentator().stop ("done", 0, "IPLR");
			return Rank;
		}

		template <class _Field>
		template <class Matrix> inline unsigned long&
		GaussDomain<_Field>::NoReordering (unsigned long &res,
						   Element       &determinant,
						   Matrix        &LigneA,
						   unsigned long  Ni,
						   unsigned long  Nj) const
		{
			// Requirements : SLA is an array of sparse rows
			// IN PLACE.
			// Without reordering (Pivot is first non-zero in row)
			//     long Ni = SLA.n_row (), Nj = SLA.n_col ();
			//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
			commentator().start ("Gaussian elimination (no reordering)",
					   "NoRe", Ni);
			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Gaussian elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

			typedef typename Matrix::Row          Vector;
			typedef typename Vector::value_type   E;
			// typedef typename Matrix::Element      Elem;

#ifdef __LINBOX_COUNT__
			long long nbelem = 0;
#endif
			Vector Vzer (0);

			field().init(determinant,1UL);
			long last = (long)Ni - 1;
			long c;
			unsigned long indcol (0);

			for (long k = 0; k < last; ++k) {
				if (!(k % 1000))
					commentator().progress (k);

				unsigned long l;

				if (!LigneA[(size_t)k].empty ()) {
					SparseFindPivot (LigneA[(size_t)k], indcol, c, determinant);
					if (c !=  -1)
						for (l = (unsigned long)k + 1; l < (unsigned long)Ni; ++l)
							eliminate (LigneA[(size_t)l], LigneA[(size_t)k], indcol, c);

#ifdef __LINBOX_COUNT__
					nbelem += LigneA[(size_t)k].size ();
#endif
					LigneA[(size_t)k] = Vzer;
				}
			}

			SparseFindPivot ( LigneA[(size_t)last], indcol, c, determinant);

#ifdef __LINBOX_COUNT__
			nbelem += LigneA[(size_t)last].size ();
			commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
			<< "Left elements : " << nbelem << std::endl;
#endif

			res = indcol;

			if ((res < Ni) || (res < Nj))
				if ((res < Ni) || (res < Nj) || (Ni == 0) || (Nj == 0))
					field().init(determinant,0UL);

			integer card;

			field().write(commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
				 << "Determinant : ", determinant)
			<< " over GF (" << field().cardinality (card) << ")" << std::endl;

			commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
			<< "Rank : " << res
			<< " over GF (" << card << ")" << std::endl;
			commentator().stop ("done", 0, "NoRe");
			return res;
		}


		template <class _Field>
		template<class Vector> inline void
		GaussDomain<_Field>::Upper (Vector        &lignecur,
					    const Vector  &lignepivot,
					    unsigned long  indcol,
					    long  indpermut) const
		{
			static typename _Field::Element zero = field().init(zero);

			long n = lignecur.size () ;
			long k = (long) indcol - 1 ;

			// permutation if one has been performed to compute the pivot
			if (indpermut != k) {
				typename Vector::value_type tmp = lignecur[k];
				lignecur[k] = lignecur[indpermut];
				lignecur[indpermut] = tmp;
			}

			typename Vector::value_type headcoeff;
			field().divin (field().neg (headcoeff, lignecur[k]), lignepivot[k]);



			// LU in place
			field().assign (lignecur[k], zero);
			for (long j = k; ++j < n;)
				field().axpyin (lignecur[j], headcoeff, lignepivot[j]) ;
		}

		template <class _Field>
		template <class Vector> inline void
		GaussDomain<_Field>::LU (Vector        &lignecur,
					 const Vector  &lignepivot,
					 unsigned long  indcol,
					 long  indpermut) const
		{
			long n = lignecur.size ();
			long k = (long) indcol - 1;

			// permutation if one has been performed to compute the pivot
			if (indpermut != k) {
				typename Vector::value_type tmp = lignecur[k];
				lignecur[k] = lignecur[indpermut];
				lignecur[indpermut] = tmp;
			}

			typename Vector::value_type headcoeff;
			// LU in place
			field().div (headcoeff, lignecur[k], lignepivot[k]);
			field().assign (lignecur[k], headcoeff);
			field().negin (headcoeff);
			for (long j = k; ++j < n;)
				field().axpyin (lignecur[j],headcoeff,lignepivot[j]);
		}


		template <class _Field>
		template <class Matrix> inline unsigned long &
		GaussDomain<_Field>::upperin (unsigned long &res, Matrix &A) const
		{
			// Requirements : A is an array of rows
			// In place (A is modified)
			// Without reordering (Pivot is first non-zero in row)
			long Ni = A.rowdim ();
			long last = Ni - 1;
			long c;
			unsigned long indcol = 0;

			for (long k = 0; k < last; ++k) {
				FindPivot (A[k], indcol, c);
				if (c != -1)
					for (long l = k + 1; l < Ni; ++l)
						Upper (A[l], A[k], indcol, c);
			}

			FindPivot (A[last], indcol, c);
			return res = indcol;
		}

		template <class _Field>
		template <class Matrix> inline unsigned long &
		GaussDomain<_Field>::LUin (unsigned long &res, Matrix &A) const
		{
			// Requirements : A is an array of rows
			// In place (A is modified)
			// Without reordering (Pivot is first non-zero in row)

			long Ni = A.rowdim ();
			long last = Ni - 1;
			long c;
			unsigned long indcol = 0;

			for (long k = 0; k < last; ++k) {
				FindPivot (A[k], indcol, c);
				if (c != -1)
					for (long l = k + 1; l < Ni; ++l)
						LU (A[l], A[k], indcol, c);
			}

			FindPivot (A[last], indcol, c);
			return res = indcol;
		}


	} // namespace LinBox

#endif // __LINBOX_gauss_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

