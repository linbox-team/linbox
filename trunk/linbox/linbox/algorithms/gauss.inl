/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/gauss.inl
 * Copyright (C) 1999 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 * -----------------------------------------------------------
 * 2003-02-02  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Ported to new matrix archetype; update interface to meet current
 * standards. Rename gauss_foo as foo and gauss_Uin as upperin
 *
 * Move function definitions to gauss.inl
 * -----------------------------------------------------------
 *
 * See COPYING for license information.
 */

// ========================================================================= //
// (C) The Linbox Group 1999
// Calcul de rang par la méthode de Gauss pivot par ligne, sur matrice creuse
// Time-stamp: <03 Nov 00 19:19:06 Jean-Guillaume.Dumas@imag.fr> 
// ========================================================================= //

#ifndef __GAUSS_INL
#define __GAUSS_INL

#include "linbox/algorithms/gauss.h"

#ifdef __LINBOX_ALL__
#define __LINBOX_COUNT__
#define __LINBOX_OFTEN__
#define __LINBOX_FILLIN__
#endif

namespace LinBox 
{

template <class _Field>
template <class Vector, class D>
void GaussDomain<_Field>::eliminate (Vector              &lignecourante,
				     const Vector        &lignepivot,
				     const unsigned long &indcol,
				     const unsigned long &indpermut,
				     D                   &columns)
{
	typedef typename Vector::value_type E;

	unsigned long k = indcol - 1;
	unsigned long nj = lignecourante.size () ;
	if (nj > 0) {
		unsigned long j_head = 0;

		for (; j_head < nj; ++j_head)
			if (lignecourante[j_head].first >= indpermut) break;

		unsigned long bjh = j_head - 1;

		if ((j_head < nj) && (lignecourante[j_head].first == indpermut)) {
			// -------------------------------------------
			// Permutation
			if ((unsigned long) indpermut != k) {
				if (lignecourante[0].first == k) {
					// non zero  <--> non zero
					Element tmp;
					_F.assign (tmp, lignecourante[0].second);
					_F.assign (lignecourante[0].second,
						   lignecourante[j_head].second);
					_F.assign (lignecourante[j_head].second, tmp);
				} else {
					// zero <--> non zero
					E tmp = lignecourante[j_head];
					--columns[tmp.first];
					++columns[k];
					tmp.first = k;

					for (long l = j_head; l > 0; l--)
						lignecourante[l] = lignecourante[l-1];

					lignecourante[0] = tmp;
				} 
				j_head = 0;
			}
			// -------------------------------------------
			// Elimination
			unsigned long npiv = lignepivot.size ();
			Vector construit (nj + npiv);

			// construit : <-- j
			// courante  : <-- m
			// pivot     : <-- l
			unsigned long j = 0;
			unsigned long m = j_head + 1;

			// A[i,k] <-- - A[i,k] / A[k,k]
			Element headcoeff;
			_F.divin (_F.neg (headcoeff, lignecourante[j_head].second),
				  lignepivot[0].second);

			--columns[lignecourante[j_head].first];
        
			// if A[k,j]=0, then A[i,j] <-- A[i,j]
			while (j < j_head) {
				construit[j] = lignecourante[j];
				j++;
			}

			unsigned long j_piv;

			unsigned long l = 0;

			for (; l < npiv; l++)
				if (lignepivot[l].first > k) break;

			// for all j such that (j>k) and A[k,j]!=0
			while (l < npiv) {
				j_piv = lignepivot[l].first;

				// if A[k,j]=0, then A[i,j] <-- A[i,j]
				while ((m < nj) && (lignecourante[m].first < j_piv))
					construit[j++] = lignecourante[m++];

				// if A[i,j]!=0, then A[i,j] <-- A[i,j] - A[i,k]*A[k,j]
				if ((m < nj) && (lignecourante[m].first == j_piv)) {
					Element tmp;

					_F.axpy (tmp, headcoeff, lignepivot[l].second,
						 lignecourante[m].second);

					if (! _F.isZero (tmp)) {
						_F.assign (lignecourante[m].second, tmp);
						construit[j++] = lignecourante[m++];
					} else
						--columns[lignecourante[m++].first];
				} else {
					Element tmp;

					_F.mul (tmp, headcoeff, lignepivot[l].second);

					if (! _F.isZero (tmp)) {
						++columns[j_piv];
						construit[j++] = E (j_piv, tmp);
					}
				}

				l++;
			}
        
			// if A[k,j]=0, then A[i,j] <-- A[i,j] 
			while (m<nj)
				construit[j++] = lignecourante[m++];
        
			construit.resize (j);
			lignecourante = construit;
		} else {
			// -------------------------------------------
			// Permutation
			if (indpermut != k) {
				unsigned long l = 0;

				for (; l < nj; ++l)
					if (lignecourante[l].first >= k) break;

				if ((l < nj) && (lignecourante[l].first == k))  {
					// non zero <--> zero
					E tmp = lignecourante[l];
					--columns[tmp.first];
					++columns[indpermut];
					tmp.first = indpermut;

					for (; l < bjh; l++)
						lignecourante[l] = lignecourante[l + 1];

					lignecourante[bjh] = tmp;
				} // else
				// zero <--> zero
			}
		}
	}
}

template <class _Field>
template <class Vector>
void GaussDomain<_Field>::eliminate (Vector              &lignecourante,
				     const Vector        &lignepivot,
				     const unsigned long &indcol,
				     const unsigned long &indpermut)
{
	typedef typename Vector::value_type E;

	unsigned long k = indcol - 1;
	unsigned long nj = lignecourante.size () ;

	if (nj > 0) {
		unsigned long j_head = 0;

		for (; j_head < nj; ++j_head)
			if (lignecourante[j_head].first >= indpermut) break;

		unsigned long bjh = j_head - 1;

		if ((j_head < nj) && (lignecourante[j_head].first == indpermut)) {
			// -------------------------------------------
			// Permutation
			if (indpermut != k) {
				if (lignecourante[0].first == k) {     
					// non zero  <--> non zero
					Element tmp;
					_F.assign (tmp, lignecourante[0].second) ;
					_F.assign (lignecourante[0].second,
						   lignecourante[j_head].second);
					_F.assign (lignecourante[j_head].second, tmp);
				} else {
					// zero <--> non zero
					E tmp = lignecourante[j_head];
					tmp.first = k;
					for (long l = j_head; l > 0; l--)
						lignecourante[l] = lignecourante[l-1];
					lignecourante[0] = tmp;
				}

				j_head = 0;
			}
			// -------------------------------------------
			// Elimination
			unsigned long npiv = lignepivot.size ();
			Vector construit (nj + npiv);
			// construit : <-- j
			// courante  : <-- m
			// pivot     : <-- l
			unsigned long j = 0;
			unsigned long m = j_head + 1;

			// A[i,k] <-- - A[i,k] / A[k,k]

			Element headcoeff;
			_F.divin (_F.neg (headcoeff, lignecourante[j_head].second),
				  lignepivot[0].second);

			// if A[k,j]=0, then A[i,j] <-- A[i,j]
			while (j < j_head) {
				construit[j] = lignecourante[j];
				j++;
			}

			unsigned long j_piv;
			unsigned long l = 0;

			for (; l < npiv; l++)
				if (lignepivot[l].first > k) break;

			// for all j such that (j>k) and A[k,j]!=0
			while (l < npiv) {
				j_piv = lignepivot[l].first;

				// if A[k,j]=0, then A[i,j] <-- A[i,j]
				while ((m < nj) && (lignecourante[m].first < j_piv))
					construit[j++] = lignecourante[m++];

				// if A[i,j]!=0, then A[i,j] <-- A[i,j] - A[i,k]*A[k,j]
				if ((m < nj) && (lignecourante[m].first == j_piv)) {
					Element tmp;
					_F.axpy (tmp, headcoeff, lignepivot[l].second,
						 lignecourante[m].second);

					if (! _F.isZero (tmp)) {
						_F.assign (lignecourante[m].second, tmp);
						construit[j++] = lignecourante[m++];
					} else
						++m;
                    
				} else {
					Element tmp;
					_F.mul (tmp, headcoeff, lignepivot[l].second);
					if (! _F.isZero (tmp))
						construit[j++] = E (j_piv, tmp);
				}
				l++;
			}
        
			// if A[k,j]=0, then A[i,j] <-- A[i,j] 
			while (m < nj)
				construit[j++] = lignecourante[m++];

			construit.resize (j);
			lignecourante = construit;
		} else {
			// -------------------------------------------
			// Permutation
			if (indpermut != k) {
				unsigned long l = 0;

				for (; l < nj; ++l)
					if (lignecourante[l].first >= k) break;

				if ((l < nj) && (lignecourante[l].first == k))  {
					// non zero <--> zero
					E tmp = lignecourante[l];
					tmp.first = indpermut;

					for (; l < bjh; l++)
						lignecourante[l] = lignecourante[l + 1];

					lignecourante[bjh] = tmp;
				} // else
				// zero <--> zero
			}
		}
	}
}

template <class _Field>
template<class Vector>
void GaussDomain<_Field>::Upper (Vector        &lignecur,
				 const Vector  &lignepivot,
				 unsigned long  indcol,
				 unsigned long  indpermut)
{
	long n = lignecur.size () ;
	long k = indcol - 1 ;

	// permutation if one has been performed to compute the pivot
	if (indpermut != k) {
		typename Vector::value_type tmp = lignecur[k];
		lignecur[k] = lignecur[indpermut];
		lignecur[indpermut] = tmp;
	}

	typename Vector::value_type headcoeff;
	_F.divin (_F.neg (headcoeff, lignecur[k]), lignepivot[k]);

	// LU in place
	_F.assign (lignecur[k], _F.zero);
	for (long j = k; ++j < n;)
		_F.axpyin (lignecur[j], headcoeff, lignepivot[j]) ;
}

template <class _Field>
template <class Vector>
void GaussDomain<_Field>::LU (Vector        &lignecur,
			      const Vector  &lignepivot,
			      unsigned long  indcol,
			      unsigned long  indpermut)
{
	long n = lignecur.size ();
	long k = indcol - 1;

	// permutation if one has been performed to compute the pivot
	if (indpermut != k) {
		typename Vector::value_type tmp = lignecur[k];
		lignecur[k] = lignecur[indpermut];
		lignecur[indpermut] = tmp;
	}

	typename Vector::value_type headcoeff;
	// LU in place
	_F.div (headcoeff, lignecur[k], lignepivot[k]);
	_F.assign (lignecur[k], headcoeff);
	_F.negin (headcoeff);
	for (long j = k; ++j < n;)
		_F.axpyin (lignecur[j],headcoeff,lignepivot[j]);
}

template <class _Field>
template <class Vector, class D>
void GaussDomain<_Field>::SparseFindPivot (Vector        &lignepivot,
					   unsigned long &indcol,
					   unsigned long &indpermut,
					   D             &columns)
{
	typedef typename Vector::value_type E;    

	long nj =  lignepivot.size ();

	if (nj > 0) {
		indpermut = lignepivot[0].first;

		long ds = --columns[indpermut], dl, p = 0;

		for (long j = 1; j < nj; ++j) {
			if ((dl = --columns[lignepivot[j].first]) < ds) {
				ds = dl;
				p = j;
			}
		}

		if (p != 0) {
			if (indpermut == indcol) {
				Element ttm;
				_F.assign (ttm, lignepivot[p].second);
				indpermut = lignepivot[p].first;
				_F.assign (lignepivot[p].second, lignepivot[0].second);
				_F.assign (lignepivot[0].second, ttm);
			} else {
				E ttm = lignepivot[p];
				indpermut = ttm.first;

				for (long m = p; m; --m)
					lignepivot[m] = lignepivot[m-1];

				lignepivot[0] = ttm;
			}
		}

		if (indpermut != indcol)
			// no need to decrement/increment, already done during the search
			lignepivot[0].first = indcol;

		indcol++ ;
	} else
		indpermut = -1;
}

template <class _Field>
template <class Vector>
void GaussDomain<_Field>::SparseFindPivot (Vector &lignepivot, unsigned long &indcol, unsigned long &indpermut)
{
	long nj = lignepivot.size ();

	if (nj > 0) {
		indpermut = lignepivot[0].first;
		if (indpermut != indcol)
			lignepivot[0].first = indcol;
		++indcol;
	} else
		indpermut = -1;
}

template <class _Field>
template <class Vector>
void GaussDomain<_Field>::FindPivot (Vector &lignepivot, unsigned long &k, unsigned long &indpermut)
{
	long n = lignepivot.size ();
	long j = k;

	for (; j < n ; ++j )
		if (!_F.isZero (lignepivot[j])) break ;

	if (j == n )
		indpermut = -1 ;
	else {
		indpermut = j ;
		if (indpermut != k) {
			typename Vector::value_type tmp = lignepivot[k] ;
			lignepivot[k] = lignepivot[j] ;
			lignepivot[j] = tmp ;
		}

		++k;
	}
}

template <class _Field>
template <class Matrix>
void GaussDomain<_Field>::rankinFullPivot (unsigned long &res,
				  Matrix        &LigneA,
				  bool           storrows)
{
	rankinFullPivot (res, LigneA, LigneA.rowdim (), LigneA.coldim (), storrows);
}

template <class _Field>
template <class Matrix>
void GaussDomain<_Field>::rankinFullPivot (unsigned long &res,
				  Matrix        &LigneA,
				  unsigned long  Ni,
				  unsigned long  Nj,
				  bool           storrows)
{
	typedef typename Matrix::Row        Vector;
	typedef typename Vector::value_type E;    

	// Requirements : LigneA is an array of sparse rows
	// In place (LigneA is modified)
	// With reordering (D is a density type. Density is allocated here)
	//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
	commentator.start ("Gaussian elimination with reordering",
			   "GaussDomain::gauss_rankin", Ni);
	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Gaussian elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

#ifdef __LINBOX_COUNT__
	long long nbelem = 0;
#endif

	Vector Vzer (0);
	if (storrows) Vzer.resize (1);

	// allocation of the column density
	std::vector<size_t> col_density (Nj);

	// assignment of LigneA with the domain object
	for (unsigned long jj = 0; jj < Ni; jj++) {
		Vector tmp = LigneA[jj];

		for (unsigned long k = 0; k < tmp.size (); k++)
			++col_density[tmp[k].first];
	}

	long last = Ni - 1;
	unsigned long c;
	unsigned long indcol = 0;

#ifdef __LINBOX_OFTEN__
	long sstep = last/40;
	if (sstep > 1000) sstep = 1000;
#else
	long sstep = 1000;
#endif
    
	// Elimination steps with reordering
	for (long k = 0; k < last; ++k) {
		unsigned long l;
		long p = k, s = LigneA[k].size (), sl;

#ifdef __LINBOX_FILLIN__  
		if ( ! (k % 100) ) {
#else          
		if ( ! (k % sstep) ) {
#endif
			commentator.progress (k);
#ifdef __LINBOX_FILLIN__            
			for (sl = 0, l = 0; l < Ni; ++l)
				sl += LigneA[l].size ();

			commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
				<< "Fillin (" << indcol << "/" << Ni << ") = "
				<< sl << std::endl;
#endif 
		}

		if (s) {
			for (l = k + 1; l < Ni; ++l)
				if (((sl = LigneA[l].size ()) < s) && (sl)) {
					s = sl;
					p = l;
				}

			if (p != k) {
				Vector vtm = LigneA[k];
				LigneA[k] = LigneA[p];
				LigneA[p] = vtm;
			}

			SparseFindPivot (LigneA[k], indcol, c, col_density);

			if (c != (unsigned long) -1)
				for (l = k + 1; l < Ni; ++l)
					eliminate (LigneA[l], LigneA[k],
						   indcol, c, col_density);
#ifdef __LINBOX_COUNT__
			nbelem += LigneA[k].size ();
#endif
			LigneA[k] = Vzer;
		}
	}

	SparseFindPivot (LigneA[last], indcol, c);

#ifdef __LINBOX_COUNT__
	nbelem += LigneA[last].size ();
	commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
		<< "Left elements : " << nbelem << std::endl;
#endif

#ifdef __LINBOX_FILLIN__  
	long sl = 0, l = 0;
	for (; l < Ni; ++l)
		sl += LigneA[l].size ();

	commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
		<< "Fillin (" << indcol << "/" << Ni << ") = " << sl << std::endl;
#endif
    
	res = indcol;

	integer card;

	commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
		<< "Rank : " << res
		<< " over GF (" << _F.cardinality (card) << ")" << std::endl;
	commentator.stop ("done", 0, "GaussDomain::gauss_rankin");
}

template <class _Field>
template <class Matrix>
void GaussDomain<_Field>::rank (unsigned long &res, const Matrix &SLA)
{
	long Ni = SLA.n_row (), Nj = SLA.n_col ();
	typedef typename Matrix::Row        Vector;
	typedef typename Vector::value_type E;    
	typedef typename Matrix::Element    Element;
    
	Vector *LigneA = new Vector[Ni];

	for (unsigned long jj = 0; jj < Ni; jj++) {
		Vector tmp = SLA[jj];
		Vector toto (tmp.size ());
		long rs = 0;
		long k = 0;

		for (; k < tmp.size (); k++) {
			Element r;

			_F.assign (r, tmp[k].second);
			if (! _F.isZero (r))
				toto[rs++] = E (tmp[k].first, r); 
		}

		toto.resize (rs);
		LigneA[jj] = toto;
	}

	long last = Ni-1;
	unsigned long c;
	unsigned long indcol = 0;
    
	for (long k = 0; k < last; ++k) {
		long l, p = k, s = LigneA[k].size (), sl;

		if (s) {
			SparseFindPivot (LigneA[k], indcol, c) ;

			if (c != -1)
				for (l = k + 1; l < Ni; ++l)
					eliminate (LigneA[l], LigneA[k], indcol, c);
		}
	}

	SparseFindPivot (LigneA[last], indcol, c);
    
	delete [] LigneA;

	res = indcol;
}

template <class _Field>
template <class Matrix>
void GaussDomain<_Field>::rankin (unsigned long &res, Matrix &LigneA)
{
	rankin (res, LigneA, LigneA.rowdim (), LigneA.coldim ());
}

template <class _Field>
template <class Matrix>
void GaussDomain<_Field>::rankin (unsigned long &res,
				  Matrix        &LigneA,
				  unsigned long  Ni,
				  unsigned long  Nj)
{
	// Requirements : SLA is an array of sparse rows
	// IN PLACE.
	// Without reordering (Pivot is first non-zero in row)
	//     long Ni = SLA.n_row (), Nj = SLA.n_col ();
	//    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
	commentator.start ("Gaussian elimination (no reordering)",
			   "GaussDomain::gauss_rankin", Ni);
	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
		<< "Gaussian elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

	typedef typename Matrix::Row          Vector;
	typedef typename Vector::value_type   E;    
	typedef typename Matrix::Element      Element;    
    
#ifdef __LINBOX_COUNT__
	long long nbelem = 0;
#endif
	Vector Vzer (0);
 

	long last = Ni - 1;
	unsigned long c;
	unsigned long indcol (0);
    
	for (long k = 0; k < last; ++k) {
		if (!(k % 1000))
			commentator.progress (k);

		unsigned long l;

		if (!LigneA[k].empty ()) {
			SparseFindPivot (LigneA[k], indcol, c);
			if (c != (unsigned long) -1)
				for (l = k + 1; l < Ni; ++l)
					eliminate (LigneA[l], LigneA[k], indcol, c);
#ifdef __LINBOX_COUNT__
			nbelem += LigneA[k].size ();
#endif
			LigneA[k] = Vzer;
		}
	}

	SparseFindPivot ( LigneA[last], indcol, c );

#ifdef __LINBOX_COUNT__
	nbelem += LigneA[last].size ();
	commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
		<< "Left elements : " << nbelem << std::endl;
#endif
    
	res = indcol;

	integer card;

	commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
		<< "Rank : " << res
		<< " over GF (" << _F.cardinality (card) << ")" << std::endl;
	commentator.stop ("done", 0, "GaussDomain::gauss_rankin");
}

template <class _Field>
template <class Matrix>
long &GaussDomain<_Field>::upperin (unsigned long &res, Matrix &A)
{
	// Requirements : A is an array of rows
	// In place (A is modified)
	// Without reordering (Pivot is first non-zero in row)
	long Ni = A.rowdim ();
	long last = Ni - 1;
	unsigned long c;
	unsigned long indcol = 0;

	for (long k = 0; k < last; ++k) {
		FindPivot (A[k], indcol, c);
		if (c != -1)
			for (long l = k + 1; l < Ni; ++l)
				Upper (A[l], A[k], indcol, c);
	}

	FindPivot (A[last], indcol, c);
	res = indcol;
}

template <class _Field>
template <class Matrix>
long &GaussDomain<_Field>::LUin (unsigned long &res, Matrix &A)
{
	// Requirements : A is an array of rows
	// In place (A is modified)
	// Without reordering (Pivot is first non-zero in row)

	long Ni = A.rowdim ();
	long last = Ni - 1;
	unsigned long c;
	unsigned long indcol = 0;

	for (long k = 0; k < last; ++k) {
		FindPivot (A[k], indcol, c);
		if (c != -1)
			for (long l = k + 1; l < Ni; ++l)
				LU (A[l], A[k], indcol, c);
	}

	FindPivot (A[last], indcol, c);
	res = indcol;
}

} // namespace LinBox

#endif // __GAUSS_INL
