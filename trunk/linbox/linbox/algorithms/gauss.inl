/* linbox/algorithms/gauss.inl
 * Copyright (C) 1999 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <11 Sep 08 14:57:31 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 */
// ========================================================================= //
// (C) The Linbox Group 1999
// Calcul de rang par la m�thode de Gauss pivot par ligne, sur matrice creuse
// ========================================================================= //

#ifndef __GAUSS_INL
#define __GAUSS_INL

#include "linbox/algorithms/gauss.h"
#include "linbox/util/commentator.h"
#include <utility>

#ifdef __LINBOX_ALL__
#define __LINBOX_COUNT__
#define __LINBOX_OFTEN__
#define __LINBOX_FILLIN__
#endif

namespace LinBox 
{
    template <class _Field>
    template <class Matrix, class Perm> unsigned long& 
    GaussDomain<_Field>::QLUPin (unsigned long &rank,
                                 Element       &determinant,
                                 Perm          &Q,
                                 Matrix        &LigneL,
                                 Matrix        &LigneA,
                                 Perm          &P,
                                 unsigned long Ni,
                                 unsigned long Nj)
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
	commentator.start ("Gaussian elimination with reordering",
			   "IPLR", Ni);
	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
            << "Gaussian QLUP elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

#ifdef __LINBOX_COUNT__
	long long nbelem = 0;
#endif

        Element Eone; _F.init(Eone,1UL);
        _F.init(determinant,1UL);
            // allocation of the column density
        std::vector<size_t> col_density (Nj);


        for(typename Matrix::RowIterator LigneL_it = LigneL.rowBegin() ;
            LigneL_it != LigneL.rowEnd(); ++LigneL_it)
            LigneL_it->reserve(16);

        std::deque<std::pair<size_t,size_t> > invQ;

            // assignment of LigneA with the domain object
	for (unsigned long jj = 0; jj < Ni; ++jj) 
            for (unsigned long k = 0; k < LigneA[jj].size (); k++)
                ++col_density[LigneA[jj][k].first];

	long last = Ni - 1;
	long c;
	rank = 0;

#ifdef __LINBOX_OFTEN__
	long sstep = last/40;
	if (sstep > 1000) sstep = 1000;
#else
	long sstep = 1000;
#endif
            // Elimination steps with reordering

        typename Matrix::RowIterator LigneA_k = LigneA.rowBegin(), LigneA_p;
	for (long k = 0; k < last; ++k, ++LigneA_k) {
            long p = k, s = 0;

#ifdef __LINBOX_FILLIN__  
            if ( ! (k % 100) ) {
#else          
            if ( ! (k % sstep) ) {
#endif
                commentator.progress (k);
#ifdef __LINBOX_FILLIN__   
                long sl;
                for (sl = 0, l = 0; l < Ni; ++l)
                    sl += LigneA[l].size ();
                
                commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
                    << "Fillin (" << rank << "/" << Ni << ") = "
                    << sl 
                    << " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
                    << double(sl)/double(Ni-k) << " avg)"
                    << std::endl;
#endif 
            }
            
            long l;
            for(l = k; l < static_cast<long>(Ni); ++l) {
                if ( (s = LigneA[l].size()) ) {
                    p = l;
                    break;
                }
            }
	    
            if (s) {
                long sl;
                    // Row permutation for the sparsest row
                for (; l < static_cast<long>(Ni); ++l)
                    if (((sl = LigneA[l].size ()) < s) && (sl)) {
                        s = sl;
                        p = l;
                    }
                
                if (p != k) {
//                         std::cerr << "Permuting rows: " << k << " <--> " << p << std::endl;
                    invQ.push_front( std::pair<size_t,size_t>(k,p) );
                    _F.negin(determinant);
                    std::swap( *LigneA_k, LigneA[p]);
                    std::swap( LigneL[k], LigneL[p]);
                }
                
                
                SparseFindPivot (*LigneA_k, rank, c, col_density, determinant);
                
                if (c != -1) {
                    long ll;
                    if ( c != (static_cast<long>(rank)-1) ) {
                        P.permute(rank-1,c);
                        for (ll=0      ; ll < k ; ++ll)
                            permute( LigneA[ll], rank, c);
                    }
                    long npiv=LigneA_k->size();
                    for (ll = k+1; ll < static_cast<long>(Ni); ++ll) {
                        E hc; hc.first=rank-1;
                        eliminate (hc.second, LigneA[ll], *LigneA_k, rank, c, npiv, col_density);
                        LigneL[ll].push_back(hc);
                    }
                }
                
//                     LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
                nbelem += LigneA_k->size ();
#endif
            }
            E one(k,Eone);
            LigneL[k].push_back(one);
//                 LigneL.write(rep << "L:= ", FORMAT_MAPLE) << std::endl;
//                 LigneA.write(rep << "U:= ", FORMAT_MAPLE) << std::endl;
	}//for k

	SparseFindPivot ( LigneA[last], rank, c, determinant);
            
        E one(last,Eone);
        LigneL[last].push_back(one);

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
            << "Fillin (" << rank << "/" << Ni << ") = " << sl 
            << std::endl;
#endif
        
        if ((rank < Ni) || (rank < Nj))
            _F.init(determinant,0UL);
        
        integer card;
        _F.write(commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
                 << "Determinant : ", determinant)
                 << " over GF (" << _F.cardinality (card) << ")" << std::endl;
        
        for(std::deque<std::pair<size_t,size_t> >::const_iterator it = invQ.begin(); it!=invQ.end();++it)
            Q.permute( it->first, it->second );
        
        
//             Q.write(rep << "Q:= ", FORMAT_MAPLE) << ':' << std::endl;
//             LigneL.write(rep << "L:= ", FORMAT_MAPLE) << ':' << std::endl;
//             LigneA.write(rep << "U:= ", FORMAT_MAPLE) << ':' << std::endl;
//             P.write(rep << "P:= ", FORMAT_MAPLE) << ':' << std::endl;
        
        commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
            << "Rank : " << rank
            << " over GF (" << card << ")" << std::endl;
        commentator.stop ("done", 0, "IPLR");
        


        return rank;
    }
        
    template <class _Field>
    template <class Matrix> unsigned long& 
    GaussDomain<_Field>::InPlaceLinearPivoting (unsigned long &res,
                                                Element        &determinant,
                                                Matrix         &LigneA,
                                                unsigned long   Ni,
                                                unsigned long   Nj)
    {
        typedef typename Matrix::Row        Vector;
        typedef typename Vector::value_type E;    
        
            // Requirements : LigneA is an array of sparse rows
            // In place (LigneA is modified)
            // With reordering (D is a density type. Density is allocated here)
            //    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
        commentator.start ("Gaussian elimination with reordering",
                           "IPLR", Ni);
        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
            << "Gaussian elimination on " << Ni << " x " << Nj << " matrix" << std::endl;
        
#ifdef __LINBOX_COUNT__
        long long nbelem = 0;
#endif

        _F.init(determinant,1UL);
        Vector Vzer (0);
            // allocation of the column density
        std::vector<size_t> col_density (Nj);        

            // assignment of LigneA with the domain object
        for (unsigned long jj = 0; jj < Ni; ++jj) 
            for (unsigned long k = 0; k < LigneA[jj].size (); k++)
                ++col_density[LigneA[jj][k].first];

        long last = Ni - 1;
        long c;
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
                    << sl 
                    << " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
                    << double(sl)/double(Ni-k) << " avg)"
                    << std::endl;
#endif 
            }
	    
            if (s) {
                    // Row permutation for the sparsest row
                for (l = k + 1; l < Ni; ++l)
                    if (((sl = LigneA[l].size ()) < s) && (sl)) {
                        s = sl;
                        p = l;
                    }

                if (p != k) {
                    _F.negin(determinant);
                    Vector vtm = LigneA[k];
                    LigneA[k] = LigneA[p];
                    LigneA[p] = vtm;
                }
		    
//                     LigneA.write(std::cerr << "BEF, k:" << k << ", indcol:" << indcol << ", c:" << c)<<std::endl;
                    
                SparseFindPivot (LigneA[k], indcol, c, col_density, determinant);
//                     LigneA.write(std::cerr << "PIV, k:" << k << ", indcol:" << indcol << ", c:" << c)<<std::endl;
                if (c != -1) {
                    for (l = k + 1; l < Ni; ++l) 
                        eliminate (LigneA[l], LigneA[k], indcol, c, col_density);
                }
                    
//                     LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
                nbelem += LigneA[k].size ();
#endif
                LigneA[k] = Vzer;
            }
	    	    
	}//for k

	SparseFindPivot (LigneA[last], indcol, c, determinant);
            
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
            << "Fillin (" << indcol << "/" << Ni << ") = " << sl 
            << std::endl;
#endif
            
        res = indcol;
            
        integer card;
            
        if ((res < Ni) || (res < Nj))
            _F.init(determinant,0UL);
        _F.write(commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
                 << "Determinant : ", determinant)
                 << " over GF (" << _F.cardinality (card) << ")" << std::endl;
            
        commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
            << "Rank : " << res
            << " over GF (" << card << ")" << std::endl;
        commentator.stop ("done", 0, "IPLR");
        return res;
    }

    template <class _Field>
    template <class Matrix> unsigned long& 
    GaussDomain<_Field>::NoReordering (unsigned long &res,
                                                       Element       &determinant,
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
                           "NoRe", Ni);
        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
            << "Gaussian elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

        typedef typename Matrix::Row          Vector;
        typedef typename Vector::value_type   E;    
        typedef typename Matrix::Element      Element;    
    
#ifdef __LINBOX_COUNT__
        long long nbelem = 0;
#endif
        Vector Vzer (0);
 
        _F.init(determinant,1UL);
        long last = Ni - 1;
        long c;
        unsigned long indcol (0);
    
        for (long k = 0; k < last; ++k) {
            if (!(k % 1000))
                commentator.progress (k);

            unsigned long l;

            if (!LigneA[k].empty ()) {
                SparseFindPivot (LigneA[k], indcol, c, determinant);
                if (c !=  -1)
                    for (l = k + 1; l < Ni; ++l)
                        eliminate (LigneA[l], LigneA[k], indcol, c);
                        
#ifdef __LINBOX_COUNT__
                nbelem += LigneA[k].size ();
#endif
                LigneA[k] = Vzer;
            }
        }

        SparseFindPivot ( LigneA[last], indcol, c, determinant);

#ifdef __LINBOX_COUNT__
        nbelem += LigneA[last].size ();
        commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
            << "Left elements : " << nbelem << std::endl;
#endif
    
        res = indcol;

        if ((res < Ni) || (res < Nj))
            _F.init(determinant,0UL);

        integer card;

        _F.write(commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
                 << "Determinant : ", determinant)
                 << " over GF (" << _F.cardinality (card) << ")" << std::endl;
                
        commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT) 
            << "Rank : " << res
            << " over GF (" << card << ")" << std::endl;
        commentator.stop ("done", 0, "NoRe");
        return res;
    }


    template <class _Field>
    template<class Vector> void 
    GaussDomain<_Field>::Upper (Vector        &lignecur,
                                const Vector  &lignepivot,
                                unsigned long  indcol,
                                long  indpermut)
    {
        static typename _Field::Element zero = _F.init(zero);
        
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
        _F.assign (lignecur[k], zero);
        for (long j = k; ++j < n;)
            _F.axpyin (lignecur[j], headcoeff, lignepivot[j]) ;
    }

    template <class _Field>
    template <class Vector> void 
    GaussDomain<_Field>::LU (Vector        &lignecur,
                             const Vector  &lignepivot,
                             unsigned long  indcol,
                             long  indpermut)
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
    template <class Matrix> unsigned long &
    GaussDomain<_Field>::upperin (unsigned long &res, Matrix &A)
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
    template <class Matrix> unsigned long &
    GaussDomain<_Field>::LUin (unsigned long &res, Matrix &A)
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

#endif // __GAUSS_INL
