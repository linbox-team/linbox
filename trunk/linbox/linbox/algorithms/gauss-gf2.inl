/* linbox/algorithms/gauss-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <15 Jun 10 16:20:16 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
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
        inline unsigned long& 
        GaussDomain<GF2>::InPlaceLinearPivoting (unsigned long &rank,
                                  bool          &determinant,
                                  SparseSeqMatrix        &LigneA,
                                  Perm           &P,
                                  unsigned long Ni,
                                  unsigned long Nj) const
    {
            // Requirements : LigneA is an array of sparse rows
            // In place (LigneA is modified)
            // With reordering (D is a density type. Density is allocated here)
            //    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
	commentator.start ("Gaussian elimination with reordering over GF2",
			   "IPLRGF2", Ni);
	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
            << "Gaussian QLUP elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

#ifdef __LINBOX_COUNT__
	long long nbelem = 0;
#endif

        determinant = true;
            // allocation of the column density
        std::vector<size_t> col_density (Nj);


            // assignment of LigneA with the domain object
	for (unsigned long jj = 0; jj < Ni; ++jj) 
            for (unsigned long k = 0; k < LigneA[jj].size (); k++)
                ++col_density[LigneA[jj][k]];

	long last = Ni - 1;
	long c;
	rank = 0;

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
            if ( ! (k % 100) ) {
#else          
            if ( ! (k % sstep) ) {
#endif
                commentator.progress (k);
#ifdef __LINBOX_FILLIN__   
                long sl(0);
                for (size_t l = 0; l < Ni; ++l)
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
                if ( (s = LigneA[l].size()) != 0 ) {
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
                    std::swap( *LigneA_k, LigneA[p]);
                }
                
                
                SparseFindPivotBinary (*LigneA_k, rank, c, col_density, determinant);
                
                if (c != -1) {
                    long ll;
                    if ( c != (static_cast<long>(rank)-1) ) {
                        P.permute(rank-1,c);
                        for (ll=0      ; ll < k ; ++ll)
                            permuteBinary( LigneA[ll], rank, c);
                    }
                    long npiv=LigneA_k->size();
                    for (ll = k+1; ll < static_cast<long>(Ni); ++ll) {
                        bool elim=false;
                        eliminateBinary (elim, LigneA[ll], *LigneA_k, rank, c, npiv, col_density);
                    }
                }
                
//                     LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
                nbelem += LigneA_k->size ();
#endif
            }
//                 LigneA.write(rep << "U:= ", FORMAT_MAPLE) << std::endl;
	}//for k

	SparseFindPivotBinary ( LigneA[last], rank, c, determinant);
        if (c != -1) {
            if ( c != (static_cast<long>(rank)-1) ) {
                P.permute(rank-1,c);
                for (long ll=0      ; ll < last ; ++ll)
                    permuteBinary( LigneA[ll], rank, c);
            }
        }
            
#ifdef __LINBOX_COUNT__
        nbelem += LigneA[last].size ();
        commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
            << "Left elements : " << nbelem << std::endl;
#endif
        
#ifdef __LINBOX_FILLIN__  
        long sl(0);
        for (size_t l=0; l < Ni; ++l)
            sl += LigneA[l].size ();
        
        commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
            << "Fillin (" << rank << "/" << Ni << ") = " << sl 
            << std::endl;
#endif
        
        if ((rank < Ni) || (rank < Nj) || (Ni == 0) || (Nj == 0))
            determinant = false;
        
        integer card;
        commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT) 
            << "Determinant : " << determinant
            << " over GF (2)" << std::endl;
        
//             LigneA.write(rep << "U:= ", FORMAT_MAPLE) << ':' << std::endl;
        
        commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT) 
            << "Rank : " << rank
            << " over GF (2)" << std::endl;
        commentator.stop ("done", 0, "IPLRGF2");
        


        return rank;
    }
        
        // Specialization over GF2
        template <class SparseSeqMatrix, class Perm> inline unsigned long& 
        GaussDomain<GF2>::QLUPin (unsigned long &rank,
                                  bool          &determinant,
                                  Perm          &Q,
                                  SparseSeqMatrix        &LigneL,
                                  SparseSeqMatrix        &LigneA,
                                  Perm          &P,
                                  unsigned long Ni,
                                  unsigned long Nj) const
    {
        linbox_check( Q.coldim() == Q.rowdim() );
        linbox_check( P.coldim() == P.rowdim() );
        linbox_check( Q.coldim() == LigneL.size() );

	typedef typename SparseSeqMatrix::value_type Vector;
	typedef typename Vector::value_type E;    

            // Requirements : LigneA is an array of sparse rows
            // In place (LigneA is modified)
            // With reordering (D is a density type. Density is allocated here)
            //    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
	commentator.start ("Gaussian elimination with reordering over GF2",
			   "IPLRGF2", Ni);
	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
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
	for (unsigned long jj = 0; jj < Ni; ++jj) 
            for (unsigned long k = 0; k < LigneA[jj].size (); k++)
                ++col_density[LigneA[jj][k]];

	long last = Ni - 1;
	long c;
	rank = 0;

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
            if ( ! (k % 100) ) {
#else          
            if ( ! (k % sstep) ) {
#endif
                commentator.progress (k);
#ifdef __LINBOX_FILLIN__   
                long sl(0);
                for (size_t l = 0; l < Ni; ++l)
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
                    std::swap( *LigneA_k, LigneA[p]);
                    std::swap( LigneL[k], LigneL[p]);
                }
                
                
                SparseFindPivotBinary (*LigneA_k, rank, c, col_density, determinant);
                
                if (c != -1) {
                    long ll;
                    if ( c != (static_cast<long>(rank)-1) ) {
                        P.permute(rank-1,c);
                        for (ll=0      ; ll < k ; ++ll)
                            permuteBinary( LigneA[ll], rank, c);
                    }
                    long npiv=LigneA_k->size();
                    for (ll = k+1; ll < static_cast<long>(Ni); ++ll) {
                        E hc; hc=rank-1; bool elim=false;
                        eliminateBinary (elim, LigneA[ll], *LigneA_k, rank, c, npiv, col_density);
                        if(elim) LigneL[ll].push_back(hc);
                    }
                }
                
//                     LigneA.write(std::cerr << "AFT " )<<std::endl;
#ifdef __LINBOX_COUNT__
                nbelem += LigneA_k->size ();
#endif
            }
            LigneL[k].push_back(k);
 //                 LigneL.write(rep << "L:= ", FORMAT_MAPLE) << std::endl;
//                 LigneA.write(rep << "U:= ", FORMAT_MAPLE) << std::endl;
	}//for k

	SparseFindPivotBinary ( LigneA[last], rank, c, determinant);
        if (c != -1) {
            if ( c != (static_cast<long>(rank)-1) ) {
                P.permute(rank-1,c);
                for (long ll=0      ; ll < last ; ++ll)
                    permuteBinary( LigneA[ll], rank, c);
            }
        }
            
        LigneL[last].push_back(last);

#ifdef __LINBOX_COUNT__
        nbelem += LigneA[last].size ();
        commentator.report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
            << "Left elements : " << nbelem << std::endl;
#endif
        
#ifdef __LINBOX_FILLIN__  
        long sl(0);
        for (size_t l=0; l < Ni; ++l)
            sl += LigneA[l].size ();
        
        commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
            << "Fillin (" << rank << "/" << Ni << ") = " << sl 
            << std::endl;
#endif
        
        if ((rank < Ni) || (rank < Nj) || (Ni == 0) || (Nj == 0))
            determinant = false;
        
        integer card;
        commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT) 
            << "Determinant : " << determinant
            << " over GF (2)" << std::endl;
        
        for(std::deque<std::pair<size_t,size_t> >::const_iterator it = invQ.begin(); it!=invQ.end();++it)
            Q.permute( it->first, it->second );
        
        
//             Q.write(rep << "Q:= ", FORMAT_MAPLE) << ':' << std::endl;
//             LigneL.write(rep << "L:= ", FORMAT_MAPLE) << ':' << std::endl;
//             LigneA.write(rep << "U:= ", FORMAT_MAPLE) << ':' << std::endl;
//             P.write(rep << "P:= ", FORMAT_MAPLE) << ':' << std::endl;
        
        commentator.report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT) 
            << "Rank : " << rank
            << " over GF (2)" << std::endl;
        commentator.stop ("done", 0, "IPLRGF2");
        


        return rank;
    }


} // namespace LinBox

#endif // __LINBOX_gauss_gf2_INL
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
