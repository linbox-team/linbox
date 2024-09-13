/* linbox/algorithms/gauss.inl
 * Copyright (C) 1999 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <09 Feb 24 14:31:48 Jean-Guillaume.Dumas@imag.fr>
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
#include <givaro/zring.h>
#include <givaro/ring-interface.h>
#include <utility>
#include <type_traits>

#ifdef __LINBOX_ALL__
#define __LINBOX_COUNT__
#define __LINBOX_OFTEN__ __LINBOX_ALL__ // BB: ???
#define __LINBOX_FILLIN__
#endif

#ifdef __LINBOX_SpD_SWITCH__
#include <linbox/matrix/dense-matrix.h>
#include <numeric>
#  ifndef __LINBOX_SpD_MAXSPARSITY__
// Sparsity less than 1% --> switch to dense
#  define __LINBOX_SpD_MAXSPARSITY__ 0.01
#  endif
#endif

namespace LinBox
{
    template <class _Field>
    template <class _Matrix, class Perm> inline size_t&
    GaussDomain<_Field>::QLUPin (size_t &Rank,
                     Element       &determinant,
                     Perm          &Q,
                     _Matrix        &LigneL,
                     _Matrix        &LigneA,
                     Perm          &P,
                     size_t Ni,
                     size_t Nj) const
    {
        linbox_check( Q.coldim() == Q.rowdim() );
        linbox_check( P.coldim() == P.rowdim() );
        linbox_check( LigneL.coldim() == LigneL.rowdim() );
        linbox_check( Q.coldim() == LigneL.rowdim() );
        linbox_check( LigneL.coldim() == LigneA.rowdim() );
        linbox_check( LigneA.coldim() == P.rowdim() );

        typedef typename _Matrix::Row        Vector;
        typedef typename Vector::value_type E;

        // Requirements : LigneA is an array of sparse rows
        // In place (LigneA is modified)
        // With reordering (D is a density type. Density is allocated here)
        //    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
        commentator().start ("QLUPin Gaussian elimination with reordering",
                     "IPLR", Ni);
        commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
        << "Gaussian QLUP elimination on " << Ni << " x " << Nj << " matrix" << std::endl;

#ifdef __LINBOX_COUNT__
        long long nbelem = 0;
#endif

        field().assign(determinant,field().one);
        // allocation of the column density
        std::vector<size_t> col_density (Nj);


        for(typename _Matrix::RowIterator LigneL_it = LigneL.rowBegin() ;
            LigneL_it != LigneL.rowEnd(); ++LigneL_it)
            LigneL_it->reserve(16);

        std::deque<std::pair<size_t,size_t> > invQ;

        // assignment of LigneA with the domain object
        for (size_t jj = 0; jj < Ni; ++jj)
            for (size_t k = 0; k < LigneA[(size_t)jj].size (); k++)
                ++col_density[LigneA[(size_t)jj][k].first];

        const long last = (long)Ni - 1;
        long c;
        Rank = 0;
        bool degeneratedense=false;

#ifdef __LINBOX_OFTEN__
        long sstep = last/40;
        if (sstep > __LINBOX_OFTEN__) sstep = __LINBOX_OFTEN__;
        if (sstep <= 0) sstep = 1;
#else
#  ifdef __LINBOX_FILLIN__
        long sstep = 100;
#  else
        long sstep = 1000;
#  endif
#endif

        // Elimination steps with reordering

        typename _Matrix::RowIterator LigneA_k = LigneA.rowBegin();
        for (long k = 0; k < last; ++k, ++LigneA_k) {

#ifdef __LINBOX_SpD_MAXSPARSITY__

            if (accumulate(col_density.begin()+Rank,col_density.end(),0)>(Ni*Nj*__LINBOX_SpD_MAXSPARSITY__)) {
#  ifdef _LB_DEBUG
                std::cerr << "Dense switch: " << accumulate(col_density.begin()+Rank,col_density.end(),0) << '>' << Ni << '*' << Nj << '*' << __LINBOX_SpD_MAXSPARSITY__ << std::endl;
#  endif
                degeneratedense=std::is_base_of<Givaro::FiniteRingInterface<Element>,_Field>::value && true; break;
            }
#endif

            long p = k, s = 0;

            if ( ! (k % sstep) )
            {
                    commentator().progress (k);
#ifdef __LINBOX_FILLIN__
                    long sl(0);
                    for (size_t l = 0; l < Ni; ++l)
                        sl +=(long) LigneA[(size_t)l].size ();

                    commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
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
                // Row permutation for the sparsest row
                for (; l < static_cast<long>(Ni); ++l) {
                    long sl;
                    if (((sl =(long) LigneA[(size_t)l].size ()) < s) && (sl)) {
                        s = sl;
                        p = l;
                    }
                }

                if (p != k) {
                    // std::cerr << "Permuting rows: " << k << " <--> " << p << std::endl;
                    invQ.emplace_front((size_t)k,(size_t)p);
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
            LigneL[(size_t)k].emplace_back((unsigned)k,field().one);


//             LigneL.write(std::cerr << "L:= ", Tag::FileFormat::Maple) << std::endl;
//             LigneA.write(std::cerr << "U:= ", Tag::FileFormat::Maple) << std::endl;

        }//for k

        Continuation<_Matrix,Perm,
            std::is_base_of<Givaro::FiniteRingInterface<Element>,_Field>::value>
            ()(*this, Rank,determinant,invQ,LigneL,LigneA,P,Ni,Nj,degeneratedense);

#ifdef __LINBOX_COUNT__
        nbelem += LigneA[(size_t)(Ni-1)].size ();
        commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
        << "Left elements : " << nbelem << std::endl;
#endif

#ifdef __LINBOX_FILLIN__
        long sl(0);
        for (size_t l=0; l < Ni; ++l)
            sl += (long) LigneA[(size_t)l].size ();

        commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
        << "Fillin (" << Rank << "/" << Ni << ") = " << sl
        << std::endl;
#endif

        commentator().progress (Ni);

        if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
            field().assign(determinant,field().zero);

        field().write(field().write(commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
                  << "Determinant : ", determinant)
        << " over ") << std::endl;

        for(std::deque<std::pair<size_t,size_t> >::const_iterator it = invQ.begin(); it!=invQ.end();++it)
            Q.permute( it->first, it->second );


        std::ostream& rep = commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT);
        Q.write(rep << "Q:= ", Tag::FileFormat::Maple) << ':' << std::endl;
        LigneL.write(rep << "L:= ", Tag::FileFormat::Maple) << ':' << std::endl;
        LigneA.write(rep << "U:= ", Tag::FileFormat::Maple) << ':' << std::endl;
        P.write(rep << "P:= ", Tag::FileFormat::Maple) << ':' << std::endl;

        field().write(
            commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
            << "Rank : " << Rank
            << " over ") << std::endl;
        commentator().stop ("done", 0, "IPLR");



        return Rank;
    }



    template <class _Field>
    template <class _Matrix, class Perm> inline size_t&
    GaussDomain<_Field>::SparseContinuation (size_t &Rank,
                     Element       &determinant,
                     std::deque<std::pair<size_t,size_t> > &invQ,
                     _Matrix        &LigneL,
                     _Matrix        &LigneA,
                     Perm          &P,
                     size_t Ni,
                     size_t Nj) const
    {
        const long last = Ni-1;
        long c;
        SparseFindPivot ( LigneA[(size_t)last], Rank, c, determinant);
        if (c != -1) {
            if ( c != (static_cast<long>(Rank)-1) ) {
                P.permute(Rank-1,(size_t)c);
                for (long ll=0      ; ll < last ; ++ll)
                    permute( LigneA[(size_t)ll], Rank, c);
            }
        }

        LigneL[(size_t)last].emplace_back((unsigned)last,field().one);
        return Rank;
    }


    template <class _Field>
    template<class _Matrix, class Perm>
    struct GaussDomain<_Field>::Continuation<_Matrix,Perm,false> {
        size_t& operator()(
            const GaussDomain<_Field>& GD,
            size_t &Rank,
            typename GaussDomain<_Field>::Element       &determinant,
            std::deque<std::pair<size_t,size_t> > &invQ,
            _Matrix        &LigneL,
            _Matrix        &LigneA,
            Perm          &P,
            size_t Ni,
            size_t Nj,
            bool degeneratedense) const
            {
                return GD.SparseContinuation(Rank,determinant,invQ,LigneL,LigneA,P,Ni,Nj);
            }
    };

    template <class _Field>
    template <class _Matrix, class Perm>
    struct GaussDomain<_Field>::Continuation<_Matrix,Perm,true> {
        size_t& operator()(
            const GaussDomain<_Field>& GD,
            size_t &Rank,
            typename GaussDomain<_Field>::Element       &determinant,
            std::deque<std::pair<size_t,size_t> > &invQ,
            _Matrix        &LigneL,
            _Matrix        &LigneA,
            Perm          &P,
            size_t Ni,
            size_t Nj,
            bool degeneratedense) const
            {
                if (degeneratedense) {
                    return GD.DenseQLUPin(Rank,determinant,invQ,LigneL,LigneA,P,Ni,Nj);
                } else {
                    return GD.SparseContinuation(Rank,determinant,invQ,LigneL,LigneA,P,Ni,Nj);
                }
            }
    };






    template <class _Field>
    template <class _Matrix, class Perm> inline size_t&
    GaussDomain<_Field>::DenseQLUPin (size_t &Rank,
                     Element       &determinant,
                     std::deque<std::pair<size_t,size_t> > &dinvQ,
                     _Matrix        &dLigneL,
                     _Matrix        &dLigneA,
                     Perm          &dP,
                     size_t Ni,
                     size_t Nj) const
    {
        linbox_check( dP.coldim()      == dP.rowdim() );
        linbox_check( dLigneL.coldim() == dLigneL.rowdim() );
        linbox_check( dLigneL.coldim() == dLigneA.rowdim() );
        linbox_check( dLigneA.coldim() == dP.rowdim() );

//         std::deque<std::pair<size_t,size_t> > dinvQ(invQ);
//         _Matrix dLigneL(LigneL, this->field());
//         _Matrix dLigneA(LigneA, this->field());
//         Perm dP(P);

//         { Perm dQ(Ni);
//           for(std::deque<std::pair<size_t,size_t> >::const_iterator it = dinvQ.begin(); it!=dinvQ.end();++it)
//               dQ.permute( it->first, it->second );
//           dQ.write(std::cerr << "dQ:= ") << ';' << std::endl;
//         }
//         dLigneL.write(std::cerr << "dL:= ", Tag::FileFormat::Maple) << ':' << std::endl;
//         dLigneA.write(std::cerr << "dU:= ", Tag::FileFormat::Maple) << ':' << std::endl;
//         dP.write(std::cerr << "dP:= ") << ';' << std::endl;


        size_t sNi=Ni-Rank, sNj=Nj-Rank;
//         std::cerr << "Dense switch: " << sNi << 'x' << sNj << std::endl;
        BlasMatrix<_Field> A(this->field(), sNi, sNj);

        for(size_t di=Rank;di<Ni;++di) {
            for(size_t dj=0;dj<dLigneA[di].size();++dj)
                A.setEntry(di-Rank,dLigneA[di][dj].first-Rank, dLigneA[di][dj].second);
            dLigneA[di].resize(0);
        }


//         A.write(std::cerr << "A:= ", Tag::FileFormat::Maple) << std::endl;

        enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
        size_t *P2 = FFLAS::fflas_new<size_t>(sNi);
        size_t *Q2 = FFLAS::fflas_new<size_t>(sNj);
        for (size_t j=0;j<sNi;j++) P2[j]=0;
        for (size_t j=0;j<sNj;j++) Q2[j]=0;
        size_t R2 = FFPACK::PLUQ(this->field(), diag, sNi, sNj, A.getPointer(), sNj, P2, Q2);
//         std::cerr << "Rank:" << Rank << '+' << R2 << '=' << (Rank+R2) << std::endl;


//         std::cerr << dLigneL.rowdim() << 'x' << dLigneL.coldim() << std::endl;

//         { Perm dQ(Ni);
//           for(std::deque<std::pair<size_t,size_t> >::const_iterator it = dinvQ.begin(); it!=dinvQ.end();++it)
//               dQ.permute( it->first, it->second );
//           dQ.write(std::cerr << "dQ:= ") << ';' << std::endl;
//         }
//         dLigneL.write(std::cerr << "dL:= ", Tag::FileFormat::Maple) << ':' << std::endl;
//         dLigneA.write(std::cerr << "dU:= ", Tag::FileFormat::Maple) << ':' << std::endl;
//         dP.write(std::cerr << "dP:= ") << ';' << std::endl;


        Givaro::ZRing<long> Z;

//         std::cerr << '['; for (size_t j=0;j<sNi;j++)
//             std::cerr << P2[j] << ' ';
//         std::cerr << ']' << std::endl;
            // Left-Trans: P2^T * G
        for (int i=sNi;--i>=0;)
            if(i != static_cast<int>(P2[i])) {
                dinvQ.emplace_front( Rank+(size_t)i, Rank+P2[i] );
                this->field().negin(determinant);
                std::swap(dLigneL[i+Rank],dLigneL[P2[i]+Rank]);
            }

            // Put L2 in bottom right corner
        for(size_t i=0; i<sNi; ++i)
            for(size_t j=0; j<i; ++j)
                if (!this->field().isZero(A.getEntry(i,j)))
                    dLigneL[Rank+i].emplace_back(Rank+j,A.getEntry(i,j));
        for(size_t i=0; i<R2; ++i)
            dLigneL[Rank+i].emplace_back(Rank+i,this->field().one);


            // Put U2 in bottom right corner
        for(size_t i=0; i<sNi; ++i)
            for(size_t j=i; j<sNj; ++j)
                if (!this->field().isZero(A.getEntry(i,j)))
                    dLigneA[Rank+i].emplace_back(Rank+j,A.getEntry(i,j));
        for(size_t i=0; i<R2; ++i)
            this->field().mulin(determinant,A.getEntry(i,i));

//         std::cerr << '['; for (size_t j=0;j<sNj;++j)
//             std::cerr << Q2[j] << ' ';
//         std::cerr << ']' << std::endl;

            // Right-Trans: H * Q2^T
        for (size_t j=0;j<sNj;++j)
            if(j != Q2[j]) {
                this->field().negin(determinant);
                for (size_t l=0; l<Rank; ++l)
                    permute( dLigneA[l], j+Rank+1, Q2[j]+Rank);
            }

        FFPACK::applyP(Z,FFLAS::FflasLeft,FFLAS::FflasNoTrans,1,0,sNj,&(*(dP.getStorage().begin()))+Rank,1,Q2);

//         { Perm dQ(Ni);
//           for(std::deque<std::pair<size_t,size_t> >::const_iterator it = dinvQ.begin(); it!=dinvQ.end();++it)
//               dQ.permute( it->first, it->second );
//           dQ.write(std::cerr << "dQ:= ") << ';' << std::endl;
//         }
//         dLigneL.write(std::cerr << "dL:= ", Tag::FileFormat::Maple) << ':' << std::endl;
//         dLigneA.write(std::cerr << "dU:= ", Tag::FileFormat::Maple) << ':' << std::endl;
//         dP.write(std::cerr << "dP:= ") << ';' << std::endl;




        return Rank+=R2;
    }


    template <class _Field>
    template <class _Matrix> inline size_t&
    GaussDomain<_Field>::InPlaceLinearPivoting (size_t &Rank,
                            Element        &determinant,
                            _Matrix         &LigneA,
                            size_t   Ni,
                            size_t   Nj) const
    {
        typedef typename _Matrix::Row        Vector;

        // Requirements : LigneA is an array of sparse rows
        // In place (LigneA is modified)
        // With reordering (D is a density type. Density is allocated here)
        //    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
        commentator().start ("IPLR Gaussian elimination with reordering",
                     "IPLR", Ni);
        field().write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                   << "Gaussian elimination on " << Ni << " x " << Nj << " matrix, over: ") << std::endl;

#ifdef __LINBOX_COUNT__
        long long nbelem = 0;
#endif

        field().assign(determinant,field().one);
        Vector Vzer(0) ;

        // allocation of the column density
        std::vector<size_t> col_density (Nj);

        // assignment of LigneA with the domain object
        for (size_t jj = 0; jj < Ni; ++jj)
            for (size_t k = 0; k < LigneA[(size_t)jj].size (); k++)
                ++col_density[LigneA[(size_t)jj][k].first];

        const long last = (long)Ni - 1;
        long c;
        Rank = 0;

#ifdef __LINBOX_OFTEN__
        long sstep = last/40;
        if (sstep > __LINBOX_OFTEN__) sstep = __LINBOX_OFTEN__;
        if (sstep <= 0) sstep = 1;
#else
        long sstep = 1000;
#endif
        // Elimination steps with reordering
        for (long k = 0; k < last; ++k) {
            long p = k, s = (long)LigneA[(size_t)k].size ();

#ifdef __LINBOX_FILLIN__
            if ( ! (k % sstep) ) {
                size_t l;
                long sl;
                commentator().progress (k);
                for (sl = 0, l = 0; l < Ni; ++l)
                    sl += (long)LigneA[(size_t)l].size ();

                commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
                << "Fillin (" << Rank << "/" << Ni << ") = "
                << sl
                << " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
                << double(sl)/double(Ni-k) << " avg)"
                << std::endl;
            }
#else
            if ( ! (k % sstep) ) {
                commentator().progress (k);
            }
#endif

            if (s) {
                size_t l;
                // Row permutation for the sparsest row
                for (l = (size_t)k + 1; l < (size_t)Ni; ++l) {
                long sl;
                    if (((sl = (long)LigneA[(size_t)l].size ()) < s) && (sl)) {
                        s = sl;
                        p = (long)l;
                    }
                }

                if (p != k) {
                    field().negin(determinant);
//                     Vector vtm = LigneA[(size_t)k];
//                     LigneA[(size_t)k] = LigneA[(size_t)p];
//                     LigneA[(size_t)p] = vtm;
                    std::swap(LigneA[(size_t)k], LigneA[(size_t)p]);
                }

                //                     LigneA.write(std::cerr << "BEF, k:" << k << ", Rank:" << Rank << ", c:" << c)<<std::endl;

                SparseFindPivot (LigneA[(size_t)k], Rank, c, col_density, determinant);
                //                     LigneA.write(std::cerr << "PIV, k:" << k << ", Rank:" << Rank << ", c:" << c)<<std::endl;
                if (c != -1) {
                    for (l = (size_t)k + 1; l < (size_t)Ni; ++l)
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

        commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
        << "Fillin (" << Rank << "/" << Ni << ") = " << sl
        << std::endl;
#endif

        integer card;

        if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
            field().assign(determinant,field().zero);

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
    template <class _Matrix, class Perm> inline size_t&
    GaussDomain<_Field>::InPlaceLinearPivoting (size_t &Rank,
                            Element        &determinant,
                            _Matrix         &LigneA,
                            Perm           &P,
                            size_t   Ni,
                            size_t   Nj) const
    {
        typedef typename _Matrix::Row        Vector;

        // Requirements : LigneA is an array of sparse rows
        // In place (LigneA is modified)
        // With reordering (D is a density type. Density is allocated here)
        //    long Ni = LigneA.n_row (), Nj = LigneA.n_col ();
        commentator().start ("IPperm Gaussian elimination with reordering",
                     "IPLR", Ni);
        field().write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                   << "Gaussian elimination on " << Ni << " x " << Nj << " matrix, over: ") << std::endl;

#ifdef __LINBOX_COUNT__
        long long nbelem = 0;
#endif

        field().assign(determinant,field().one);
        // allocation of the column density
        std::vector<size_t> col_density (Nj);

        // assignment of LigneA with the domain object
        for (size_t jj = 0; jj < Ni; ++jj)
            for (size_t k = 0; k < LigneA[(size_t)jj].size (); k++)
                ++col_density[LigneA[(size_t)jj][k].first];

        const long last = (long)Ni - 1;
        long c;
        Rank = 0;

#ifdef __LINBOX_OFTEN__
        long sstep = last/40;
        if (sstep > __LINBOX_OFTEN__) sstep = __LINBOX_OFTEN__;
        if (sstep <= 0) sstep = 1;
#else
        long sstep = 1000;
#endif
        // Elimination steps with reordering
        for (long k = 0; k < last; ++k) {
            long p = k, s =(long) LigneA[(size_t)k].size ();

#ifdef __LINBOX_FILLIN__
            if ( ! (k % sstep) )
            {
                size_t l;
                long sl;
                commentator().progress (k);
                for (sl = 0, l = 0; l < Ni; ++l)
                    sl +=(long) LigneA[(size_t)l].size ();

                commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
                << "Fillin (" << Rank << "/" << Ni << ") = "
                << sl
                << " (" << double(sl)*100.0/double(Ni-k)/double(Nj-k) << "%, "
                << double(sl)/double(Ni-k) << " avg)"
                << std::endl;

            }
#else
            if ( ! (k % sstep) )
            {
                commentator().progress (k);
            }
#endif


            if (s) {
                size_t l;
                // Row permutation for the sparsest row
                for (l = (size_t)k + 1; l < (size_t)Ni; ++l) {
                long sl;
                    if (((sl =(long) LigneA[(size_t)l].size ()) < s) && (sl)) {
                        s = sl;
                        p = (long)l;
                    }
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

                    for (l = (size_t)k + 1; l < (size_t)Ni; ++l)
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

        commentator().report (Commentator::LEVEL_NORMAL, PARTIAL_RESULT)
        << "Fillin (" << Rank << "/" << Ni << ") = " << sl
        << std::endl;
#endif

        integer card;

        if ((Rank < Ni) || (Rank < Nj) || (Ni == 0) || (Nj == 0))
            field().assign(determinant,field().zero);

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
    template <class _Matrix> inline size_t&
    GaussDomain<_Field>::NoReordering (size_t &res,
                       Element       &determinant,
                       _Matrix        &LigneA,
                       size_t  Ni,
                       size_t  Nj) const
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

        typedef typename _Matrix::Row          Vector;
        // typedef typename Vector::value_type   E;
        // typedef typename _Matrix::Element      Elem;

#ifdef __LINBOX_COUNT__
        long long nbelem = 0;
#endif
        Vector Vzer (0);

        field().assign(determinant,field().one);
        const long last = (long)Ni - 1;
        long c;
        size_t indcol (0);

        for (long k = 0; k < last; ++k) {
            if (!(k % 1000))
                commentator().progress (k);


            if (!LigneA[(size_t)k].empty ()) {
                SparseFindPivot (LigneA[(size_t)k], indcol, c, determinant);
                if (c !=  -1) {
                    size_t l;
                    for (l = (size_t)k + 1; l < (size_t)Ni; ++l)
                        eliminate (LigneA[(size_t)l], LigneA[(size_t)k], indcol, c);
                }

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
                field().assign(determinant,field().zero);

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
                    size_t  indcol,
                    long  indpermut) const
    {

        long n = lignecur.size () ;
        long k = (long) indcol - 1 ;

        // permutation if one has been performed to compute the pivot
        if (indpermut != k) {
            std::swap(lignecur[k],lignecur[indpermut]);
        }

        typename Vector::value_type headcoeff;
        field().divin (field().neg (headcoeff, lignecur[k]), lignepivot[k]);



        // LU in place
        field().assign (lignecur[k], field().zero);
        for (long j = k; ++j < n;)
            field().axpyin (lignecur[j], headcoeff, lignepivot[j]) ;
    }

    template <class _Field>
    template <class Vector> inline void
    GaussDomain<_Field>::LU (Vector        &lignecur,
                 const Vector  &lignepivot,
                 size_t  indcol,
                 long  indpermut) const
    {
        long n = lignecur.size ();
        long k = (long) indcol - 1;

        // permutation if one has been performed to compute the pivot
        if (indpermut != k) {
           std::swap(lignecur[k],lignecur[indpermut]);
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
    template <class _Matrix> inline size_t &
    GaussDomain<_Field>::upperin (size_t &res, _Matrix &A) const
    {
        // Requirements : A is an array of rows
        // In place (A is modified)
        // Without reordering (Pivot is first non-zero in row)
        long Ni = A.rowdim ();
        const long last = Ni - 1;
        long c;
        size_t indcol = 0;

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
    template <class _Matrix> inline size_t &
    GaussDomain<_Field>::LUin (size_t &res, _Matrix &A) const
    {
        // Requirements : A is an array of rows
        // In place (A is modified)
        // Without reordering (Pivot is first non-zero in row)

        long Ni = A.rowdim ();
        const long last = Ni - 1;
        long c;
        size_t indcol = 0;

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

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
