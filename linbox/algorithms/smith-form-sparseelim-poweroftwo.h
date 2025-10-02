/* algorithms/smith-form-sparseelim-poweroftwo.h
 * Copyright (C) LinBox
 * Written by JG Dumas
 * Time-stamp: <11 Sep 25 17:27:18 Jean-Guillaume.Dumas@imag.fr>
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


#ifndef __LINBOX_pp_gauss_poweroftwo_H
#define __LINBOX_pp_gauss_poweroftwo_H
#include <map>
#include <givaro/givconfig.h> // for Signed_Trait
#include "linbox/algorithms/smith-form-sparseelim-local.h"

#ifdef LINBOX_DEBUG
#  ifndef LINBOX_pp_gauss_intermediate_OUT
#    define LINBOX_pp_gauss_intermediate_OUT
#  endif
#endif

// LINBOX_pp_gauss_intermediate_OUT outputs intermediate matrices
#ifdef LINBOX_pp_gauss_intermediate_OUT
#  ifndef LINBOX_pp_gauss_steps_OUT
#    define LINBOX_pp_gauss_steps_OUT
#  endif
#endif

// LINBOX_pp_gauss_steps_OUT outputs elimination steps
#ifdef LINBOX_pp_gauss_steps_OUT
// LINBOX_PRANK_OUT outputs intermediate ranks
#  ifndef LINBOX_PRANK_OUT
#    define LINBOX_PRANK_OUT
#  endif
#endif


namespace LinBox
{

        /** \brief Repository of functions for rank modulo 
         * a prime power by elimination on sparse matrices.
         * Specialization for powers of 2
         */
    template<typename UnsignedIntType>
    class PowerGaussDomainPowerOfTwo  {
        typedef UnsignedIntType UInt_t;
        typedef UnsignedIntType Element;
        const Element zero;
        const Element one;
    public:

            /** \brief The field parameter is the domain
             * over which to perform computations
             */
        PowerGaussDomainPowerOfTwo () : zero(0U), one(1U) {}

            //Copy constructor
            ///
        PowerGaussDomainPowerOfTwo (const PowerGaussDomainPowerOfTwo &M) {}



            // --------------------------------------------
            // Modulo operators
        bool isNZero(const UInt_t& a ) const { return (bool)a ;}
        bool isZero(const UInt_t& a ) const { return a == 0U;}
        bool isOne(const UInt_t& a ) const { return a == 1U;}
            /// @todo use Givaro isOdd
        bool isOdd(const UInt_t& b) const {
            return (bool)(b & 1U);
        }

        bool MY_divides(const UInt_t& a, const UInt_t& b) const {
            return (!(b%a));
        }

            // [On Newton-Raphson iteration for multiplicative
            //  inverses modulo prime powers. J-G. Dumas.
            //  IEEE Trans. on Computers, 63(8), pp 2106-2109, 2014]
            // http://doi.org/10.1109/TC.2013.94
        UInt_t& MY_Zpz_inv (UInt_t& u1, const UInt_t& a, const size_t exponent, const UInt_t& TWOTOEXPMONE) const {
            static const UInt_t ttep2(TWOTOEXPMONE+3U);
            if (this->isOne(a)) return u1=this->one;
            REQUIRE( (one<<exponent) == (TWOTOEXPMONE+1U) );
            REQUIRE( a <= TWOTOEXPMONE );
            REQUIRE( (a & 1U) );
            u1=ttep2-a; // 2-a
            UInt_t xmone(a-1);
            for(size_t i=2; i<exponent; i<<=1) {
                UInt_t tmp(xmone); xmone *= tmp;
                xmone &= TWOTOEXPMONE;
                u1 *= ++xmone; --xmone;
                u1 &= TWOTOEXPMONE;
            }
            ENSURE( ((a * u1) & TWOTOEXPMONE) == 1U );
            return u1;
        }

            // ------------------------------------------------
            // Pivot Searchers and column strategy
            // ------------------------------------------------
        template<class Vecteur, class D>
        void SameColumnPivoting(const Vecteur& lignepivot, size_t& indcol, long& indpermut, D& columns, Boolean_Trait<true>::BooleanType ) {
                // Try first in the same column
            size_t nj = (size_t) lignepivot.size() ;
            if (nj && (indcol == lignepivot[0].first) && (this->isOdd((UInt_t)lignepivot[0].second) ) ) {
                indpermut = (long)indcol;
                for(size_t j=nj;j--;)
                    --columns[ lignepivot[(size_t)j].first ];
                ++indcol;
            }
        }

        template<class BB, class Mmap, class D>
        bool SameColumnPivotingTrait(size_t& p, const BB& LigneA, const Mmap& psizes, size_t& indcol, long& indpermut, D& columns, Boolean_Trait<false>::BooleanType ) {
                // Do not try first in the same column
            return false;
        }

        template<class BB, class Mmap, class D>
        bool SameColumnPivotingTrait(size_t& p, const BB& LigneA, const Mmap& psizes, size_t& indcol, long& c, D& columns, Boolean_Trait<true>::BooleanType truetrait) {
            c=-2;
            for( typename Mmap::const_iterator iter = psizes.begin(); iter != psizes.end(); ++iter) {
                p = (size_t) (*iter).second;
                SameColumnPivoting(LigneA[(size_t)p], indcol, c, columns, truetrait ) ;
                if (c > -2 ) break;
            }
            if (c > -2)
                return true;
            else
                return false;

        }

        template<class Vecteur, class D>
        void CherchePivot(Vecteur& lignepivot, size_t& indcol , long& indpermut, D& columns ) {
            typedef typename Vecteur::value_type E;
            long nj =  (long) lignepivot.size() ;
            if (nj) {
                indpermut = (long)lignepivot[0].first;
                long pp=0;
                for(;pp<nj;++pp)
                    if (this->isOdd((UInt_t)lignepivot[(size_t)pp].second) ) break;

                if (pp < nj) {
                    long ds = (long)columns[ lignepivot[(size_t)pp].first ],p=pp,j=pp;
                    for(++j;j<nj;++j){
                        long dl;
                        if ( ( (dl=(long)columns[(size_t)lignepivot[(size_t)j].first] ) < ds ) && (this->isOdd((UInt_t)lignepivot[(size_t)j].second) ) ) {
                            ds = dl;
                            p = j;
                        }
                    }
                    if (p != 0) {
                        if (indpermut == (long)indcol) {
                            UInt_t ttm = (UInt_t)lignepivot[(size_t)p].second;
                            indpermut = (long)lignepivot[(size_t)p].first;
                            lignepivot[(size_t)p].second = (lignepivot[0].second);
                            lignepivot[0].second = (UInt_t)(ttm);
                        }
                        else {
                            E ttm = (E) lignepivot[(size_t)p];
                            indpermut = (long)ttm.first;
                            for(long m=p;m;--m)
                                lignepivot[(size_t)m] = lignepivot[(size_t)m-1];
                            lignepivot[0] = ttm;
                        }
                    }
                    for(j=nj;j--;)
                        --columns[ lignepivot[(size_t)j].first ];
                    if (indpermut != (long)indcol)
                        lignepivot[0].first = (indcol);
                    indcol++ ;
                }
                else
                    indpermut = -2;
            }
            else
                indpermut = -1;
        }

        template<class Vecteur>
        void PermuteColumn(Vecteur& lignecourante,
                           const size_t& nj,
                           const size_t& k,
                           const long& indpermut) {
            REQUIRE( nj > 0 );
            REQUIRE( indpermut > (long)k );

                // Find first non-zero element whose index is
                // greater than required permutation, if it exists
            size_t j_head(0);
            for(; j_head<nj; ++j_head)
                if (long(lignecourante[(size_t)j_head].first) >= indpermut) break;

            if ((j_head<nj) && (long(lignecourante[(size_t)j_head].first) == indpermut)) {
                size_t l(0);
                for(; l<j_head; ++l)
                    if (lignecourante[(size_t)l].first >= k) break;
                if (l<j_head && lignecourante[l].first == k) {
                        // non zero  <--> non zero
                    auto tmp = lignecourante[l].second ;
                    lignecourante[l].second = lignecourante[(size_t)j_head].second;
                    lignecourante[(size_t)j_head].second = tmp;
                } else {
                        // zero <--> non zero
                    auto tmp = lignecourante[(size_t)j_head];
                    tmp.first = (k);
                    for(size_t ll=j_head; ll>l; ll--)
                        lignecourante[ll] = lignecourante[ll-1];
                    lignecourante[l] = tmp;
                }
            } else {
                    // -------------------------------------------
                    // Permutation
                size_t l(0);
                for(; l<nj; ++l)
                    if (lignecourante[(size_t)l].first >= k) break;
                if ((l<nj) && (lignecourante[(size_t)l].first == k)) {
                        // non zero <--> zero
                    auto tmp = lignecourante[(size_t)l];
                    tmp.first = (size_t)(indpermut);
					const size_t bjh(j_head-1); // Position before j_head
                    for(;l<bjh;l++)
                        lignecourante[(size_t)l] = lignecourante[(size_t)l+1];
                    lignecourante[bjh] = tmp;
                } // else
                    // zero <--> zero
            }
        }

        template<class Vecteur>
        void PreserveUpperMatrixRow(Vecteur&, Boolean_Trait<true>::BooleanType ) {}

        template<class Vecteur>
        void PreserveUpperMatrixRow(Vecteur& ligne, Boolean_Trait<false>::BooleanType ) {
            ligne = Vecteur(0);
        }

        template<class SpMat>
        void PermuteSubMatrix(SpMat& LigneA,
							  const size_t & start,
                              const size_t & stop,
							  const size_t & currentrank,
                              const long & c) {
			for(size_t l=start; l < stop; ++l)
				if ( LigneA[(size_t)l].size() )
					PermuteColumn(LigneA[(size_t)l], LigneA[(size_t)l].size(), currentrank, c);
		}

        template<class SpMat>
        void PermuteUpperMatrix(SpMat& LigneA, const size_t & k, const size_t & currentrank, const long & c, Boolean_Trait<true>::BooleanType ) {
			PermuteSubMatrix(LigneA, 0, k, currentrank, c);
		}

        template<class SpMat>
        void PermuteUpperMatrix(SpMat&, const size_t &, const size_t &, const long &, Boolean_Trait<false>::BooleanType ) {}

        template<class Vecteur, class De>
        void FaireElimination( const size_t EXPONENT, const UInt_t& TWOK, const UInt_t& TWOKMONE,
                               Vecteur& lignecourante,
                               const Vecteur& lignepivot,
                               const UInt_t& invpiv,
                               const size_t& k,
                               const long& indpermut,
                               De& columns) {

            typedef typename Vecteur::value_type E;

            size_t nj = (size_t) lignecourante.size() ;

            if (nj) {
                if (lignecourante[0].first == k) {
                        // -------------------------------------------
                        // Head non-zero ==> Elimination
                    size_t npiv = (size_t) lignepivot.size();
                    Vecteur construit(nj + npiv);
                        // construit : <-- ci
                        // courante  : <-- m
                        // pivot     : <-- l
                    auto ci = construit.begin();
                    size_t m=1;
                    size_t l(0);

                        // A[(size_t)i,k] <-- A[(size_t)i,k] / A[(size_t)k,k]
                    UInt_t headcoeff = TWOK-(UInt_t)(lignecourante[0].second);
                    headcoeff *= invpiv;
                    headcoeff &= TWOKMONE ;

                    --columns[ lignecourante[0].first ];

                    for(;l<npiv;++l)
                        if (lignepivot[(size_t)l].first > k) break;
                        // for all j such that (j>k) and A[(size_t)k,j]!=0
                    for(;l<npiv;++l) {
                        size_t j_piv;
                        j_piv = (size_t) lignepivot[(size_t)l].first;
                            // if A[(size_t)k,j]=0,
                            // then A[(size_t)i,j] <-- A[(size_t)i,j]
                        for (;(m<nj) && (lignecourante[(size_t)m].first < j_piv);)
                            *ci++ = lignecourante[(size_t)m++];
                            // if A[(size_t)i,j]!=0, then A[(size_t)i,j]
                            // <-- A[(size_t)i,j] - A[(size_t)i,k]*A[(size_t)k,j]
                        if ((m<nj) && (lignecourante[(size_t)m].first == j_piv)) {
                            STATE( UInt_t lcs = lignecourante[(size_t)m].second);

                            lignecourante[(size_t)m].second += ( headcoeff  *  (UInt_t)lignepivot[(size_t)l].second );
                            lignecourante[(size_t)m].second &= TWOKMONE;

                            if (isNZero((UInt_t)(lignecourante[(size_t)m].second)))
                                *ci++ = lignecourante[(size_t)m++];
                            else {
                                --columns[ lignecourante[(size_t)m++].first ];
							}
                        } else {
                            UInt_t tmp(headcoeff);
                            tmp *= (UInt_t)lignepivot[(size_t)l].second;
                            tmp &= TWOKMONE;
                            if (isNZero(tmp)) {
                                ++columns[(size_t)j_piv];
                                *ci++ =  E(j_piv, (UInt_t)tmp );
                            }
                        }
                    }
                        // if A[(size_t)k,j]=0,
                        // then A[(size_t)i,j] <-- A[(size_t)i,j]
                    for (;m<nj;)
                        *ci++ = lignecourante[(size_t)m++];

                    construit.erase(ci,construit.end());
                    lignecourante = construit;
                }
            }
        }

            // ------------------------------------------------------
            // Rank calculators, defining row strategy
            // ------------------------------------------------------

        template<class BB, class D, class Container, class Perm, bool PrivilegiateNoColumnPivoting, bool PreserveUpperMatrix>
        void gauss_rankin(size_t EXPONENTMAX, Container& ranks, BB& LigneA, Perm& Q, const size_t Ni, const size_t Nj, const D& density_trait)
            {
                commentator().start ("Gaussian elimination with reordering modulo a prime power of 2",
                                     "PRGEPo2", Ni);

                linbox_check( Q.coldim() == Q.rowdim() );
                linbox_check( Q.coldim() == Nj );

                ranks.resize(0);

                typedef typename BB::Row Vecteur;
                uint64_t EXPONENT = EXPONENTMAX;
                UInt_t TWOK(1U); TWOK <<= EXPONENT;
                UInt_t TWOKMONE(TWOK); --TWOKMONE;
                ENSURE( TWOK == (UInt_t(1U) << EXPONENT) );
                ENSURE( TWOKMONE == (TWOK - 1U) );



#ifdef LINBOX_PRANK_OUT
                std::cerr << "Elimination mod " << TWOK << std::endl;
#endif

                D col_density(Nj);

                    // assignment of LigneA with the domain object
                    // and computation of the actual density
                for(size_t jj=0; jj<Ni; ++jj) {
                    Vecteur toto;
                    for(auto const & iter : LigneA[(size_t)jj]) {
                        UInt_t r = ((UInt_t)iter.second) & TWOKMONE;
                        if (isNZero(r)) {
                            ++col_density[ iter.first ];
                            toto.emplace_back(iter.first,r);
                        }
                    }
                    LigneA[(size_t)jj] = toto;
                }

                size_t last = Ni-1;
                long c(0);
                size_t indcol(0);
#ifdef LINBOX_PRANK_OUT
                size_t ind_pow = 1;
#endif
                size_t maxout = Ni/100; maxout = (maxout<10 ? 10 : (maxout>1000 ? 1000 : maxout) );
                size_t thres = Ni/maxout; thres = (thres >0 ? thres : 1);


                for (size_t k=0; k<last;++k) {
                    if ( ! (k % maxout) ) commentator().progress ((long)k);

                        // Look for invertible pivot
                    size_t p=k;
                    for(;;) {
                            // Order the rows 
                        std::multimap< long, long > psizes;
                        for(p=k; p<Ni; ++p)
                            psizes.insert( psizes.end(), std::pair<long,long>( (long)LigneA[(size_t)p].size(), (long)p) );

#ifdef  LINBOX_pp_gauss_steps_OUT
                        std::cerr << "------------ ordered rows " << k << " -----------" << std::endl;
                        for(auto const& iter: psizes)
                            std::cerr << iter.second << " : #" << iter.first << std::endl;
                        std::cerr << "---------------------------------------" << std::endl;
#endif

                        if ( SameColumnPivotingTrait(p, LigneA, psizes, indcol, c, col_density, typename Boolean_Trait<PrivilegiateNoColumnPivoting>::BooleanType() ) )
                            break;

                        for(auto const& iter: psizes) {
                            p = (size_t) iter.second;

                            CherchePivot( LigneA[(size_t)p], indcol, c , col_density) ;
                            if (c > -2 ) break;
                        }

                        if (c > -2) break;

                            // No invertible pivot found
                            // reduce everything by one power of 2
                        for(size_t ii=k;ii<Ni;++ii)
                            for(size_t jjj=LigneA[(size_t)ii].size();jjj--;)
                                LigneA[(size_t)ii][(size_t)jjj].second >>= 1;

                        --EXPONENT;
                        TWOK >>= 1;
                        TWOKMONE >>=1;

                        ENSURE( TWOK == (UInt_t(1U) << EXPONENT) );
                        ENSURE( TWOKMONE == (TWOK - 1U) );

                        ranks.push_back( indcol );
#ifdef LINBOX_PRANK_OUT
                        ++ind_pow;
                        std::cerr << "Rank mod 2^" << ind_pow << " : " << indcol << std::endl;
                        if (TWOK == 1) std::cerr << "wattadayada inhere ?" << std::endl;
#endif

                    }
                    if (p != k) {
#ifdef  LINBOX_pp_gauss_steps_OUT
                        std::cerr << "------------ permuting rows " << p << " and " << k << " ---" << std::endl;
#endif
                        Vecteur vtm = LigneA[(size_t)k];
                        LigneA[(size_t)k] = LigneA[(size_t)p];
                        LigneA[(size_t)p] = vtm;
                    }
                    if (c != -1) {
                            // Pivot has been found
                        REQUIRE( indcol > 0);
                        const size_t currentrank(indcol-1); 

                        if (c != (long)currentrank) {
#ifdef  LINBOX_pp_gauss_steps_OUT
                            std::cerr << "------------ permuting cols " << currentrank << " and " << c << " ---" << std::endl;
#endif
                            Q.permute(currentrank,c);
                            std::swap(col_density[currentrank],col_density[c]);
							PermuteUpperMatrix(LigneA, k, currentrank, c, typename Boolean_Trait<PreserveUpperMatrix>::BooleanType());
							PermuteSubMatrix(LigneA, k+1, Ni, currentrank, c);
                        }

                            // Compute the inverse of the found pivot
                        UInt_t invpiv;
                        MY_Zpz_inv(invpiv, (UInt_t) (LigneA[(size_t)k][0].second), EXPONENT, TWOKMONE);
                        for(size_t l=k + 1; (l < Ni) && (col_density[currentrank]); ++l)
                            FaireElimination(EXPONENT, TWOK, TWOKMONE, LigneA[(size_t)l], LigneA[(size_t)k], invpiv, currentrank, c, col_density);
                    }
                    
#ifdef  LINBOX_pp_gauss_steps_OUT
                    std::cerr << "step[" << k << "], pivot: " << c << std::endl;
#  ifdef  LINBOX_pp_gauss_intermediate_OUT
                    LigneA.write(std::cerr, Tag::FileFormat::Maple ) << std::endl;
#  endif
#endif

                    PreserveUpperMatrixRow(LigneA[(size_t)k], typename Boolean_Trait<PreserveUpperMatrix>::BooleanType());
                }

                c = -2;
                SameColumnPivoting(LigneA[(size_t)last], indcol, c, col_density, typename Boolean_Trait<PrivilegiateNoColumnPivoting>::BooleanType() );
                if (c == -2) CherchePivot( LigneA[(size_t)last], indcol, c, col_density );
                while( c == -2) {
                    ranks.push_back( indcol );
                    for(long jjj=(long)LigneA[(size_t)last].size();jjj--;)
                        LigneA[(size_t)last][(size_t)jjj].second >>= 1;
                    TWOK >>= 1;
                    CherchePivot( LigneA[(size_t)last], indcol, c, col_density );
                }
                if (c != -1) {
                    const size_t currentrank(indcol-1);
                    if (c != (long)currentrank) {
#ifdef  LINBOX_pp_gauss_steps_OUT
                        std::cerr << "------------ permuting cols " << (indcol-1) << " and " << c << " ---" << std::endl;
#endif
                        Q.permute(currentrank,c);
						PermuteUpperMatrix(LigneA, last, currentrank, c, typename Boolean_Trait<PreserveUpperMatrix>::BooleanType());
                    }
                }
                while( TWOK > 1) {
                    TWOK >>= 1;
                    ranks.push_back( indcol );
                }

#ifdef LINBOX_pp_gauss_steps_OUT
                std::cerr << "step[" << Ni-1 << "], pivot: " << c << std::endl;
#  ifdef LINBOX_pp_gauss_intermediate_OUT
                LigneA.write(std::cerr, Tag::FileFormat::Maple ) << std::endl;
#  endif
#endif
#ifdef LINBOX_PRANK_OUT
                std::cerr << "Rank mod 2^" << EXPONENTMAX << " : " << indcol << std::endl;
#endif
                commentator().stop ("done", 0, "PRGEPo2");

            }

        template<class BB, class D, class Container, class Perm>
        void prime_power_rankin (size_t EXPONENT, Container& ranks, BB& SLA, Perm& Q, const size_t Ni, const size_t Nj, const D& density_trait, int StaticParameters=PRIVILEGIATE_NO_COLUMN_PIVOTING) {
            if (PRIVILEGIATE_NO_COLUMN_PIVOTING & StaticParameters) {
                if (PRESERVE_UPPER_MATRIX & StaticParameters) {
                    gauss_rankin<BB,D,Container,Perm,true,true>(EXPONENT,ranks, SLA, Q, Ni, Nj, density_trait);
                } else {
                    gauss_rankin<BB,D,Container,Perm,true,false>(EXPONENT,ranks, SLA, Q, Ni, Nj, density_trait);
                }
            } else {
                if (PRESERVE_UPPER_MATRIX & StaticParameters) {
                    gauss_rankin<BB,D,Container,Perm,false,true>(EXPONENT,ranks, SLA, Q, Ni, Nj, density_trait);
                } else {
                    gauss_rankin<BB,D,Container,Perm,false,false>(EXPONENT,ranks, SLA, Q, Ni, Nj, density_trait);
                }
            }
        }


        template<class Matrix, class Perm, template<class, class> class Container, template<class> class Alloc>
        Container<std::pair<UInt_t,size_t>, Alloc<std::pair<UInt_t,size_t> > >& operator()(Container<std::pair<UInt_t,size_t>, Alloc<std::pair<UInt_t,size_t> > >& L, Matrix& A, Perm& Q, size_t EXPONENT, int StaticParameters=PRIVILEGIATE_NO_COLUMN_PIVOTING) {
            Container<size_t, Alloc<size_t> > ranks;
            prime_power_rankin( EXPONENT, ranks, A, Q, A.rowdim(), A.coldim(), std::vector<size_t>(),StaticParameters);
            L.resize( 0 ) ;
            UInt_t MOD(1);
            size_t num = 0;
            for( typename Container<size_t, Alloc<size_t> >::const_iterator it = ranks.begin(); it != ranks.end(); ++it) {
                size_t diff = *it-num;
                if (diff > 0) L.emplace_back(MOD,diff);
                MOD <<= 1;
                num = *it;
            }
            return L;
        }

    };



} // end of LinBox namespace

#endif  //__LINBOX_pp_gauss_poweroftwo_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
