/* Copyright (C) Givaro Team 1999
 * Copyright (C) LinBox
 * Written by JG Dumas
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_pp_gauss_H
#define __LINBOX_pp_gauss_H

#include <map>
#include <givaro/givconfig.h> // for Signed_Trait
#include "linbox/solutions/smith-form.h"
#include "linbox/algorithms/gauss.h"

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

    template <bool Boolean> struct Boolean_Trait;
    template <> struct Boolean_Trait<true> {
        typedef int BooleanType; // int does not matter, only that it differs from float
    };
    template <> struct Boolean_Trait<false> {
        typedef float BooleanType;// float does not matter, only that it differs from int
    };

    enum {
            // Combine these in binary for use in StaticParameters
        PRIVILEGIATE_NO_COLUMN_PIVOTING	= 1,
        PRIVILEGIATE_REDUCING_FILLIN	= 2,
        PRESERVE_UPPER_MATRIX		= 4
    };

        /** \brief Repository of functions for rank modulo 
         * a prime power by elimination on sparse matrices.
         */
	template <class _Field>
	class PowerGaussDomain : public GaussDomain<_Field> {
		typedef GaussDomain<_Field> Father_t;
		typedef _Field Field;
		typedef typename Field::Element Element;
	public:

            /** \brief The field parameter is the domain
             * over which to perform computations
             */
		PowerGaussDomain (const Field &F) :
                Father_t(F)
            {}

            //Copy constructor
            ///
		PowerGaussDomain (const PowerGaussDomain &M) :
                Father_t(M)
            {}



            // --------------------------------------------
            // Modulo operators
		template<class Modulu>
		bool isNZero(const Modulu& a ) const { return (bool)a ;}
		template<class Modulu>
		bool isZero(const Modulu& a ) const { return a == 0U;}

		template<class Modulo, class Modulo2, class Modulo3>
		Modulo& MY_Zpz_inv_classic (Modulo& u1, const Modulo2 a, const Modulo3 _p) const
            {
                u1 = Modulo(1U);
                Modulo r0((Modulo)_p), r1((Modulo)a); //! clang complains for examples/smith.C and examples/smithvalence.C
                Modulo q(r0/r1);

                r0 -= q * r1;
                if ( isZero(r0) ) return u1;
                Modulo u0 = q;

                q = r1/r0;
                r1 -= q * r0;

                while ( isNZero(r1) ) {
                    u1 += q * u0;

                    q = r0/r1;
                    r0 -= q * r1;
                    if ( isZero(r0) ) return u1;
                    u0 += q * u1;

                    q = r1/r0;
                    r1 -= q * r0;

                }

                return u1=(Modulo)_p-u0;
            }

            // [On Newton-Raphson iteration for multiplicative
            //  inverses modulo prime powers. J-G. Dumas.
            //  IEEE Trans. on Computers, 63(8), pp 2106-2109, 2014]
            // http://doi.org/10.1109/TC.2013.94
		template<class Modulo, class Modulo2, class Modulo3>
		Modulo& MY_Zpz_inv (Modulo& v1, const Modulo2& a, const Modulo3& _p, const Modulo3& _q, uint64_t exponent) const
            {
                Modulo2 ta(a%_p);
                MY_Zpz_inv_classic(v1, ta, _p);					// b=1/a % p
                Modulo3 ttep2(_q); ttep2 += 2U;					// q+2
                Modulo2 atimeb(a); atimeb*= v1; atimeb %= _q;	// ab
                ttep2 -= atimeb;								// 2-ab
                v1 *= ttep2; v1 %= _q;							// b(2-ab)
                Modulo abmone(atimeb); --abmone; 				// ab-1
                for(size_t i=2; i<exponent; i<<=1) {
                    Modulo tmp(abmone); abmone *= tmp;			// (ab-1)^(2^i)
                    abmone %= _q;
                    v1 *= ++abmone; --abmone;					// in place
                    v1 %= _q;
                }
                return v1;
            }

		template<class Ring1, class Ring2>
		bool MY_divides(Ring1 a, Ring2 b) const
            {
                return (!(b%a));
            }

            // ------------------------------------------------
            // Pivot Searchers and column strategy
            // ------------------------------------------------
        template<class Modulo, class Vecteur, class D>
        void SameColumnPivoting(Modulo PRIME,  const Vecteur& lignepivot, size_t& indcol, long& indpermut, D& columns, Boolean_Trait<false>::BooleanType ) {}


        template<class Modulo, class Vecteur, class D>
        void SameColumnPivoting(Modulo PRIME,  const Vecteur& lignepivot, size_t& indcol, long& indpermut, D& columns, Boolean_Trait<true>::BooleanType ) {
                // Try first in the same column
			size_t nj = (size_t)  lignepivot.size() ;
			if (nj && (indcol == lignepivot[0].first) && (! this->MY_divides(PRIME,lignepivot[0].second) ) ) {
                indpermut = (long) indcol;
                for(size_t j=nj;j--;)
                    --columns[ lignepivot[(size_t)j].first ];
                ++indcol;
            }
        }

        template<class Modulo, class BB, class Mmap, class D>
        bool SameColumnPivotingTrait(Modulo PRIME, size_t& p, const BB& LigneA, const Mmap& psizes, size_t& indcol, long& indpermut, D& columns, Boolean_Trait<false>::BooleanType ) {
                // Do not try first in the same column
            return false;
        }

        template<class Modulo, class BB, class Mmap, class D>
        bool SameColumnPivotingTrait(Modulo PRIME, size_t& p, const BB& LigneA, const Mmap& psizes, size_t& indcol, long& c, D& columns, Boolean_Trait<true>::BooleanType truetrait) {
            c=-2;
            for( typename Mmap::const_iterator iter = psizes.begin(); iter != psizes.end(); ++iter) {
                p = (size_t)(*iter).second;
                SameColumnPivoting(PRIME, LigneA[(size_t)p], indcol, c, columns, truetrait ) ;
                if (c > -2 ) break;
            }
            if (c > -2)
                return true;
            else
                return false;

        }

		template<class Modulo, class Vecteur, class D>
		void CherchePivot(Modulo PRIME, Vecteur& lignepivot, size_t& indcol , long& indpermut, D& columns ) {
            long nj =(long)  lignepivot.size() ;
            if (nj) {
                indpermut = (long)lignepivot[0].first;
                long pp=0;
                for(;pp<nj;++pp)
                    if (! this->MY_divides(PRIME,lignepivot[(size_t)pp].second) ) break;

                if (pp < nj) {

                    long ds = (long)columns[ lignepivot[(size_t)pp].first ],p=pp,j=pp;
                    for(++j;j<nj;++j){
                        long dl;
                        if ( ( (dl=(long)columns[(size_t)lignepivot[(size_t)j].first] ) < ds ) && (! MY_divides(PRIME,lignepivot[(size_t)j].second) ) ) {
                            ds = dl;
                            p = j;
                        }
                    }
                    if (p != 0) {
                        if (indpermut == (long)indcol) {
                            auto ttm = lignepivot[(size_t)p].second;
                            indpermut = (long)lignepivot[(size_t)p].first;
                            lignepivot[(size_t)p].second = (lignepivot[0].second);
                            lignepivot[0].second = (ttm);
                        }
                        else {
                            auto ttm = lignepivot[(size_t)p];
                            indpermut = (long)ttm.first;
                            for(long m=p;m;--m)
                                lignepivot[(size_t)m] = lignepivot[(size_t)m-1];
                            lignepivot[0] = ttm;
                        }
                    }
                    for(j=nj;j--;)
                        --columns[ lignepivot[(size_t)j].first ];
                    if (indpermut != (long)indcol) {
#ifdef  LINBOX_pp_gauss_steps_OUT
                        std::cerr << "------CP---- permuting cols " << indcol << " and " << indpermut << " ---" << std::endl;
#endif

                        lignepivot[0].first = (indcol);
                    }
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
                    // -------------------------------------------
                    // Permutation
                if (l<j_head && lignecourante[l].first == k) {
                        // non zero  <--> non zero
                    auto tmp = lignecourante[l].second ;
                    lignecourante[l].second = (lignecourante[(size_t)j_head].second );
                    lignecourante[(size_t)j_head].second = (tmp);
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
                if ((l<nj) && (lignecourante[(size_t)l].first == k))  {
                        // non zero <--> zero
                    auto tmp = lignecourante[(size_t)l];
                    tmp.first = (size_t) (indpermut);
                    const size_t bjh(j_head-1);  // Position before j_head
                    for(;l<bjh;l++)
                        lignecourante[(size_t)l] = lignecourante[(size_t)l+1];
                    lignecourante[bjh] = tmp;
                } // else
                    // zero <--> zero
            }
        }
        

        template<class SpMat>
        void PermuteSubMatrix(SpMat& LigneA,
							  const size_t & start, const size_t & stop,
							  const size_t & currentrank, const long & c) {
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


        template<class Vecteur>
        void PreserveUpperMatrixRow(Vecteur& ligne, Boolean_Trait<true>::BooleanType ) {}

        template<class Vecteur>
        void PreserveUpperMatrixRow(Vecteur& ligne, Boolean_Trait<false>::BooleanType ) {
            ligne = Vecteur(0);
        }


		template<class Modulo, class Vecteur, class De>
		void FaireElimination( Modulo MOD,
                               Vecteur& lignecourante,
                               const Vecteur& lignepivot,
                               const typename Signed_Trait<Modulo>::unsigned_type& invpiv,
                               const size_t& k,
                               const long& indpermut,
                               De& columns) {

                //     typedef typename Vecteur::coefficientSpace F;
                //     typedef typename Vecteur::value_types E;
			typedef typename Field::Element F;
			typedef typename Vecteur::value_type E;

			typedef typename Signed_Trait<Modulo>::unsigned_type UModulo;

			size_t nj = (size_t)  lignecourante.size() ;
			if (nj) {
                if (lignecourante[0].first == k) {
                        // -------------------------------------------
                        // Elimination
					size_t npiv = (size_t) lignepivot.size();
					Vecteur construit(nj + npiv);
                        // construit : <-- ci
                        // courante  : <-- m
                        // pivot     : <-- l
					typedef typename Vecteur::iterator Viter;
					Viter ci = construit.begin();
					size_t m=1;
					size_t l(0);
 
					F headcoeff = MOD-(lignecourante[0].second);
					headcoeff *= invpiv;
					headcoeff %= (UModulo)MOD ;
					lignecourante[0].second = headcoeff;
					--columns[ lignecourante[0].first ];

					for(;l<npiv;++l)
						if (lignepivot[(size_t)l].first > k) break;
                        // for all j such that (j>k) and A[(size_t)k,j]!=0
					for(;l<npiv;++l) {
						size_t j_piv;
						j_piv = (size_t) lignepivot[(size_t)l].first;
                            // if A[(size_t)k,j]=0, then A[(size_t)i,j] <-- A[(size_t)i,j]
						for (;(m<nj) && (lignecourante[(size_t)m].first < j_piv);)
							*ci++ = lignecourante[(size_t)m++];
                            // if A[(size_t)i,j]!=0, then A[(size_t)i,j] <-- A[(size_t)i,j] - A[(size_t)i,k]*A[(size_t)k,j]
						if ((m<nj) && (lignecourante[(size_t)m].first == j_piv)) {
                                //lignecourante[(size_t)m].second = ( ((UModulo)( headcoeff  *  lignepivot[(size_t)l].second  + lignecourante[(size_t)m].second ) ) % (UModulo)MOD );
							lignecourante[(size_t)m].second += ( headcoeff  *  lignepivot[(size_t)l].second );
                            lignecourante[(size_t)m].second %= (UModulo)MOD;
							if (isNZero(lignecourante[(size_t)m].second))
								*ci++ = lignecourante[(size_t)m++];
							else
								--columns[ lignecourante[(size_t)m++].first ];
                                //                         m++;
						}
						else {
							F tmp(headcoeff);
							tmp *= lignepivot[(size_t)l].second;
							tmp %= (UModulo)MOD;
							if (isNZero(tmp)) {
								++columns[(size_t)j_piv];
								*ci++ =  E(j_piv, tmp );
							}
						}
					}
                        // if A[(size_t)k,j]=0, then A[(size_t)i,j] <-- A[(size_t)i,j]
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

		template<class Modulo, class BB, class D, class Container, class Perm, bool PrivilegiateNoColumnPivoting, bool PreserveUpperMatrix>
		void gauss_rankin(Modulo FMOD, Modulo PRIME, Container& ranks, BB& LigneA, Perm& Q, const size_t Ni, const size_t Nj, const D& density_trait)
            {
                linbox_check( Q.coldim() == Q.rowdim() );
                linbox_check( Q.coldim() == Nj );

                commentator().start ("Gaussian elimination with reordering modulo a prime power",
                                     "PRGE", Ni);

                ranks.resize(0);

                typedef typename BB::Row Vecteur;
                typedef typename Signed_Trait<Modulo>::unsigned_type UModulo;

                Modulo MOD = FMOD;
                uint64_t exponent(1);
                Modulo tq(FMOD);
                while(tq > PRIME) {
                    tq /= PRIME;
                    ++exponent;
                }
#ifdef LINBOX_PRANK_OUT
                std::cerr << "Elimination mod " << MOD << " (=" << PRIME << '^' << exponent << ',' << PrivilegiateNoColumnPivoting << ',' << PreserveUpperMatrix << ')' << std::endl;
#endif

                D col_density(Nj);

                    // assignment of LigneA with the domain object
                size_t jj;
                for(jj=0; jj<Ni; ++jj) {
                    Vecteur tmp = LigneA[(size_t)jj];
                    Vecteur toto(tmp.size());
                    size_t k=0,rs=0;
                    for(; k<tmp.size(); ++k) {
                        Modulo r = tmp[(size_t)k].second;
                        if ((r <0) || (r >= MOD)) r %= MOD ;
                        if (r <0) r += MOD ;
                        if (isNZero(r)) {
                            ++col_density[ tmp[(size_t)k].first ];
                            toto[rs] =tmp[(size_t)k];
                            toto[rs].second = ( r );
                            ++rs;
                        }
                    }
                    toto.resize(rs);
                    LigneA[(size_t)jj] = toto;
                        //                 LigneA[(size_t)jj].reactualsize(Nj);

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


                    size_t p=k;
                    for(;;) {


                        std::multimap< long, long > psizes;
                        for(p=k; p<Ni; ++p)
                            psizes.insert( psizes.end(), std::pair<long,long>( (long) LigneA[(size_t)p].size(), (long) p) );

#ifdef  LINBOX_pp_gauss_steps_OUT
                        std::cerr << "------------ ordered rows " << k << " -----------" << std::endl;
                        for(auto const& iter: psizes)
                            std::cerr << iter.second << " : #" << iter.first << std::endl;
                        std::cerr << "---------------------------------------" << std::endl;
#endif

                        if ( SameColumnPivotingTrait(PRIME, p, LigneA, psizes, indcol, c, col_density, typename Boolean_Trait<PrivilegiateNoColumnPivoting>::BooleanType() ) )
                            break;

                        for( typename std::multimap< long, long >::const_iterator iter = psizes.begin(); iter != psizes.end(); ++iter) {
                            p = (size_t)(*iter).second;

                            CherchePivot( PRIME, LigneA[(size_t)p], indcol, c , col_density) ;
                            if (c > -2 ) break;
                        }

                        if (c > -2) break;
                        for(size_t ii=k;ii<Ni;++ii)
                            for(size_t jjj=LigneA[(size_t)ii].size();jjj--;)
                                LigneA[(size_t)ii][(size_t)jjj].second /= PRIME;
                        MOD /= PRIME;
                        --exponent;

#ifdef LINBOX_PRANK_OUT
                        std::cerr << "Rank mod " << PRIME << "^" << ind_pow << " : " << indcol << std::endl;
                        if (MOD == 1) std::cerr << "wattadayada inhere ?" << std::endl;
                        ++ind_pow;
#endif
                        ranks.push_back( indcol );

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
                        REQUIRE( indcol > 0);
                        const size_t currentrank(indcol-1); 
                        if (c != (long)currentrank) {
#ifdef  LINBOX_pp_gauss_steps_OUT
                            std::cerr << "------------ permuting cols " << (indcol-1) << " and " << c << " ---" << std::endl;
#endif
                            Q.permute(indcol-1,c);
                            std::swap(col_density[currentrank],col_density[c]);
							PermuteUpperMatrix(LigneA, k, currentrank, c, typename Boolean_Trait<PreserveUpperMatrix>::BooleanType());
							PermuteSubMatrix(LigneA, k+1, Ni, currentrank, c);
                        }

                        UModulo invpiv; 
                        MY_Zpz_inv(invpiv, LigneA[(size_t)k][0].second, PRIME, MOD, exponent);

                        for(size_t l=k + 1; (l < Ni) && (col_density[currentrank]); ++l)
                            FaireElimination(MOD, LigneA[(size_t)l], LigneA[(size_t)k], invpiv, currentrank, c, col_density);
                    }
                

#ifdef  LINBOX_pp_gauss_steps_OUT
                    LigneA.write(std::cerr << "step[" << k << "], pivot: " << c << std::endl, Tag::FileFormat::Maple) << std::endl;
#endif

                    PreserveUpperMatrixRow(LigneA[(size_t)k], typename Boolean_Trait<PreserveUpperMatrix>::BooleanType());
                }

                c = -2;
                SameColumnPivoting(PRIME, LigneA[(size_t)last], indcol, c, col_density, typename Boolean_Trait<PrivilegiateNoColumnPivoting>::BooleanType() );
                if (c == -2) CherchePivot( PRIME, LigneA[(size_t)last], indcol, c, col_density );
                while( c == -2) {
                    ranks.push_back( indcol );
                    for(long jjj=(long)LigneA[(size_t)last].size();jjj--;)
                        LigneA[(size_t)last][(size_t)jjj].second /= PRIME;
                    MOD /= PRIME;
                    CherchePivot( PRIME, LigneA[(size_t)last], indcol, c, col_density );
                }
                if (c != -1) {
                    const size_t currentrank(indcol-1);
                    if (c != (long)currentrank) {
#ifdef  LINBOX_pp_gauss_steps_OUT
                        std::cerr << "------------ permuting cols " << (indcol-1) << " and " << c << " ---" << std::endl;
#endif
                        Q.permute(indcol-1,c);
						PermuteUpperMatrix(LigneA, last, currentrank, c, typename Boolean_Trait<PreserveUpperMatrix>::BooleanType());
                    }
                }
                while( MOD > 1) {
                    MOD /= PRIME;
                    ranks.push_back( indcol );
                }

#ifdef LINBOX_pp_gauss_steps_OUT
                LigneA.write(std::cerr << "step[" << Ni-1 << "], pivot: " << c << std::endl, Tag::FileFormat::Maple) << std::endl;
#endif
#ifdef LINBOX_PRANK_OUT
                std::cerr << "Rank mod " << FMOD << " : " << indcol << std::endl;
#endif
                commentator().stop ("done", 0, "PRGE");

            }

		template<class Modulo, class BB, class D, class Container, class Perm>
		void prime_power_rankin (Modulo FMOD, Modulo PRIME, Container& ranks, BB& SLA, Perm& Q, const size_t Ni, const size_t Nj, const D& density_trait, int StaticParameters=PRIVILEGIATE_NO_COLUMN_PIVOTING)
            {
                if (PRIVILEGIATE_NO_COLUMN_PIVOTING & StaticParameters) {
                    if (PRESERVE_UPPER_MATRIX & StaticParameters) {
                        gauss_rankin<Modulo,BB,D,Container,Perm,true,true>(FMOD,PRIME,ranks, SLA, Q, Ni, Nj, density_trait);
                    } else {
                        gauss_rankin<Modulo,BB,D,Container,Perm,true,false>(FMOD,PRIME,ranks, SLA, Q, Ni, Nj, density_trait);
                    }
                } else {
                    if (PRESERVE_UPPER_MATRIX & StaticParameters) {
                        gauss_rankin<Modulo,BB,D,Container,Perm,false,true>(FMOD,PRIME,ranks, SLA, Q, Ni, Nj, density_trait);
                    } else {
                        gauss_rankin<Modulo,BB,D,Container,Perm,false,false>(FMOD,PRIME,ranks, SLA, Q, Ni, Nj, density_trait);
                    }
                }
            }


		template<class Modulo, class Matrix, class Perm, template<class, class> class Container, template<class> class Alloc>
		Container<std::pair<Modulo,size_t>, Alloc<std::pair<Modulo,size_t> > >& operator()(Container<std::pair<Modulo,size_t>, Alloc<std::pair<Modulo,size_t> > >& L, Matrix& A, Perm& Q, Modulo FMOD, Modulo PRIME, int StaticParameters=PRIVILEGIATE_NO_COLUMN_PIVOTING)
            {
                Container<size_t, Alloc<size_t> > ranks;
                prime_power_rankin( FMOD, PRIME, ranks, A, Q, A.rowdim(), A.coldim(), std::vector<size_t>(), StaticParameters);
                L.resize( 0 ) ;
                Modulo MOD = 1;
                size_t num = 0;
                for( typename Container<size_t, Alloc<size_t> >::const_iterator it = ranks.begin(); it != ranks.end(); ++it) {
                    size_t diff(*it-num);
                    if (diff > 0)
                        L.emplace_back(MOD,diff);
                    MOD *= PRIME;
                    num = *it;
                }
                return L;
            }

	};


} // end of LinBox namespace

#endif  //__LINBOX_pp_gauss_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
