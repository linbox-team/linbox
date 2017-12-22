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
#include "linbox/algorithms/gauss.h"

#ifdef DEBUG
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
		Modulo& MY_Zpz_inv (Modulo& u1, const Modulo2 a, const Modulo3 _p) const
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
		template<class Modulo, class Modulo2>
		Modulo MY_Zpz_inv (const Modulo a, const Modulo2 _p) const
            {
                Modulo u1; return MY_Zpz_inv(u1,a,_p);
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
        void SameColumnPivoting(Modulo PRIME,  const Vecteur& lignepivot, unsigned long& indcol, long& indpermut, D& columns, Boolean_Trait<false>::BooleanType ) {}


        template<class Modulo, class Vecteur, class D>
        void SameColumnPivoting(Modulo PRIME,  const Vecteur& lignepivot, unsigned long& indcol, long& indpermut, D& columns, Boolean_Trait<true>::BooleanType ) {
                // Try first in the same column
			unsigned long nj = (unsigned long)  lignepivot.size() ;
			if (nj && (indcol == lignepivot[0].first) && (! this->MY_divides(PRIME,lignepivot[0].second) ) ) {
                indpermut = (long) indcol;
                for(unsigned long j=nj;j--;)
                    --columns[ lignepivot[(size_t)j].first ];
                ++indcol;
            }
        }

        template<class Modulo, class BB, class Mmap, class D>
        bool SameColumnPivotingTrait(Modulo PRIME, unsigned long& p, const BB& LigneA, const Mmap& psizes, unsigned long& indcol, long& indpermut, D& columns, Boolean_Trait<false>::BooleanType ) {
                // Do not try first in the same column
            return false;
        }

        template<class Modulo, class BB, class Mmap, class D>
        bool SameColumnPivotingTrait(Modulo PRIME, unsigned long& p, const BB& LigneA, const Mmap& psizes, unsigned long& indcol, long& c, D& columns, Boolean_Trait<true>::BooleanType truetrait) {
            c=-2;
            for( typename Mmap::const_iterator iter = psizes.begin(); iter != psizes.end(); ++iter) {
                p = (unsigned long)(*iter).second;
                SameColumnPivoting(PRIME, LigneA[(size_t)p], indcol, c, columns, truetrait ) ;
                if (c > -2 ) break;
            }
            if (c > -2)
                return true;
            else
                return false;

        }

		template<class Modulo, class Vecteur, class D>
		void CherchePivot(Modulo PRIME, Vecteur& lignepivot, unsigned long& indcol , long& indpermut, D& columns ) {
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
                           const unsigned long& nj,
                           const unsigned long& k,
                           const long& indpermut) {

            REQUIRE( nj > 0 );
            REQUIRE( indpermut > (long)k );

                // Find first non-zero element whose index is
                // greater than required permutation, if it exists
            unsigned long j_head(0);
            for(; j_head<nj; ++j_head)
                if (long(lignecourante[(size_t)j_head].first) >= indpermut) break;
            if ((j_head<nj) && (long(lignecourante[(size_t)j_head].first) == indpermut)) {
                unsigned long l(0);
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
                    for(unsigned long ll=j_head; ll>l; ll--)
                        lignecourante[ll] = lignecourante[ll-1];
                    lignecourante[l] = tmp;
                }                
            } else {
                    // -------------------------------------------
                    // Permutation
                unsigned long l(0);
                for(; l<nj; ++l)
                    if (lignecourante[(size_t)l].first >= k) break;
                if ((l<nj) && (lignecourante[(size_t)l].first == k))  {
                        // non zero <--> zero
                    auto tmp = lignecourante[(size_t)l];
                    tmp.first = (unsigned long) (indpermut);
                    const unsigned long bjh(j_head-1);  // Position before j_head
                    for(;l<bjh;l++)
                        lignecourante[(size_t)l] = lignecourante[(size_t)l+1];
                    lignecourante[bjh] = tmp;
                } // else
                    // zero <--> zero
            }
        }
        

        template<class SpMat>
        void PermuteSubMatrix(SpMat& LigneA,
							  const unsigned long & start, const unsigned long & stop,
							  const unsigned long & currentrank, const long & c) {
			for(unsigned long l=start; l < stop; ++l)
				if ( LigneA[(size_t)l].size() )
					PermuteColumn(LigneA[(size_t)l], LigneA[(size_t)l].size(), currentrank, c);
		}

        template<class SpMat>
        void PermuteUpperMatrix(SpMat& LigneA, const unsigned long & k, const unsigned long & currentrank, const long & c, Boolean_Trait<true>::BooleanType ) {
			PermuteSubMatrix(LigneA, 0, k, currentrank, c);
		}

        template<class SpMat>
        void PermuteUpperMatrix(SpMat&, const unsigned long &, const unsigned long &, const long &, Boolean_Trait<false>::BooleanType ) {}


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
                               const unsigned long& k,
                               const long& indpermut,
                               De& columns) {

                //     typedef typename Vecteur::coefficientSpace F;
                //     typedef typename Vecteur::value_types E;
			typedef typename Field::Element F;
			typedef typename Vecteur::value_type E;

			typedef typename Signed_Trait<Modulo>::unsigned_type UModulo;

			unsigned long nj = (unsigned long)  lignecourante.size() ;
			if (nj) {
                if (lignecourante[0].first == k) {
                        // -------------------------------------------
                        // Elimination
					unsigned long npiv = (unsigned long) lignepivot.size();
					Vecteur construit(nj + npiv);
                        // construit : <-- ci
                        // courante  : <-- m
                        // pivot     : <-- l
					typedef typename Vecteur::iterator Viter;
					Viter ci = construit.begin();
					unsigned long m=1;
					unsigned long l(0);
 
					F headcoeff = MOD-(lignecourante[0].second);
					headcoeff *= invpiv;
					headcoeff %= (UModulo)MOD ;
					lignecourante[0].second = headcoeff;
					--columns[ lignecourante[0].first ];

					for(;l<npiv;++l)
						if (lignepivot[(size_t)l].first > k) break;
                        // for all j such that (j>k) and A[(size_t)k,j]!=0
					for(;l<npiv;++l) {
						unsigned long j_piv;
						j_piv = (unsigned long) lignepivot[(size_t)l].first;
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
#ifdef LINBOX_PRANK_OUT
                std::cerr << "Elimination mod " << MOD << " (" << PrivilegiateNoColumnPivoting << ',' << PreserveUpperMatrix << ')' << std::endl;
#endif

                D col_density(Nj);

                    // assignment of LigneA with the domain object
                size_t jj;
                for(jj=0; jj<Ni; ++jj) {
                    Vecteur tmp = LigneA[(size_t)jj];
                    Vecteur toto(tmp.size());
                    unsigned long k=0,rs=0;
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

                unsigned long last = Ni-1;
                long c(0);
                unsigned long indcol(0);
                unsigned long ind_pow = 1;
                unsigned long maxout = Ni/100; maxout = (maxout<10 ? 10 : (maxout>1000 ? 1000 : maxout) );
                unsigned long thres = Ni/maxout; thres = (thres >0 ? thres : 1);


                for (unsigned long k=0; k<last;++k) {
                    if ( ! (k % maxout) ) commentator().progress ((long)k);


                    unsigned long p=k;
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
                            p = (unsigned long)(*iter).second;

                            CherchePivot( PRIME, LigneA[(size_t)p], indcol, c , col_density) ;
                            if (c > -2 ) break;
                        }

                        if (c > -2) break;
                        for(unsigned long ii=k;ii<Ni;++ii)
                            for(unsigned long jjj=LigneA[(size_t)ii].size();jjj--;)
                                LigneA[(size_t)ii][(size_t)jjj].second /= PRIME;
                        MOD /= PRIME;
                        ranks.push_back( indcol );
                        ++ind_pow;
#ifdef LINBOX_PRANK_OUT
                        std::cerr << "Rank mod " << PRIME << "^" << ind_pow << " : " << indcol << std::endl;
                        if (MOD == 1) std::cerr << "wattadayada inhere ?" << std::endl;
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
                        REQUIRE( indcol > 0);
                        const unsigned long currentrank(indcol-1); 
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
                        MY_Zpz_inv(invpiv, LigneA[(size_t)k][0].second, MOD);

                        for(unsigned long l=k + 1; l < Ni; ++l)
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
                    const unsigned long currentrank(indcol-1);
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

                    //             ranks.push_back(indcol);
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
		Container<std::pair<size_t,Modulo>, Alloc<std::pair<size_t,Modulo> > >& operator()(Container<std::pair<size_t,Modulo>, Alloc<std::pair<size_t,Modulo> > >& L, Matrix& A, Perm& Q, Modulo FMOD, Modulo PRIME, int StaticParameters=PRIVILEGIATE_NO_COLUMN_PIVOTING)
            {
                Container<size_t, Alloc<size_t> > ranks;
                prime_power_rankin( FMOD, PRIME, ranks, A, Q, A.rowdim(), A.coldim(), std::vector<size_t>(), StaticParameters);
                L.resize( 0 ) ;
                Modulo MOD = 1;
                size_t num = 0;
                for( typename Container<size_t, Alloc<size_t> >::const_iterator it = ranks.begin(); it != ranks.end(); ++it) {
                    size_t diff(*it-num);
                    if (diff > 0)
                        L.push_back( std::pair<size_t,Modulo>(diff,MOD) );
                    MOD *= PRIME;
                    num = *it;
                }
                return L;
            }

	};


} // end of LinBox namespace

#endif  //__LINBOX_pp_gauss_H

// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:

