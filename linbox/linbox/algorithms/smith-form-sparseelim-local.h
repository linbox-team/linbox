/* Copyright (C) Givaro Team 1999
 * Copyright (C) LinBox
 * Written by JG Dumas
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


#ifndef __LINBOX_pp_gauss_H
#define __LINBOX_pp_gauss_H

#include <map>
#include <givaro/givconfig.h> // for Signed_Trait
#include <linbox/algorithms/gauss.h>

namespace LinBox 
{
  
#if 0
	template<class T, template <class X> class Container>
		 std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
			 for(typename Container<T>::const_iterator refs =  C.begin();
					 refs != C.end() ;
					 ++refs )
				 o << (*refs) << " " ;
			 return o << std::endl;
		 }
#endif

/** \brief Repository of functions for rank modulo a prime power by elimination 
    on sparse matrices.
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
	PowerGaussDomain (const Field &F) : Father_t(F) {}

            //Copy constructor
            /// 
	PowerGaussDomain (const PowerGaussDomain &M) : Father_t(M) {}


        
            // --------------------------------------------
            // Modulo operators
        template<class Modulu>
        bool isNZero(const Modulu& a ) { return (bool)a ;}
        
#ifdef GIVARO_PRANK_OUT
#ifndef GIVARO_JRANK_OUT
#define GIVARO_JRANK_OUT
#endif
#endif

        
        template<class Ring>
        Ring MY_Zpz_bezout(Ring a, Ring b, Ring& u)
            {
                Ring u1,u3 ;
                Ring v1,v3 ;
                u1 = 1 ; u3 = a ;
                v1 = 0 ; v3 = b ;
                while (v3 != 0)
                {
                    Ring q , t1 ,t3 ;
                    q = u3 / v3 ;
                    t1 = u1 - q * v1 ; t3 = u3 - q * v3 ;
                    u1 = v1 ; u3 = v3 ; v1 = t1 ; v3 = t3 ;
                }
                u = u1 ; 
                return u3 ;
            } ;
        
        template<class Modulo, class Modulo2>
        Modulo MY_Zpz_inv (const Modulo a, const Modulo2 pp)
            {
                typedef Modulo Ring;
                if (a == 1) return a;
                if (a == -1) return a;
                Ring u, d;
                d = MY_Zpz_bezout(Ring(a),Ring(pp),u) ;
                
                if (d == -1) { u = -u; }
                if (u <0) u += pp ;
                return u ;
            }
        

        template<class Ring, class Ring2>
        Ring MY_gcd(Ring a, Ring2 b) {
            Ring q,r,d;
            Ring ma=a, mb=b;
            while (mb != 0) {
                q = ma / mb;
                r = ma - q * mb;
                ma = mb;
                mb = r;
            }
            return ma;
        };

        template<class Ring1, class Ring2>
        bool MY_divides(Ring1 a, Ring2 b) {
            return (!(b%a));
        };

// ------------------------------------------------
// Pivot Searchers and column strategy
// ------------------------------------------------
        template<class Vecteur>
        void CherchePivot( Vecteur& lignepivot, unsigned long& indcol , long& indpermut ) {
            unsigned long nj =  lignepivot.size() ;
            if (nj) {
                indpermut= lignepivot[0].first;
                if (indpermut != indcol)
                    lignepivot[0].first = (indcol);
                indcol++ ;
            } else
                indpermut = -1;
        }



        template<class Modulo, class Vecteur, class D>
        void CherchePivot(Modulo PRIME, Vecteur& lignepivot, unsigned long& indcol , long& indpermut, D& columns ) {
            typedef typename Vecteur::value_type E;
            typedef typename Field::Element F;
            long nj =  lignepivot.size() ;
            if (nj) {
                indpermut = lignepivot[0].first;
                long pp=0;
                for(;pp<nj;++pp)
                    if (! this->MY_divides(PRIME,lignepivot[pp].second) ) break;

                if (pp < nj) {
           
                    long ds = columns[ lignepivot[pp].first ],dl,p=pp,j=pp;
                    for(++j;j<nj;++j)
                        if ( ( (dl=columns[lignepivot[j].first] ) < ds ) && (! MY_divides(PRIME,lignepivot[j].second) ) ) {
                            ds = dl;
                            p = j;
                        }
                    if (p != 0) {
                        if (indpermut == (long)indcol) {
                            F ttm = lignepivot[p].second;
                            indpermut = lignepivot[p].first;
                            lignepivot[p].second = (lignepivot[0].second);
                            lignepivot[0].second = (ttm);
                        } else {
                            E ttm = lignepivot[p];
                            indpermut = ttm.first;
                            for(long m=p;m;--m)
                                lignepivot[m] = lignepivot[m-1];
                            lignepivot[0] = ttm;
                        }
                    } 
                    if (indpermut != (long)indcol)
                        lignepivot[0].first = (indcol);
                    indcol++ ;
                    for(j=nj;j--;)
                        --columns[ lignepivot[j].first ];
                } else
                    indpermut = -2;
            } else
                indpermut = -1;
        }




        template<class Modulo, class Vecteur, class De>
        void FaireElimination( Modulo MOD,
                               Vecteur& lignecourante,
                               const Vecteur& lignepivot,
                               const long& indcol,
                               const long& indpermut,
                               De& columns) {
    
//     typedef typename Vecteur::coefficientSpace F;
//     typedef typename Vecteur::value_types E;
            typedef typename Field::Element F;
            typedef typename Vecteur::value_type E;
    
	    typedef typename Signed_Trait<Modulo>::unsigned_type UModulo;

            unsigned long k = indcol - 1;
            unsigned long nj =  lignecourante.size() ;
            if (nj) {
                unsigned long j_head(0);
                for(; j_head<nj; ++j_head)
                    if (long(lignecourante[j_head].first) >= indpermut) break;
                unsigned long bjh(j_head-1);
                if ((j_head<nj) && (long(lignecourante[j_head].first) == indpermut)) {
                        // -------------------------------------------
                        // Permutation
                    if (indpermut != (long)k) {
                        if (lignecourante[0].first == k) {     
                                // non zero  <--> non zero
                            F tmp = lignecourante[0].second ;
                            lignecourante[0].second = (lignecourante[j_head].second );
                            lignecourante[j_head].second = (tmp);
                        } else {
                                // zero <--> non zero
                            E tmp = lignecourante[j_head];
                            --columns[ tmp.first ];
                            ++columns[k];
                            tmp.first = (k);
                            for(long l=j_head; l>0; l--)
                                lignecourante[l] = lignecourante[l-1];
                            lignecourante[0] = tmp;
                        } 
                        j_head = 0;
                    }
                        // -------------------------------------------
                        // Elimination
                    unsigned long npiv = lignepivot.size();
                    Vecteur construit(nj + npiv);
                        // construit : <-- ci
                        // courante  : <-- m
                        // pivot     : <-- l
                    typedef typename Vecteur::iterator Viter;
                    Viter ci = construit.begin();
                    unsigned long m=1;
                    unsigned long l(0);
                        // A[i,k] <-- A[i,k] / A[k,k]
                    lignecourante[0].second = (  ((UModulo)( ( MOD-(lignecourante[0].second) ) * ( MY_Zpz_inv( lignepivot[0].second, MOD) ) ) ) % (UModulo)MOD ) ;
                    F headcoeff = lignecourante[0].second ;
                    --columns[ lignecourante[0].first ];
           
                    unsigned long j_piv;
                    F tmp;
                    for(;l<npiv;++l)
                        if (lignepivot[l].first > k) break;
                        // for all j such that (j>k) and A[k,j]!=0
                    for(;l<npiv;++l) {
                        j_piv = lignepivot[l].first;
                            // if A[k,j]=0, then A[i,j] <-- A[i,j]
                        for (;(m<nj) && (lignecourante[m].first < j_piv);) 
                            *ci++ = lignecourante[m++];
                            // if A[i,j]!=0, then A[i,j] <-- A[i,j] - A[i,k]*A[k,j]
                        if ((m<nj) && (lignecourante[m].first == j_piv)) {
                            lignecourante[m].second = ( ((UModulo)( headcoeff  *  lignepivot[l].second  + lignecourante[m].second ) ) % (UModulo)MOD );
                            if (isNZero(lignecourante[m].second))
                                *ci++ = lignecourante[m++];
                            else 
                                --columns[ lignecourante[m++].first ];
//                         m++;
                        } else if (isNZero(tmp = ((UModulo)(headcoeff * lignepivot[l].second)) %(UModulo)MOD)) {
                            ++columns[j_piv];
                            *ci++ =  E(j_piv, tmp );
                        }
                    }
                        // if A[k,j]=0, then A[i,j] <-- A[i,j] 
                    for (;m<nj;)  
                        *ci++ = lignecourante[m++];        
        
                    construit.erase(ci,construit.end());
                    lignecourante = construit;
                } else
                        // -------------------------------------------
                        // Permutation
                    if (indpermut != (long)k) {
                        unsigned long l(0);
                        for(; l<nj; ++l)
                            if (lignecourante[l].first >= k) break;
                        if ((l<nj) && (lignecourante[l].first == k))  {
                                // non zero <--> zero
                            E tmp = lignecourante[l];
                            --columns[tmp.first ];
                            ++columns[indpermut];
                            tmp.first = (indpermut);
                            for(;l<bjh;l++)
                                lignecourante[l] = lignecourante[l+1];
                            lignecourante[bjh] = tmp;
                        } // else
                            // zero <--> zero
                    }
            }
        }
 
// ------------------------------------------------------
// Rank calculators, defining row strategy
// ------------------------------------------------------

        template<class Modulo, class BB, class D, class Container>
        void gauss_rankin(Modulo FMOD, Modulo PRIME, Container& ranks, BB& LigneA, const size_t Ni, const size_t Nj, const D& density_trait) {
        commentator.start ("Gaussian elimination with reordering modulo a prime power",
                           "PRGE", Ni);

            ranks.resize(0);

            typedef typename BB::Row Vecteur;
    
            Modulo MOD = FMOD;
#ifdef GIVARO_PRANK_OUT
            std::cerr << "Elimination mod " << MOD << std::endl;
#endif
    
            D col_density(Nj);

                // assignment of LigneA with the domain object
            size_t jj;
            for(jj=0; jj<Ni; ++jj) {
                Vecteur tmp = LigneA[jj];
                Vecteur toto(tmp.size());
                unsigned long k=0,rs=0;
                for(; k<tmp.size(); ++k) {
                    Modulo r = tmp[k].second;
                    if ((r <0) || (r >= MOD)) r = r % MOD ;
                    if (r <0) r = r + MOD ;
                    if (isNZero(r)) {
                        ++col_density[ tmp[k].first ];
                        toto[rs] =tmp[k];
                        toto[rs].second = ( r );
                        ++rs;
                    } 
                }
                toto.resize(rs);
                LigneA[jj] = toto;
//                 LigneA[jj].reactualsize(Nj);
        
            }
            Vecteur Vzer(0);

            unsigned long last = Ni-1;
            long c(0);
            unsigned long indcol(0);
            unsigned long ind_pow = 1;
            unsigned long maxout = Ni/100; maxout = (maxout<10 ? 10 : (maxout>1000 ? 1000 : maxout) );
            unsigned long thres = Ni/maxout; thres = (thres >0 ? thres : 1);


            for (unsigned long k=0; k<last;++k) {
            if ( ! (k % maxout) ) commentator.progress (k);


                unsigned long p=k;
                for(;;) {
            
                   
                    std::multimap< long, long > psizes; 
                    for(p=k; p<Ni; ++p)
                        psizes.insert( psizes.end(), std::pair<long,long>( LigneA[p].size(), p) );
                        
/*
#ifdef  GIVARO_PRANK_OUT   
                    std::cerr << "------------  ordered rows -----------" << std::endl;
                    for( std::multimap< long, long >::const_iterator iter = psizes.begin(); iter != psizes.end(); ++iter) {
                        std::cerr << (*iter).first << " : " <<  (*iter).second << std::endl;
                    }
                    std::cerr << "--------------------------------------" << std::endl;
#endif
*/

                    for( typename std::multimap< long, long >::const_iterator iter = psizes.begin(); iter != psizes.end(); ++iter) {
                        p = (*iter).second;
                        CherchePivot( PRIME, LigneA[p], indcol, c , col_density) ;
                        if (c > -2 ) break;
                    }

                    if (c > -2) break;
                    for(unsigned long ii=k;ii<Ni;++ii)
                        for(unsigned long jj=LigneA[ii].size();jj--;)
                            LigneA[ii][jj].second = ( LigneA[ii][jj].second / PRIME);
                    MOD = MOD / PRIME;
                    ranks.push_back( indcol );
                    ++ind_pow;
#ifdef GIVARO_PRANK_OUT
                    std::cerr << "Rank mod " << (unsigned long)PRIME << "^" << ind_pow << " : " << indcol << std::endl;
                    if (MOD == 1) std::cerr << "wattadayada inhere ?" << std::endl;
#endif

                }
                if (p != k) {
                    Vecteur vtm = LigneA[k];
                    LigneA[k] = LigneA[p];
                    LigneA[p] = vtm;
                }
                if (c != -1)
                    for(unsigned long l=k + 1; l < Ni; ++l)
                        FaireElimination(MOD, LigneA[l], LigneA[k], indcol, c, col_density);


//                 LigneA.write(cout << "step[" << k << "], pivot: " << c << std::endl) << endl;
    
#ifndef GIVARO_PRANK_OUT
                LigneA[k] = Vzer;
#endif
            }
            CherchePivot( PRIME, LigneA[last], indcol, c, col_density );
            while( c == -2) {
                ranks.push_back( indcol );
                for(long jj=LigneA[last].size();jj--;)
                    LigneA[last][jj].second = ( LigneA[last][jj].second / PRIME);
                MOD /= PRIME;
                CherchePivot( PRIME, LigneA[last], indcol, c, col_density );
            }
            while( MOD > 1) {
                MOD /= PRIME;
                ranks.push_back( indcol );
            }
                
//             ranks.push_back(indcol);
#ifdef GIVARO_JRANK_OUT
            std::cerr << "Rank mod " << (unsigned long)FMOD << " : " << indcol << std::endl;
#endif
	commentator.stop ("done", 0, "PRGE");

        }

        template<class Modulo, class BB, class D, class Container>
        void prime_power_rankin (Modulo FMOD, Modulo PRIME, Container& ranks, BB& SLA, const size_t Ni, const size_t Nj, const D& density_trait){
            gauss_rankin(FMOD,PRIME,ranks, SLA, Ni, Nj, density_trait);
        }


        template<class Modulo, class Matrix, template<class, class> class Container, template<class> class Alloc>
       Container<std::pair<size_t,size_t>, Alloc<std::pair<size_t,size_t> > >& operator()(Container<std::pair<size_t,size_t>, Alloc<std::pair<size_t,size_t> > >& L, Matrix& A, Modulo FMOD, Modulo PRIME) { 
           Container<size_t, Alloc<size_t> > ranks;
           prime_power_rankin( FMOD, PRIME, ranks, A, A.rowdim(), A.coldim(), std::vector<size_t>());
           L.resize( 0 ) ;
           size_t MOD = 1;
           size_t num = 0, diff;
           for( typename Container<size_t, Alloc<size_t> >::const_iterator it = ranks.begin(); it != ranks.end(); ++it) {
               diff = *it-num;
               if (diff > 0)
                   L.push_back( std::pair<size_t,size_t>(*it-num,MOD) );
               MOD *= PRIME;
               num = *it;
           }
           return L;
	}
 
    };
    


} // end of LinBox namespace

#endif  //__LINBOX_pp_gauss_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
