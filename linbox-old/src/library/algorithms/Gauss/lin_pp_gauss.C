#ifndef __LINBOX_GAUSS_C__
#define __LINBOX_GAUSS_C__
// ========================================================================= //
// (C) Givaro Team 1999
// Time-stamp: <03 Jul 00 12:28:12 Jean-Guillaume.Dumas@imag.fr> 
// ========================================================================= //

// --------------------------------------------
// Modulo operators
// int iszero(const Modulo& a ) { return !a ;};

typedef long MYLONGTYPE;


template<class Ring>
Ring MY_Zpz_bezout(Ring a, Ring b, Ring& u, Ring& v)
{
    Ring u1,u2,u3 ;
    Ring v1,v2,v3 ;
    u1 = 1 ; u2 = 0 ; u3 = a ;
    v1 = 0 ; v2 = 1 ; v3 = b ;
    while (v3 != 0)
    {
        Ring q , t1, t2 ,t3 ;
        q = u3 / v3 ;
        t1 = u1 - q * v1 ; t2 = u2 - q * v2 ; t3 = u3 - q * v3 ;
        u1 = v1 ; u2 = v2 ; u3 = v3 ; v1 = t1 ; v2 = t2 ; v3 = t3 ;
    }
    u = u1 ; v = u2 ;
    return u3 ;
} ;

template<class Modulo>
Modulo MY_Zpz_inv (const Modulo a, const Modulo pp)
{
    typedef typename Signed_Trait<Modulo>::signed_type Ring;
    if (a == 1) return a;
    if (a == -1) return a;
    Ring u,v, d ;
    d = MY_Zpz_bezout(Ring(a),Ring(pp),u,v) ;

    if (d == -1) { d = -d ;  u = -u ; v = -v ; }
//  if (d != 1) {
//      cerr << "*** division by zero in: Modulo div " << endl;
//      cerr << a << '*'<< u << '+' << v << '*' << pp << '=' << u*a+v*pp << endl ;
//      cerr << "and gcd = " << d << endl ;
//      return 0;
//    }
    if (u <0) u += pp ;
    return u ;
}


template<class Ring>
Ring MY_gcd(Ring a, Ring b) {
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

template<class Ring>
short MY_divides(Ring a, Ring b) {
    return iszero(b%a);
};

// ------------------------------------------------
// Pivot Searchers and column strategy
// ------------------------------------------------
template<class Vecteur>
void CherchePivot( Vecteur& lignepivot, long& indcol , long& indpermut ) {
    long nj =  lignepivot.size() ;
    if (nj) {
       indpermut= lignepivot[0].j();
       if (indpermut != indcol)
           lignepivot[0].change_j(indcol);
       indcol++ ;
    } else
        indpermut = -1;
}



template<class Modulo, class Vecteur, class D>
void CherchePivot(Modulo PRIME, Vecteur& lignepivot, long& indcol , long& indpermut, D& columns ) {
    typedef typename Vecteur::element E;
    typedef typename Vecteur::Type_t F;
    long nj =  lignepivot.size() ;
    if (nj) {
       indpermut = lignepivot[0].j();
       long pp=0;
       for(;pp<nj;++pp)
           if (! MY_divides(PRIME,lignepivot[pp].getvalue()) ) break;

       if (pp < nj) {
           
           long ds = columns.getvalue( lignepivot[pp].j() ),dl,p=pp,j=pp;
           for(++j;j<nj;++j)
               if ( ( (dl=columns.getvalue(lignepivot[j].j()) ) < ds ) && (! MY_divides(PRIME,lignepivot[j].getvalue()) ) ) {
                   ds = dl;
                   p = j;
               }
           if (p != 0) {
               if (indpermut == indcol) {
                   F ttm = lignepivot[p].getvalue();
                   indpermut = lignepivot[p].j();
                   lignepivot[p].affect(lignepivot[0].getvalue());
                   lignepivot[0].affect(ttm);
               } else {
                   E ttm = lignepivot[p];
                   indpermut = ttm.j();
                   for(long m=p;m;--m)
                       lignepivot[m] = lignepivot[m-1];
                   lignepivot[0] = ttm;
               }
           } 
           if (indpermut != indcol)
               lignepivot[0].change_j(indcol);
           indcol++ ;
           for(j=nj;j--;)
               columns.decr(lignepivot[j].j()); 
       } else
           indpermut = -2;
    } else
	indpermut = -1;
}



    
// ---------------------------------------------------------
// Line Elimination techniques with local permutations
// In the dense case L can be computed in place via LU
// ---------------------------------------------------------

template<class T>
void LU( vector<T>& lignecur,
         const vector<T>& lignepivot,
         long indcol ,
         long indpermut ) {

    long n =  lignecur.size() ;
    long k = indcol - 1 ;

        // permutation if one has been performed to compute the pivot
    if (indpermut != k) {
        T tmp = lignecur[k] ;
        lignecur[k] = lignecur[indpermut] ;
        lignecur[indpermut] = tmp ;
    }

    lignecur[k] = -lignecur[k] / lignepivot[k] ;
    T headcoeff = lignecur[k] ;
    long j = 0 ;
    for (; j < k ; ++j ) 
        lignecur[j] += headcoeff * lignepivot[j] ;
    for (++j ; j < n ; ++j )
        lignecur[j] += headcoeff * lignepivot[j] ;
}



template<class Modulo, class Vecteur, class De>
void FaireElimination( Modulo MOD,
                       Vecteur& lignecourante,
                       const Vecteur& lignepivot,
                       const long& indcol,
                       const long& indpermut,
                       De& columns) {
    
//     typedef typename Vecteur::coefficientSpace F;
//     typedef typename Vecteur::elements E;
    typedef typename Vecteur::Type_t F;
    typedef typename Vecteur::element E;
    

    long k = indcol - 1;
    long nj =  lignecourante.size() ;
    if (nj) {
        long j_head(0);
        for(; j_head<nj; ++j_head)
            if (lignecourante[j_head].j() >= indpermut) break;
        long bjh(j_head-1);
        if ((j_head<nj) && (lignecourante[j_head].j() == indpermut)) {
                // -------------------------------------------
                // Permutation
            if (indpermut != k) {
                if (lignecourante[0].j() == k) {     
                        // non zero  <--> non zero
                    F tmp = lignecourante[0].getvalue() ;
                    lignecourante[0].affect(lignecourante[j_head].getvalue() );
                    lignecourante[j_head].affect(tmp);
                } else {
                        // zero <--> non zero
                    E tmp = lignecourante[j_head];
                    columns.decr(tmp.j());
                    columns.incr(k);
                    tmp.change_j(k);
                    for(long l=j_head; l>0; l--)
                        lignecourante[l] = lignecourante[l-1];
                    lignecourante[0] = tmp;
                } 
                j_head = 0;
            }
                // -------------------------------------------
                // Elimination
            long npiv = lignepivot.size();
            Vecteur construit(nj + npiv);
                // construit : <-- ci
                // courante  : <-- m
                // pivot     : <-- l
            typedef typename Vecteur::iterator Viter;
            Viter ci = construit.begin();
            long m=1;
            long l(0);
                // A[i,k] <-- A[i,k] / A[k,k]
            lignecourante[0].affect(  ((MYLONGTYPE)( ( MOD-(lignecourante[0].getvalue()) ) * ( MY_Zpz_inv( lignepivot[0].getvalue(), MOD) ) ) ) % (MYLONGTYPE)MOD ) ;
            F headcoeff = lignecourante[0].getvalue() ;
            columns.decr(lignecourante[0].j());
            
            long j_piv;
            F tmp;
            for(;l<npiv;++l)
                if (lignepivot[l].j() > k) break;
                // for all j such that (j>k) and A[k,j]!=0
            for(;l<npiv;++l) {
                j_piv = lignepivot[l].j();
                    // if A[k,j]=0, then A[i,j] <-- A[i,j]
                for (;(m<nj) && (lignecourante[m].j() < j_piv);) 
                    *ci++ = lignecourante[m++];
                    // if A[i,j]!=0, then A[i,j] <-- A[i,j] - A[i,k]*A[k,j]
                if ((m<nj) && (lignecourante[m].j() == j_piv)) {
                    lignecourante[m].affect( ((MYLONGTYPE)( headcoeff  *  lignepivot[l].getvalue()  + lignecourante[m].getvalue() ) ) % (MYLONGTYPE)MOD );
                    if (! iszero(lignecourante[m].getvalue()))
                        *ci++ = lignecourante[m++];
                    else 
                        columns.decr(lignecourante[m++].j());
//                         m++;
                } else if (! iszero(tmp = ((MYLONGTYPE)(headcoeff * lignepivot[l].getvalue())) %(MYLONGTYPE)MOD)) {
                    columns.incr(j_piv);
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
            if (indpermut != k) {
                long l(0);
                for(; l<nj; ++l)
                    if (lignecourante[l].j() >= k) break;
                if ((l<nj) && (lignecourante[l].j() == k))  {
                        // non zero <--> zero
                    E tmp = lignecourante[l];
                    columns.decr(tmp.j());
                    columns.incr(indpermut);
                    tmp.change_j(indpermut);
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

template<class Modulo, class BB, class D>
void gauss_rankin(Modulo FMOD, Modulo PRIME, unsigned long& rank, BB& LigneA, long Ni, long Nj, const D& density_trait) {

    typedef typename BB::value_type Vecteur;
    
    Modulo MOD = FMOD;
    
    D col_density(Nj);

        // assignment of LigneA with the domain object
    long jj;
    for(jj=0; jj<Ni; jj++) {
        Vecteur tmp = LigneA[jj];
        Vecteur toto(tmp.size());
        long k=0,rs=0;
        for(; k<tmp.size(); ++k) {
            Modulo r = tmp[k].getvalue();
            if ((r <0) || (r >= MOD)) r = r % MOD ;
            if (r <0) r = r + MOD ;
            if (! iszero(r)) {
                col_density.incr(tmp[k].j());
                toto[rs] =tmp[k];
                toto[rs].affect( r );
                ++rs;
            } 
        }
        toto.resize(rs);
        LigneA[jj] = toto;
        
    }

    long last = Ni-1;
    long c;
    long indcol(0);

    for (long k=0; k<last;++k) {
        long p=k;
        for(;;) {
            long M=Ni,s,sl;
            for(p=k,s=LigneA[k].size(); M > k;) {
                if (s) {
                    for(long l=k+1; l < M; ++l)
                        if (((sl=LigneA[l].size()) < s) && (sl)) {
                            s = sl;
                            p = l;
                        }
                }
                CherchePivot( PRIME, LigneA[p], indcol, c , col_density) ;
                if (c > -2) break;
                else {
                    Vecteur vtm = LigneA[p];
                    LigneA[p] = LigneA[--M];
                    LigneA[M] = vtm;
                }
            }
            if (c > -2) break;
            for(long ii=k;ii<Ni;++ii)
                for(long jj=LigneA[ii].size();jj--;)
                    LigneA[ii][jj].affect( LigneA[ii][jj].getvalue() / PRIME);
            MOD = MOD / PRIME;
//             if (MOD == 1) cerr << "wattadayada inhere ?" << endl;
        }
        if (p != k) {
            Vecteur vtm = LigneA[k];
            LigneA[k] = LigneA[p];
            LigneA[p] = vtm;
        }
        if (c != -1)
            for(long l=k + 1; l < Ni; ++l)
                FaireElimination(MOD, LigneA[l], LigneA[k], indcol, c, col_density);
    }
    CherchePivot( LigneA[last], indcol, c );
    
    rank = indcol;
}

template<class Modulo, class BB, class D>
void prime_power_rankin (Modulo FMOD, Modulo PRIME, unsigned long& rank, BB& SLA, long Ni, long Nj, const D& density_trait){
    gauss_rankin(FMOD,PRIME,rank, SLA, Ni, Nj, density_trait);
}



#endif __LINBOX_GAUSS_C__
