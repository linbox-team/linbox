#ifndef __LINBOX_BERLEKAMP_C__
#define __LINBOX_BERLEKAMP_C__
// ================================================================================== //
// Givaro / Athapascan-1
// Factorisation de polynomes : méthode de Berlekamp
// Time-stamp: <07 Sep 99 17:30:58 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================================== //
#include <poly1.h>
#include <vector.h>
#include "lin_gauss.C"


// ------------------------------------------------
// Here we define a threshold in order to switch
// from a small prime and a big prime algorithm.
// see :
// S.R. Czapor, K.O. Geddes, and G. Labahn, 
// Algorithms for Computer Algebra.
// Kluwer Academic Publishers, 1992.
// pp 359.
// ------------------------------------------------
#ifndef _BERLEKAMP_THRESHOLD_
#define _BERLEKAMP_THRESHOLD_ 1000
#endif _BERLEKAMP_THRESHOLD_

// ------------------------------------------------
// Null-spaces
template<class T, class MM>
void null_space (vector< Poly1<T> >& Res, const Poly1<T>& P, MM Mod) {

    Indice n=P.degree();
    Indice rank;
    T pivot,un(one),Zero(zero);
    Poly1<T> pxun(P.indeter(),Degree(1),un),pun(P.indeter(),Degree(0),un);
    
    vector<T> * LigneA = new vector<T> [n];
    vector<T> sfq(n);
    sfq[0] = un;
    for(Indice j=1; j<n; j++){
        sfq[j] = Zero;
    }
    LigneA[0] = sfq;
        // Initialisation : formQ 
    if (Mod < _BERLEKAMP_THRESHOLD_) {
        vector<T> rfq(n);
        Indice row;
        sfq[0] = un;
        for(Indice j=1; j<n; ++j)
            sfq[j] = Zero;
        Indice endi = (n-1)*Mod;
        for(Indice m=1; m<=endi; ++m){    
            rfq = sfq;

            sfq[0] = -rfq[n-1] * P[0];        
            for(Indice j=1; j<n; ++j)
                sfq[j] = rfq[j-1] + -rfq[n-1] * P[j];

            if ((m%Mod) == 0) {       // les x^ip, autres lignes
                row=Indice(m/Mod);
                LigneA[row]= sfq;
            }
        }
    } else
    {

        Poly1<T> xq= powmod(pxun,Mod,P);
        Poly1<T> xpk = pun;

        for(Indice kb=1;kb<n;kb++) {
            xpk = (xq * xpk) % P;        
            Indice dx=xpk.degree();
            Indice j=n;
            for(--j;j>dx;--j)
                sfq[j] = 0;
            for(;j>=0;--j)
                sfq[j] = xpk[j];
            LigneA[kb] = sfq;
        }
    }
        //  A-I
    for(Indice i=0; i<n; i++)      
        LigneA[i][i] -= un;

        // end formQ

        // Null spaces
    Indice * c = new Indice [n];

    Indice indcol(0),k;
    
    for (k=0; k<n;k++) {        
        CherchePivot( LigneA[k], indcol, c[k] ) ;
        if (c[k] != -1)
            for(Indice l=k + 1; l < n; l++)
                LU(LigneA[l], LigneA[k], indcol, c[k]);
    }
    
    rank = n - indcol;

    // end Null spaces
    // transformation des vecteurs en polynomes

    Indice r(0);
    vector< Poly1<T> > tab1(rank);
    Poly1<T> cour;


    for(k=0; k<n; k++){         
        if (c[k] == -1) {
            Indice jv(0);
            cour = Poly1<T>(P.indeter(),Degree(k), un);
            for(Indice alpha = 0; alpha<k; alpha++)
                if (c[alpha] != -1)
                    cour[alpha] = LigneA[k][jv++];
            tab1[r++] = cour;
        }
    }
        // end transformation en polynomes

    ::delete [] LigneA;
	
    // Pgcds
    vector< Poly1<T> > factors(rank + 1);
    factors[1] = P;
    Indice nb(1),nbfor,jj,j;
    Poly1<T> u,g,vp;
 
    if (Mod < _BERLEKAMP_THRESHOLD_) {
        r = Indice(1);
        MM sm;
        while (nb<rank) {
            nbfor = nb;
            for(j=1; j<=nbfor; j++) {
                vp = tab1[r];
                u = factors[j];
                for(sm=MM(zero); sm<Mod; ++sm) {
                    vp[0] = T(sm);
                    g = gcd(vp,u);
                    if ((g!=pun) && (g!=u)) {
                        u = u / g;
                        factors[j] = u;
                        factors[++nb] = g;
                    }
                    if (nb>=rank) break;
                }      
                if (nb>=rank) break;
            }
            r++;
        }
    } else {

        srand48(BaseTimer::seed());
        while (nb<rank) {
            nbfor = nb;
            for(j=1; j<=nbfor; j++) {
                u = factors[j];
                vp = Poly1<T>(zero);
                for(r=0; r<rank; r++) 
                    vp += tab1[r] * T(drand48()*Mod);
                if (Mod != 2)
                    vp = powmod(vp,(Mod-1)/2,u) - pun;
                g = gcd(vp,u);
                if ((g!=pun) && (g!=u)) {
                    u = u / g;
                    factors[j] = u;
                    for(jj=++nb; jj>j;jj--)
                        factors[jj] = factors[jj-1];
                    factors[j] = g;
                    j--;
                }
                if (nb>=rank) break;
            }
        }
    }
    
        // end pgcd's
    Res = factors;
}
 
template<class T, class MM>
void berlekamp (vector<Poly1<T> >& factors, vector<Degree>& degrees, const Poly1<T>& ent, MM Mod) {

    vector< Poly1<T> > tab;
    Poly1<T> * Fact = new Poly1<T>[ent.degree()+1];
    int nb;
    sqrfree(ent, nb, Fact);
    for(Degree count=0; count<=ent.degree(); count++) {
        if (Fact[count].degree() > 1) {
            null_space(tab, Fact[count], Mod);
            vector< Poly1<T> >::iterator ti=tab.begin();
            for(++ti;ti!=tab.end();++ti) {
                factors.push_back( *ti );
                degrees.push_back( count );
            }
        } else if (Fact[count].degree() > deginfty) {
            factors.push_back( Fact[count] );
            degrees.push_back( count );
        }
    }
    ::delete [] Fact;
}

#endif __LINBOX_BERLEKAMP_C__
