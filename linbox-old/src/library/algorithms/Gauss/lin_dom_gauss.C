#ifndef __LINBOX_DOM_GAUSS_C__
#define __LINBOX_DOM_GAUSS_C__
// ========================================================================= //
// (C) The Linbox Group 1999
// Calcul de rang par la méthode de Gauss pivot par ligne, sur matrice creuse
// Time-stamp: <06 Sep 00 15:29:26 Jean-Guillaume.Dumas@imag.fr> 
// ========================================================================= //

#include <commentator.C>


template <class Domain>
class GaussDom {
private:
    Domain _domain;
    Commentator _comm;
public:

    typedef          Domain           Domain_t;
    typedef          GaussDom<Domain> Self_t;

        //-- Default cstors:
    GaussDom() {};
    GaussDom(const Domain& D, const Commentator& Co = 0) : _domain(D), _comm(Co) {}
    
        //-- Cstor of recopy: compiler's generated
    GaussDom(const Self_t& M, const Commentator& Co = 0) : _domain(M._domain), _comm(Co) {}

    const Domain& getdomain() { return _domain; }

protected:
    


//-----------------------------------------
// Sparse elimination using a pivot row :
// lc <-- lc - lc[k]/lp[0] * lp 
// D is the number of elements per column
//   it is updated and used for reordering
// Vecteur is a vector of Pair (lin_pair.h)
//-----------------------------------------
template<class Vecteur, class D>
void FaireElimination( Vecteur& lignecourante,
                       const Vecteur& lignepivot,
                       const long& indcol,
                       const long& indpermut,
                       D& columns) {
    
    typedef typename Vecteur::value_type E;
    typedef typename Vecteur::Type_t     Type_t;

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
                    Type_t tmp = lignecourante[0].getvalue() ;
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
                // construit : <-- j
                // courante  : <-- m
                // pivot     : <-- l
            long j=0;
            long m=j_head+1;
                // A[i,k] <-- - A[i,k] / A[k,k]

//             lignecourante[j_head].affect(  -lignecourante[j_head].getvalue() / lignepivot[0].getvalue() ) ;
//             Type_t headcoeff = lignecourante[j_head].getvalue() ;
            Type_t headcoeff;
            _domain.divin( _domain.neg(headcoeff, lignecourante[j_head].getvalue()), lignepivot[0].getvalue());
            columns.decr(lignecourante[j_head].j());
        
                // if A[k,j]=0, then A[i,j] <-- A[i,j]
            while(j < j_head) {
                construit[j] = lignecourante[j];
                j++;
            }
            long j_piv;

            long l(0);
            for(; l<npiv; l++)
                if (lignepivot[l].j() > k) break;
                    // for all j such that (j>k) and A[k,j]!=0
                while(l < npiv) {
                    j_piv = lignepivot[l].j();
                        // if A[k,j]=0, then A[i,j] <-- A[i,j]
                     while ((m<nj) && (lignecourante[m].j() < j_piv))
                        construit[j++] = lignecourante[m++];
                        // if A[i,j]!=0, then A[i,j] <-- A[i,j] - A[i,k]*A[k,j]
                    if ((m<nj) && (lignecourante[m].j() == j_piv)) {
                        Type_t tmp; 
//                         tmp.muladd(headcoeff,lignepivot[l].getvalue(),lignecourante[m].getvalue());
                        _domain.axpy(tmp, headcoeff,lignepivot[l].getvalue(),lignecourante[m].getvalue());
                        if (! _domain.iszero(tmp)) {
                            lignecourante[m].affect( tmp );
                            construit[j++] = lignecourante[m++];
                        } else
                            columns.decr(lignecourante[m++].j());
//                             m++;
                    
                    } else {
//                         Type_t tmp = headcoeff * lignepivot[l].getvalue();
                        Type_t tmp; _domain.mul(tmp,headcoeff,lignepivot[l].getvalue());
                        if (! _domain.iszero(tmp)) {
                            columns.incr(j_piv);
                            construit[j++] = E(j_piv, tmp );
                        }
                    }
                    l++;
                }
        
                // if A[k,j]=0, then A[i,j] <-- A[i,j] 
            while (m<nj)
                construit[j++] = lignecourante[m++];
        
            construit.resize(j);
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

//-----------------------------------------
// Sparse elimination using a pivot row :
// lc <-- lc - lc[k]/lp[0] * lp 
// No density update
// Vecteur is a vector of Pair (lin_pair.h)
//-----------------------------------------
template<class Vecteur>
void FaireElimination( Vecteur& lignecourante,
                       const Vecteur& lignepivot,
                       const long& indcol,
                       const long& indpermut) {
    
    typedef typename Vecteur::Type_t Type_t;
    typedef typename Vecteur::value_type E;

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
                    Type_t tmp = lignecourante[0].getvalue() ;
                    lignecourante[0].affect(lignecourante[j_head].getvalue() );
                    lignecourante[j_head].affect(tmp);
                } else {
                        // zero <--> non zero
                    E tmp = lignecourante[j_head];
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
                // construit : <-- j
                // courante  : <-- m
                // pivot     : <-- l
            long j=0;
            long m=j_head+1;
                // A[i,k] <-- - A[i,k] / A[k,k]

//             lignecourante[j_head].affect(  -lignecourante[j_head].getvalue() / lignepivot[0].getvalue() ) ;
//             Type_t headcoeff = lignecourante[j_head].getvalue() ;
            Type_t headcoeff;
            _domain.divin( _domain.neg(headcoeff, lignecourante[j_head].getvalue()), lignepivot[0].getvalue());
        
                // if A[k,j]=0, then A[i,j] <-- A[i,j]
            while(j < j_head) {
                construit[j] = lignecourante[j];
                j++;
            }
            long j_piv;

            long l(0);
            for(; l<npiv; l++)
                if (lignepivot[l].j() > k) break;
                    // for all j such that (j>k) and A[k,j]!=0
                while(l < npiv) {
                    j_piv = lignepivot[l].j();
                        // if A[k,j]=0, then A[i,j] <-- A[i,j]
                     while ((m<nj) && (lignecourante[m].j() < j_piv))
                        construit[j++] = lignecourante[m++];
                        // if A[i,j]!=0, then A[i,j] <-- A[i,j] - A[i,k]*A[k,j]
                    if ((m<nj) && (lignecourante[m].j() == j_piv)) {
                        Type_t tmp; 
//                         tmp.muladd(headcoeff,lignepivot[l].getvalue(),lignecourante[m].getvalue());
                        _domain.axpy(tmp, headcoeff,lignepivot[l].getvalue(),lignecourante[m].getvalue());
                        if (! _domain.iszero(tmp)) {
                            lignecourante[m].affect( tmp );
                            construit[j++] = lignecourante[m++];
                        } else
                            ++m;
                    
                    } else {
//                         Type_t tmp = headcoeff * lignepivot[l].getvalue();
                        Type_t tmp; _domain.mul(tmp,headcoeff,lignepivot[l].getvalue());
                        if (! _domain.iszero(tmp)) {
                            construit[j++] = E(j_piv, tmp );
                        }
                    }
                    l++;
                }
        
                // if A[k,j]=0, then A[i,j] <-- A[i,j] 
            while (m<nj)
                construit[j++] = lignecourante[m++];
        
            construit.resize(j);
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
                    tmp.change_j(indpermut);
                    for(;l<bjh;l++)
                        lignecourante[l] = lignecourante[l+1];
                    lignecourante[bjh] = tmp;
                } // else
                    // zero <--> zero
            }
    }
}



//-----------------------------------------
// Dense elimination using a pivot row :
// lc <-- lc - lc[k]/lp[0] * lp 
// Computing only for k to n (and not 0 to n in LU)
//-----------------------------------------
template<class Vect_t>
void Upper( Vect_t& lignecur,
         const Vect_t& lignepivot,
         long indcol ,
         long indpermut ) {

    long n =  lignecur.size() ;
    long k = indcol - 1 ;

        // permutation if one has been performed to compute the pivot
    if (indpermut != k) {
        typename Vect_t::value_type tmp = lignecur[k] ;
        lignecur[k] = lignecur[indpermut] ;
        lignecur[indpermut] = tmp ;
    }

    typename Vect_t::value_type headcoeff;
    _domain.divin( _domain.neg(headcoeff, lignecur[k]), lignepivot[k]);
// LU in place
    _domain.assign( lignecur[k], _domain.zero);
    for (long j=k ; ++j < n ;)
        _domain.axpyin(lignecur[j],headcoeff,lignepivot[j]) ;
}


//-----------------------------------------
// Dense elimination using a pivot row :
// lc <-- lc - lc[k]/lp[0] * lp 
//-----------------------------------------
template<class Vect_t>
void LU( Vect_t& lignecur,
         const Vect_t& lignepivot,
         long indcol ,
         long indpermut ) {

    long n =  lignecur.size() ;
    long k = indcol - 1 ;

        // permutation if one has been performed to compute the pivot
    if (indpermut != k) {
        typename Vect_t::value_type tmp = lignecur[k] ;
        lignecur[k] = lignecur[indpermut] ;
        lignecur[indpermut] = tmp ;
    }

    typename Vect_t::value_type headcoeff;
// LU in place
    _domain.div(headcoeff, lignecur[k], lignepivot[k]);
    _domain.assign(lignecur[k], headcoeff);
    for (long j = k ; ++j < n ; )
        _domain.axmyin(lignecur[j],headcoeff,lignepivot[j]) ;
}




//------------------------------------------
// Looking for a non-zero pivot in a row 
// Using the column density for reordering
// Pivot is chosen as to :
// 1. Row density is minimum
// 2. Column density is minimum for this row
//------------------------------------------
template<class Vecteur, class D>
void SparseCherchePivot( Vecteur& lignepivot, long& indcol , long& indpermut, D& columns ) {
    typedef typename Vecteur::value_type             E;    
    typedef typename Vecteur::Type_t                 Type_t;    
//     typedef typename Vecteur::element          E;
//     typedef typename Vecteur::coefficientSpace Type_t;
    long nj =  lignepivot.size() ;
    if (nj) {
       indpermut= lignepivot[0].j();


       long ds = columns.decr(indpermut),dl,p=0;
       for(long j=1;j<nj;++j) {
           if ( (dl=columns.decr(lignepivot[j].j()) ) < ds ) {
               ds = dl;
               p = j;
           }
       }
       if (p != 0) {
           if (indpermut == indcol) {
               Type_t ttm = lignepivot[p].getvalue();
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
// no need to decrement/increment, already done during the search
           lignepivot[0].change_j(indcol);
       indcol++ ;
    } else
	indpermut = -1;

//     _comm.report(LVL_IMP, INTERNAL_DESCRIPTION) << "Pivot: " << indpermut << endl;
}


//------------------------------------------
// Looking for a non-zero pivot in a row 
// No reordering
//------------------------------------------
template<class Vecteur>
void SparseCherchePivot( Vecteur& lignepivot, long& indcol , long& indpermut ) {
    long nj =  lignepivot.size() ;
    if (nj) {
       indpermut= lignepivot[0].j();
       if (indpermut != indcol)
           lignepivot[0].change_j(indcol);
       indcol++ ;
    } else
	indpermut = -1;
}

//------------------------------------------
// Looking for a non-zero pivot in a row 
// Dense search
//------------------------------------------
template<class Vect_t>
void CherchePivot(Vect_t& lignepivot, unsigned long& k, long& indpermut ) {
     long n =  lignepivot.size();
     long j = k;
     for (; j< n ; ++j )
        if (! _domain.iszero(lignepivot[j])) break ;
     if (j == n )
       indpermut = -1 ;
     else {
       indpermut = j ;
       if (indpermut != k) {
         typename Vect_t::value_type tmp = lignepivot[k] ;
         lignepivot[k] = lignepivot[j] ;
         lignepivot[j] = tmp ;
       }
       ++k;
     }
}


public:


//--------------------------------------------
// Sparse in place
// Gaussian elimination with reordering
// Uses : SparseCherchePivot
//        FaireElimination(..., density)
//--------------------------------------------
template<class SparseM, class D>
void gauss_rankin(unsigned long& rank, SparseM& LigneA, const D& density_trait) {
	gauss_rankin(rank, LigneA, LigneA.n_row(),  LigneA.n_col(), density_trait);
}

template<class SparseM, class D>
void gauss_rankin(unsigned long& rank, SparseM& LigneA, unsigned long Ni, unsigned long Nj, const D& density_trait) {
//    typedef typename SparseM::Row_t                  Vecteur;
    typedef typename SparseM::value_type                  Vecteur;
    typedef typename Vecteur::value_type             E;    
    typedef typename Vecteur::Type_t                 Type_t;    
// Requirements : LigneA is an array of sparse rows
// In place (LigneA is modified)
// With reordering (D is a density type. Density is allocated here)
//    long Ni = LigneA.n_row(), Nj = LigneA.n_col();
    _comm.start("Gauss Reordering",LVL_NORMAL,INTERNAL_DESCRIPTION) 
        << Ni << " x " << Nj << endl;

    Vecteur Vzer(0);

        // allocation of the column density
    D col_density(Nj);

        // assignment of LigneA with the domain object
    long jj;
    for(jj=0; jj<Ni; jj++) {
        Vecteur tmp = LigneA[jj];
//         Vecteur toto(tmp.size());
//         long rs=0;
	long k=0;
        for(; k<tmp.size(); k++) {
//             Type_t r;
//             _domain.assign(r,tmp[k].getvalue());
//             if (! _domain.iszero(r)) {
            col_density.incr(tmp[k].j());
//                 toto[rs++] = E(tmp[k].j(), r); 
//             }
        }
//         toto.resize(rs);
//         LigneA[jj] = toto;
    }

    long last = Ni-1;
    long c;
    long indcol(0);
    

        // Elimination steps with reordering
    for (long k=0; k<last;++k) {
        if ( ! (k % 1000) ) _comm.progress("row steps",LVL_IMP,k,Ni);

        long l,p=k,s=LigneA[k].size(),sl;
        if (s) {
            for(l=k+1; l < Ni; ++l)
                if (((sl=LigneA[l].size()) < s) && (sl)) {
                    s = sl;
                    p = l;
                }
            if (p != k) {
                Vecteur vtm = LigneA[k];
                LigneA[k] = LigneA[p];
                LigneA[p] = vtm;
            }       
            SparseCherchePivot( LigneA[k], indcol, c , col_density) ;
            if (c != -1)
                for(l=k + 1; l < Ni; ++l)
                    FaireElimination(LigneA[l], LigneA[k], indcol, c, col_density);
		LigneA[k] = Vzer;
        }
    }
    SparseCherchePivot( LigneA[last], indcol, c );
    
    rank = indcol;

    _comm.stop(LVL_NORMAL,PARTIAL_RESULT) 
        << "Rank : " << rank
        << " over GF(" << _domain.size() << ")" << endl;
}
    

//--------------------------------------------
// Copy (NOT IN PLACE) and Sparse
// Gaussian elimination without reordering
// Uses : SparseCherchePivot
//        FaireElimination
//--------------------------------------------
template<class SparseM>
void gauss_rank(unsigned long& rank, const SparseM& SLA) {
// Requirements : SLA is an array of sparse rows
// WARNING : NOT IN PLACE, THERE IS A COPY.
// Without reordering (Pivot is first non-zero in row)
    long Ni = SLA.n_row(), Nj = SLA.n_col();
//    typedef typename SparseM::Row_t                  Vecteur;
    typedef typename SparseM::value_type               Vecteur;
    typedef typename Vecteur::value_type               E;    
    typedef typename Vecteur::Type_t       Type_t;    
    
    Vecteur * LigneA = new Vecteur[Ni];
    long jj;
    for(jj=0; jj<Ni; jj++) {
        Vecteur tmp = SLA[jj];
        Vecteur toto(tmp.size());
        long rs=0;
	long k=0;
        for(; k<tmp.size(); k++) {
            Type_t r;
            _domain.assign(r,tmp[k].getvalue());
            if (! _domain.iszero(r)) {
                toto[rs++] = E(tmp[k].j(), r); 
            }
        }
        toto.resize(rs);
        LigneA[jj] = toto;
    }

    long last = Ni-1;
    long c;
    long indcol(0);
    
    for (long k=0; k<last;++k) {
        long l,p=k,s=LigneA[k].size(),sl;
        if (s) {
            SparseCherchePivot( LigneA[k], indcol, c) ;
            if (c != -1)
                for(l=k + 1; l < Ni; ++l)
                    FaireElimination(LigneA[l], LigneA[k], indcol, c);
        }
    }
    SparseCherchePivot( LigneA[last], indcol, c );
    
    delete [] LigneA;
        
    rank = indcol;
}

template<class SparseM>
void gauss_rankin(unsigned long& rank, SparseM& LigneA) {
	gauss_rankin(rank, LigneA, LigneA.n_row(), LigneA.n_col());
}

template<class SparseM>
void gauss_rankin(unsigned long& rank, SparseM& LigneA, unsigned long Ni, unsigned long Nj ) {
// Requirements : SLA is an array of sparse rows
// IN PLACE.
// Without reordering (Pivot is first non-zero in row)
//     long Ni = SLA.n_row(), Nj = SLA.n_col();
//    long Ni = LigneA.n_row(), Nj = LigneA.n_col();
    _comm.start("Gauss",LVL_NORMAL,INTERNAL_DESCRIPTION) 
        << Ni << " x " << Nj << endl;

//    typedef typename SparseM::Row_t                  Vecteur;
    typedef typename SparseM::value_type               Vecteur;
    typedef typename Vecteur::value_type               E;    
    typedef typename Vecteur::Type_t       Type_t;    
    
//     Vecteur * LigneA = new Vecteur[Ni];
//     long jj;
//     for(jj=0; jj<Ni; jj++) {
//         Vecteur tmp = SLA[jj];
//         Vecteur toto(tmp.size());
//         long rs=0;
// 	long k=0;
//         for(; k<tmp.size(); k++) {
//             Type_t r;
//             _domain.assign(r,tmp[k].getvalue());
//             if (! _domain.iszero(r)) {
//                 toto[rs++] = E(tmp[k].j(), r); 
//             }
//         }
//         toto.resize(rs);
//         LigneA[jj] = toto;
//     }

    long last = Ni-1;
    long c;
    long indcol(0);
    
    for (long k=0; k<last;++k) {
        if ( ! (k % 1000) ) _comm.progress("row steps",LVL_IMP,k,Ni);
        long l,p=k,s=LigneA[k].size(),sl;
        if (s) {
            SparseCherchePivot( LigneA[k], indcol, c) ;
            if (c != -1)
                for(l=k + 1; l < Ni; ++l)
                    FaireElimination(LigneA[l], LigneA[k], indcol, c);
        }
    }
    SparseCherchePivot( LigneA[last], indcol, c );
    
    rank = indcol;
    _comm.stop(LVL_NORMAL,PARTIAL_RESULT) 
        << "Rank : " << rank
        << " over GF(" << _domain.size() << ")" << endl;
}


//--------------------------------------------
// Dense in place
// Gaussian elimination without reordering
// Uses : CherchePivot
//        LU
//--------------------------------------------
template<class VoV>
long& gauss_Uin(unsigned long& rank, VoV & A) {
// Requirements : A is an array of rows
// In place (A is modified)
// Without reordering (Pivot is first non-zero in row)
    long Ni = A.size();
    long last = Ni-1;
    long c;
    unsigned long indcol(0);
    for(long k=0; k<last; ++k) {
        CherchePivot(A[k], indcol, c);
        if (c != -1)
            for(long l=k + 1; l < Ni; ++l)
                Upper(A[l], A[k], indcol, c);
    }
    CherchePivot(A[last], indcol, c);
    rank = indcol;
}

//--------------------------------------------
// Dense in place
// LU factorization without reordering
// Uses : CherchePivot
//        LU
//--------------------------------------------
template<class VoV>
long& gauss_LUin(unsigned long& rank, VoV & A) {
// Requirements : A is an array of rows
// In place (A is modified)
// Without reordering (Pivot is first non-zero in row)
    long Ni = A.size();
    long last = Ni-1;
    long c;
    unsigned long indcol(0);
    for(long k=0; k<last; ++k) {
        CherchePivot(A[k], indcol, c);
        if (c != -1)
            for(long l=k + 1; l < Ni; ++l)
                LU(A[l], A[k], indcol, c);
    }
    CherchePivot(A[last], indcol, c);
    rank = indcol;
}



};



#endif __LINBOX_DOM_GAUSS_C__
