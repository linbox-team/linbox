// ======================================================================= //
// Givaro / Athapascan-1
// Valence computation
// Time-stamp: <09 Dec 99 12:06:18 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_DOM_NULLSPACE_H_
#define __LINBOX_DOM_NULLSPACE_H_

// #include <inttcr.C>
// #include <matrix.h>
#include <commentator.h>
#include <math.h>
#include "lin_berlekamp.C"
#include "lin_dom_gauss.C" 
#include "lin_dom_wiedemann.C"     // Symetric Wiedemann


template<class BB>
class FactorMinPolyAlg {
public:
    template<class Polynomial>
    Polynomial& operator() (Polynomial& phi, UIndice& r, ulong p1, char * Matrix_File_Name, Commentator& Comm) {

        Zpz::setZpz(p1);
        typename BB::Domain_t F(p1,1);

        BB TheMat( F );
        typename BB::Storage_t StorMat;
        long ni = StorMat.get_n_row(Matrix_File_Name);
        if (ni > StorMat.n_col())
            TheMat.init(Matrix_File_Name);
        else
            TheMat.init_transpose(Matrix_File_Name);
        TheMat.assign();
        WiedemannDom< BB > WD( Comm, &TheMat );
        WD.minpolyin(phi, r);

        for(Indice i=phi.size(); i--;)
            phi[i] = F.access(phi[i]);
  
        {

            Poly1<Zpz> P(Indeter("X"),Degree(phi.size()-1));
            for(Indice i=phi.size(); i--;)
                P[i] = phi[i];

            vector< Poly1< Zpz > > factors;
            vector<Degree> degrees;
    
            Timer tim; tim.clear(); tim.start();
            berlekamp(factors, degrees, P, p1);
            tim.stop();
    
            vector< Poly1< Zpz > >::const_iterator iF = factors.begin();
            vector< Degree >::const_iterator iD = degrees.begin();
            Degree degmax = 0;
            for(;iF!=factors.end();++iF, ++iD) {
                if (*iD > 0)
                    if (*iD > 1)
                        cerr << "(" << *iF << ")^" << *iD << endl;
                    else
                        cerr << *iF << endl;
                degmax = GIVMAX(degmax, (*iF).degree());
            }

            ++degmax;
            
            cerr << tim << endl;


            typedef typename BB::PreferredInMatrix_t           Right;
            typedef typename BB::PreferredOutMatrix_t          Left;
            Right diag_right(TheMat.n_col());
            for(Indice ll=0; ll<TheMat.n_col(); ++ll)
                diag_right[ll] = F.zero;
            Left inter(TheMat.n_row());
            Right v1(TheMat.n_col()), v2(TheMat.n_col());

            typedef TransposeDom< BB > SPBBT; SPBBT AT(&TheMat);
            typedef CompositionDom< SPBBT, BB, Left > SPBBATA; SPBBATA ATA(&AT,&TheMat,&inter);
            typedef PolynomialCompositionDom< SPBBATA, Right, Polynomial, Right> UsedBB;
            Right u(TheMat.n_col()), v(TheMat.n_col());
            vector< Right > NullSpace;
            NullSpace.reserve(degmax);
            GaussDom< typename BB::Domain_t > GD(F);

            srand48(BaseTimer::seed());

            for(iF = factors.begin(), iD = degrees.begin();iF!=factors.end();++iF, ++iD) {
                if (( *iD > 0 ) && ( (*iF).degree() > 0 )) {
                    if ( *iD > 1) 
                        cerr << "Error : matrix is not symetric, its minimum polynomial has a non linear factor" << endl;
            
                    Poly1< Zpz > f1 = P / *iF;
                    
//                     cerr << "iF : " << *iF << endl;

                    vector< typename BB::Type_t > pol(f1.degree()+2);
                    for(Indice ll=1;ll<pol.size();++ll)
                        F.assign(pol[ll], f1[ll-1].residu() );
                    pol[0] = F.zero;
                    
                    UIndice dim = (*iF).degree();
                    UIndice nb_vect = dim + 1;
                    
                    UsedBB C(&ATA, &v1, &v2, &pol, &diag_right, F);
                    NullSpace.resize(nb_vect);
                    for(Indice i=0; i<nb_vect; ++i) {
                        for(Indice ll=0; ll<TheMat.n_col(); ++ll)
                            F.random(u[ll]);
                        NullSpace[i].resize(TheMat.n_col());
                        C.Apply(NullSpace[i],u);
                    }

                    UIndice rrak;
//                     GD.gauss_rank(rrak, NullSpace);
                    GD.gauss_LUin(rrak, NullSpace);
       
//                     cerr << "Rank : " << rrak << endl;

                    if (dim < rrak) 
                        cerr << "Pb with : " << *iF << endl;
                    else
                        cerr << "To suppress : " << *iF << endl;
                }
            }
        }
        
            
    
        return phi;
    }

};
#endif __LINBOX_DOM_NULLSPACE_H_
