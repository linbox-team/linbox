// ================================================================
// LinBox Project 1999
// Symmetrizing iterator (for rank computations)
// Time-stamp: <26 May 00 18:46:24 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __Symmetrize_Lanczos_H__
#define __Symmetrize_Lanczos_H__
#include <LinBox/lin_lanczos_symmetric_bbit.h>
#include <LinBox/lin_bb_aat.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t, class RandIter = Random>
class Symmetrize_Lanczos : public Symmetric_Lanczos< AAT<BlackBoxDomain, Vecteur> , Vecteur > {
public:
    Symmetrize_Lanczos() {} 

    Symmetrize_Lanczos(BlackBoxDomain * D, const Vecteur& u0) 
            : Symmetric_Lanczos< AAT<BlackBoxDomain, Vecteur> , Vecteur>( new AAT<BlackBoxDomain, Vecteur>(D), u0) {}
    
    Symmetrize_Lanczos(BlackBoxDomain * D, RandIter& g ) 
            : Symmetric_Lanczos< AAT<BlackBoxDomain, Vecteur> , Vecteur>( new AAT<BlackBoxDomain, Vecteur>(D) ,g) {}

};

#endif // __Symmetrize_Lanczos_H__
