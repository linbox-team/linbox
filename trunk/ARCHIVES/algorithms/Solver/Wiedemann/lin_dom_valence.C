// ======================================================================= //
// Givaro / Athapascan-1
// Valence computation
// Time-stamp: <09 Dec 99 12:03:47 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_DOM_VALENCE_C__
#define __LINBOX_DOM_VALENCE_C__

#include <inttcr.C>
#include <commentator.h>
#include <math.h>
#include "lin_dom_wiedemann.C"     // Symetric Wiedemann
#define MAXPRIMESIZE 65536
#define MINLOWERSIZE 5000
#define valBOUND 3

  
template<class BB>
class OneValenceAlg {
protected:
    Commentator _Comm;
public:
    OneValenceAlg() : _Comm(PRINT_NOTHING,PRINT_NOTHING) {}
    OneValenceAlg(const Commentator& C) : _Comm(C) {}

    typename BB::Type_t& operator() (typename BB::Type_t& valence, UIndice& rank, ulong p1, char * Matrix_File_Name) {
        typename BB::Domain_t F(p1,1);
        BB MF( F );
        typename BB::Storage_t TheMat;
        long ni = TheMat.get_n_row(Matrix_File_Name);
        if (ni > TheMat.n_col())
            MF.init(Matrix_File_Name);
        else
            MF.init_transpose(Matrix_File_Name);
        MF.assign();
        WiedemannDom< BB > WD( _Comm, &MF );
        WD.valencein(valence, rank);
        return valence = F.access(valence);
    }

};

template<class BB>
class OneMinPolyAlg {
protected:
    Commentator _Comm;
public:
    OneMinPolyAlg() : _Comm(PRINT_NOTHING,PRINT_NOTHING) {}
    OneMinPolyAlg(const Commentator& C) : _Comm(C) {}
public:
    template<class Polynomial>
    Polynomial& operator() (Polynomial& phi, UIndice& rank, ulong p1, char * Matrix_File_Name) {
        typename BB::Domain_t F(p1,1);
        BB MF( F );
        typename BB::Storage_t TheMat;
        long ni = TheMat.get_n_row(Matrix_File_Name);
        if (ni > TheMat.n_col())
            MF.init(Matrix_File_Name);
        else
            MF.init_transpose(Matrix_File_Name);
        MF.assign();
        WiedemannDom< BB > WD( _Comm, &MF );
        WD.minpolyin(phi, rank);

        for(Indice i=phi.size(); i--;)
            phi[i] = F.access(phi[i]);
        return phi;
    }
};

template<class IntAbsBB>
class CassiniBound {
protected:
    Commentator _Comm;
public:
    CassiniBound() : _Comm(PRINT_NOTHING,PRINT_NOTHING) {}
    CassiniBound(const Commentator& C) : _Comm(C) {}
public:
    double& operator() (double& bound, char * Matrix_File_Name) {
            
        _Comm.start("Cassini",LVL_NORMAL,INTERNAL_DESCRIPTION) << endl;
        IntAbsBB GershB;
        GershB.init(Matrix_File_Name);
        GershB.assign();
        bound = (GershB.square_cassini()).Int2double();
        _Comm.stop(LVL_IMP,PARTIAL_RESULT) << "disks bound : " << bound << endl;
        return bound;
    }

};

template<class BB, class IntAbsBB>
class nbValenceAlg : public OneValenceAlg<BB> {
public:
    nbValenceAlg() : OneValenceAlg<BB>() {}
    nbValenceAlg(const Commentator& C) : OneValenceAlg<BB>(C) {}
    ulong& operator() (ulong& nbloops, typename BB::Type_t& valence, UIndice& rank, ulong p1, char * Matrix_File_Name) {

        double bound;
        CassiniBound<IntAbsBB>(_Comm).operator()(bound,Matrix_File_Name);
        OneValenceAlg<BB>::operator()(valence,rank,p1,Matrix_File_Name);
        return nbloops = (unsigned long)ceil(  (double)rank * log(2*bound) / log(MAXPRIMESIZE-MINLOWERSIZE) );

    }

};

template<class BB, class IntAbsBB>
class nbMinPolyAlg : public OneMinPolyAlg<BB> {
public:
    nbMinPolyAlg() : OneMinPolyAlg<BB>() {}
    nbMinPolyAlg(const Commentator& C) : OneMinPolyAlg<BB>(C) {}
    template<class Polynomial>
    Polynomial& operator() (Polynomial& phi, ulong& nbloops, ulong p1, char * Matrix_File_Name) {
            
        double bound;
        CassiniBound<IntAbsBB>(_Comm).operator()(bound,Matrix_File_Name);
        UIndice rank;
        OneMinPolyAlg<BB>::operator()(phi,rank,p1,Matrix_File_Name);
        nbloops = (unsigned long)ceil(  (double)rank * log(2*bound) / log( MAXPRIMESIZE-MINLOWERSIZE) );
        _Comm.report(LVL_IMP,PARTIAL_RESULT) << "loops : " << nbloops << endl;
        return phi;
    }

};


template<class BB, class IntAbsBB>
class IntValenceAlg : public nbValenceAlg<BB,IntAbsBB> {
public:
    IntValenceAlg() : nbValenceAlg<BB,IntAbsBB>() {}
    IntValenceAlg(const Commentator& C) : nbValenceAlg<BB,IntAbsBB>(C) {}
    template<class Int>
    Int& operator() (Int& V, char * Matrix_File_Name) {

        _Comm.start("IntV",LVL_NORMAL,INTERNAL_DESCRIPTION) << endl;
        UIndice rank;
        ulong p;
        do
            p = prevprime(MAXPRIMESIZE - (ulong)lrand48()%(ulong)MINLOWERSIZE).Int2long();
        while (! isprime(p)) ;
        typename BB::Type_t valence;
        ulong nbloops;
        nbValenceAlg<BB,IntAbsBB>::operator()(nbloops,valence,rank,p,Matrix_File_Name);
        V = Int(valence);
        Integer prod = p;

        for(;--nbloops;) {
            do
                p = prevprime(MAXPRIMESIZE - (ulong)lrand48()%(ulong)MINLOWERSIZE).Int2long();
            while ((! isprime(p)) && (! isone(gcd(prod,p)))) ;
            OneValenceAlg<BB>::operator()(valence,rank,p,Matrix_File_Name);
            tcr(V,prod,valence,p);
            prod *= p;
        }
        if (V > prod/2) V -= prod;
        _Comm.stop(LVL_NORMAL,PARTIAL_RESULT) << "Valence : " << V << endl;
        return V;
    }

};


template<class BB, class IntAbsBB>
class IntMinPolyAlg : public nbMinPolyAlg<BB,IntAbsBB> {
public:
    IntMinPolyAlg() : nbMinPolyAlg<BB,IntAbsBB>() {}
    IntMinPolyAlg(const Commentator& C) : nbMinPolyAlg<BB,IntAbsBB>(C) {}
    template<class Polynomial>
    Polynomial& operator() (Polynomial& P, char * Matrix_File_Name) {
        _Comm.start("IntMP",LVL_NORMAL,INTERNAL_DESCRIPTION) << endl;
        ulong p;
        do
            p = prevprime(MAXPRIMESIZE - (ulong)lrand48()%(ulong)MINLOWERSIZE).Int2long();
        while (! isprime(p)) ;
        typename WiedemannDom<BB>::InternalPolynomial_t phi;
        ulong nbloops;
        nbMinPolyAlg<BB,IntAbsBB>::operator()(phi,nbloops,p,Matrix_File_Name);
        P = Polynomial(phi.size());
        for(Indice i=phi.size();i--;)
            P[i] = phi[i];
        Integer prod = p;

        for(;--nbloops;) {
            do
                p = prevprime(MAXPRIMESIZE - (ulong)lrand48()%(ulong)MINLOWERSIZE).Int2long();
            while ((! isprime(p)) && (! isone(gcd(prod,p)))) ;
            UIndice rank;
            OneMinPolyAlg<BB>::operator()(phi,rank,p,Matrix_File_Name);
            for(Indice i=P.size();i--;)
                tcr(P[i],prod,phi[i],p);
            prod *= p;
        }
        
        Integer po2 = prod/2;
        for(Indice i=P.size();i--;) {
            if( P[i] > po2) P[i] = P[i] - prod;
            if (P[i] < -po2) P[i] = P[i] + prod;
        }
        
        _Comm.stop(LVL_NORMAL,PARTIAL_RESULT) << "Valence : " << P[0] << endl;
        return P;
    }

};


     


#endif __LINBOX_DOM_VALENCE_C__
