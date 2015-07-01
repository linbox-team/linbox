// ======================================================================= //
// Linbox project 1999
// Domain Massey
// Multiply computations are stopped when the polynomials remains the same
// for more than EARLY_TERM_THRESOLD
// Time-stamp: <20 Dec 99 16:24:16 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef _LIN_DOM_MASSEY_C_
#define _LIN_DOM_MASSEY_C_
#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif

#include "commentator.C"
#include <vector.h>

#define EARLY_TERM_THRESHOLD 10

template<class ThunkDomain> class MasseyDom {
public:
    typedef          ThunkDomain                 ThunkDomain_t;
    typedef typename ThunkDomain::Domain_t       Domain_t;
    typedef typename ThunkDomain::Type_t         Type_t;
    typedef 	     vector<Type_t>              PreferredPolynomial_t;
    typedef          MasseyDom< ThunkDomain >    Self_t;
private:
    ThunkDomain_t * _BB_domain;
    Domain_t _domain;
    Commentator _Comm;

public:
        //-- Constructors
    MasseyDom() 
            : _BB_domain(), 
              _domain(), 
              _Comm(PRINT_NOTHING,PRINT_NOTHING) 
        {}

    MasseyDom(const Self_t& M) 
            : _BB_domain(M._BB_domain), 
              _domain(M._domain), 
              _Comm(M._Comm) 
        {}

    MasseyDom(ThunkDomain_t * D) 
            : _BB_domain(D), 
              _domain(D->getdomain()), 
              _Comm(PRINT_NOTHING,PRINT_NOTHING) 
        {}
  
    MasseyDom(ThunkDomain_t * MD, const Domain_t& D) 
            : _BB_domain(MD), 
              _domain(D), 
              _Comm(PRINT_NOTHING,PRINT_NOTHING) 
        {}  

       //-- with Commentator
    MasseyDom(const Commentator& C) 
            : _BB_domain(), 
              _domain(), 
              _Comm(C) 
        {}

    MasseyDom(const Commentator& C, ThunkDomain_t * D) 
            : _BB_domain(D), 
              _domain(D->getdomain()), 
              _Comm(C) 
        {}
  
    MasseyDom(const Commentator& C, ThunkDomain_t * MD, const Domain_t& D) 
            : _BB_domain(MD), 
              _domain(D), 
              _Comm(C) 
        {}  

        //-- Principal method
    template<class Polynomial, class InMatrix>
    void operator() (Polynomial& C) {
        massey(C);
    };
    
        //-- Domains access
    Domain_t getdomain() const { return _domain; }
    ThunkDomain_t * getBBdomain() const { return _BB_domain; }
    

private:
// -----------------------------------------------
// Polynomial emulation
// Only container aspects of polynomials
// AND degree and valuation are needed !
// -----------------------------------------------
#define _DEGINFTY_ -1

// Degree of v
template <class V>
long v_degree(V& v) {
   long i = v.size()-1;
   if (i == _DEGINFTY_) return _DEGINFTY_ ;
   else if (!_domain.iszero(v[i]) ) return i ;
   // We must re-compute the degree :
   for (long j=i-1; j>=0 ; j--)
      if (!_domain.iszero((v)[j]) ) {
          v.resize(j+1);
          return j ;
      }
   return _DEGINFTY_ ;
}

// Valuation of v
template <class V>
long v_val(V& v) {
   long i = v.size()-1;
   if (i == _DEGINFTY_) return _DEGINFTY_ ;
   else if (!_domain.iszero(v[0]) ) return 0;
   // We must compute the valuation :
   for (long j=1; j<=i ; j++)
      if (!_domain.iszero((v)[j]) ) return j ;
   return _DEGINFTY_ ;
}



// -------------------------------------------------------------------
// Berlekamp/Massey algorithm with Massey's Sequence generation
// -------------------------------------------------------------------

    template<class Polynomial>
    void massey(Polynomial& C) { 
        const long ni = _BB_domain->n_row(), nj = _BB_domain->n_col();
        const long n = GIVMIN(ni,nj);
        const long END = 2*n;

        _Comm.start("Massey",LVL_NORMAL,INTERNAL_DESCRIPTION) 
            << ni << " x " << nj << endl;

// ====================================================
// Sequence initialization
//
        Polynomial S(END+1);

// -----------------------------------------------
// Preallocation. No further allocation.
//
        C.reserve(n+1); C.resize(1); C[0] = _domain.one; 
        Polynomial B(n+1); B.resize(1); B[0] = _domain.one;

        long L=0;
        Type_t b = _domain.one, d, Ds;
        long x = 1, b_deg = 0, c_deg = 0, l_deg;

        for (long N=0; N<END && x<EARLY_TERM_THRESHOLD; ++N) {
            if ( ! (N % 1000) ) 
                _Comm.progress("matrix-vector products",LVL_IMP,N,END);
// ====================================================
// Next coefficient in the sequence
// Discrepancy computation
// 
            d = _BB_domain->next( S[N] );
            for(long i=GIVMIN(L,c_deg);i;--i) _domain.axpyin(d,C[i],S[N-i]);
            if (_domain.iszero(d)) {
                ++x;
            } else {
                if (L > (N>>1)) {
// -----------------------------------------------
// C = C + (Polynome(X,x,-d/b) * B);
// 
                    _domain.divin( _domain.neg(Ds,d), b);
                    long i = l_deg = x+b_deg;
                    if (l_deg > c_deg) {
                        C.resize(l_deg+1);
                        if (x>c_deg) {
                            for(;i>=x;--i) _domain.mul(C[i], Ds, B[i-x]);
                            for(;i>c_deg;--i) C[i] = _domain.zero;
                        } else {
                            for(;i>c_deg;--i) _domain.mul(C[i], Ds, B[i-x]);
                            for(;i>=x;--i) _domain.axpyin(C[i], Ds, B[i-x]) ;
                        }
                    } else {
                        for(;i>=x;--i) _domain.axpyin(C[i],Ds,B[i-x]);
                    }
// -----------------------------------------------
                    c_deg = v_degree(C);
                    ++x;
                } else {
// -----------------------------------------------
// C = C + (Polynome(X,x,-d/b) * B); // B = C;
// 
                    _domain.divin( _domain.neg(Ds,d), b);
                    long i = l_deg = x+b_deg;
                    B.resize(C.size());
                    if (l_deg > c_deg) {
                        C.resize(l_deg+1);
                        if (x>c_deg) {
                            for(;i>=x;--i) _domain.mul(C[i], Ds, B[i-x]);
                            for(;i>c_deg;--i) C[i] = _domain.zero;
                        } else {
                            for(;i>c_deg;--i) _domain.mul(C[i], Ds, B[i-x]);
                            for(;i>=x;--i) _domain.axpy(C[i],Ds,B[i-x], B[i] = C[i] ) ;
                        }
                    } else {
                        for(i=c_deg;i>l_deg;--i) B[i] = C[i];
                        for(;i>=x;--i) _domain.axpy(C[i],Ds,B[i-x], B[i] = C[i] );
                    }
                    for(;i>=0;--i) B[i] = C[i];
// -----------------------------------------------
                    L = N+1-L;
                    b_deg = c_deg;
                    c_deg = v_degree( C );
                    b = d;
                    x = 1;
                }
            }
// ====================================================
        }
        _Comm.stop(LVL_NORMAL,INTERNAL_DESCRIPTION) 
            << "Rank : " << v_degree(C)-v_val(C) 
            << " over GF(" << _domain.size() << ")" << endl;
    }

public:
// ---------------------------------------------
// Massey
//
    void pseudo_rank(unsigned long& rank) {
        PreferredPolynomial_t phi;
        massey(phi);
        rank = v_degree(phi) - v_val(phi);
    };
 
    void valence(Type_t& valence, ulong& rank) {
        _Comm.start("Valence",LVL_NORMAL,INTERNAL_DESCRIPTION) << endl;
        PreferredPolynomial_t phi;
        massey(phi);
        rank = v_degree(phi) - v_val(phi);
        valence = phi[v_degree(phi)] ;
        _Comm.stop(LVL_NORMAL,PARTIAL_RESULT) 
            << _domain.access(valence) 
            << ", rank : " << rank
            << endl;
    }

    template<class Polynomial>
    void pseudo_minpoly(Polynomial& phi, ulong& rank) {
//         _Comm.start("MinPoly",LVL_NORMAL,INTERNAL_DESCRIPTION) << endl;
        massey(phi);
        rank = v_degree(phi) - v_val(phi);
        if (phi.size()) {
            long dp = v_degree(phi);
            for(long i = dp >> 1;i > 0; --i) {
                phi[0] = phi[i];
                phi[i] = phi[dp-i];
                phi[dp-i] = phi[0];
            }
            phi[0] = phi[dp];
            phi[dp] = _domain.one;
        }
//         _Comm.stop(LVL_NORMAL,PARTIAL_RESULT) 
//             << _domain.access(phi[0]) 
//             << ", rank : " << v_degree(phi) - v_val(phi)
//             << endl;
    }

};


    
#endif _LIN_DOM_MASSEY_C_
