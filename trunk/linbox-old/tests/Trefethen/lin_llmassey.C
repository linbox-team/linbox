// ======================================================================= //
// Linbox project 1999
// Domain Massey
// - Computation is stopped when the polynomials remain the same
//   for more than EARLY_TERM_THRESOLD
// - When minimal polynomial equals characteristic polynomial,
//   2 additional iterations are needed to compute it 
//   (parameter DEFAULT_ADDITIONAL_ITERATION), but those
//   iterations are not needed for the rank
// Time-stamp: <27 Aug 01 18:18:12 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef _LIN_DOM_MASSEY_C_
#define _LIN_DOM_MASSEY_C_

#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif

#include <LinBox/commentator.C>
#include <vector.h>

#define DEFAULT_EARLY_TERM_THRESHOLD 20
#define DEFAULT_ADDITIONAL_ITERATION 2

template<class Sequence> class MasseyDom {
public:
    typedef          Sequence                 Sequence_t;
    typedef typename Sequence::Domain_t       Domain_t;
    typedef typename Sequence::Type_t         Type_t;
    typedef 	     vector<Type_t>           PreferredPolynomial_t;
    typedef          MasseyDom< Sequence >    Self_t;
private:
    Sequence_t * _container;
    Domain_t _domain;
    Commentator _Comm;
    unsigned long long EARLY_TERM_THRESHOLD;

public:
        //-- Constructors
    MasseyDom(unsigned long long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
            : _container(), 
              _domain(), 
              _Comm(PRINT_NOTHING,PRINT_NOTHING),
              EARLY_TERM_THRESHOLD( ett_default )
        {}

    MasseyDom(const Self_t& M, unsigned long long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
            : _container(M._container), 
              _domain(M._domain), 
              _Comm(M._Comm) ,
              EARLY_TERM_THRESHOLD( ett_default )
        {}

    MasseyDom(Sequence_t * D, unsigned long long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
            : _container(D), 
              _domain(D->getdomain()), 
              _Comm(PRINT_NOTHING,PRINT_NOTHING),
              EARLY_TERM_THRESHOLD( ett_default )
        {}
  
    MasseyDom(Sequence_t * MD, const Domain_t& D, unsigned long long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
            : _container(MD), 
              _domain(D), 
              _Comm(PRINT_NOTHING,PRINT_NOTHING),
              EARLY_TERM_THRESHOLD( ett_default ) 
        {}  

       //-- with Commentator
    MasseyDom(const Commentator& C, unsigned long long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
            : _container(), 
              _domain(), 
              _Comm(C) ,
              EARLY_TERM_THRESHOLD( ett_default )
        {}

    MasseyDom(const Commentator& C, Sequence_t * D, unsigned long long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
            : _container(D), 
              _domain(D->getdomain()), 
              _Comm(C) ,
              EARLY_TERM_THRESHOLD( ett_default )
        {}
  
    MasseyDom(const Commentator& C, Sequence_t * MD, const Domain_t& D, unsigned long long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
            : _container(MD), 
              _domain(D), 
              _Comm(C) ,
              EARLY_TERM_THRESHOLD( ett_default )
        {}  

        //-- Principal method
    template<class Polynomial>
    void operator() (Polynomial& C, bool full_poly = 0) {
        massey(C, full_poly);
    };
    
        //-- Domains access
    const Domain_t& getdomain() const { return _domain; }
    Sequence_t * getsequence() const { return _container; }
    

private:
// -----------------------------------------------
// Polynomial emulation
// Only container aspects of polynomials
// AND degree and valuation are needed !
// -----------------------------------------------
#define _DEGINFTY_ -1

// Degree of v
template <class V>
long long v_degree(V& v) {
   long long i = v.size()-1;
   if (i == _DEGINFTY_) return _DEGINFTY_ ;
   else if (!_domain.iszero(v[i]) ) return i ;
   // We must re-compute the degree :
   for (long long j=i-1; j>=0 ; j--)
      if (!_domain.iszero((v)[j]) ) {
          v.resize(j+1);
          return j ;
      }
   return _DEGINFTY_ ;
}

// Valuation of v
template <class V>
long long v_val(V& v) {
   long long i = v.size()-1;
   if (i == _DEGINFTY_) return _DEGINFTY_ ;
   else if (!_domain.iszero(v[0]) ) return 0;
   // We must compute the valuation :
   for (long long j=1; j<=i ; j++)
      if (!_domain.iszero((v)[j]) ) return j ;
   return _DEGINFTY_ ;
}



// -------------------------------------------------------------------
// Berlekamp/Massey algorithm with Massey's Sequence generation
// -------------------------------------------------------------------

    template<class Polynomial>
    void massey(Polynomial& C, bool full_poly = 0) { 
//         const long long ni = _container->n_row(), nj = _container->n_col();
//         const long long n = GIVMIN(ni,nj);
        const long long END = _container->size() + (full_poly? DEFAULT_ADDITIONAL_ITERATION:0);
        const long long n = END >> 1;

        _Comm.start("Massey",LVL_NORMAL,INTERNAL_DESCRIPTION) 
            << END << endl;

// ====================================================
// Sequence and iterator initialization
//
        typename Sequence_t::const_iterator _iter( _container->begin() );
        Polynomial S(END+1);

// -----------------------------------------------
// Preallocation. No further allocation.
//
        C.reserve(n+1); C.resize(1); C[0] = _domain.one; 
        Polynomial B(n+1); B.resize(1); B[0] = _domain.one;

        long long L=0;
        Type_t b = _domain.one, d, Ds;
        long long x = 1, b_deg = 0, c_deg = 0, l_deg;
	long long COMMOD = (END>25?END/25:1);

        for (long long N=0; N<END && x<EARLY_TERM_THRESHOLD; ++N, ++_iter) {
            if ( ! (N % COMMOD) ) 
                _Comm.progress("m-v prods",LVL_IMP,N,END);
// ====================================================
// Next coefficient in the sequence
// Discrepancy computation
// 
            d = S[N] = *_iter; 
//             d = _container->next( S[N] );
            for(long long i=GIVMIN(L,c_deg);i;--i) _domain.axpyin(d,C[i],S[N-i]);
            if (_domain.iszero(d)) {
                ++x;
            } else {
                if (L > (N>>1)) {
// -----------------------------------------------
// C = C + (Polynome(X,x,-d/b) * B);
// 
                    _domain.divin( _domain.neg(Ds,d), b);
                    long long i = l_deg = (x+b_deg);
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
                    long long i = l_deg = x+b_deg;
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


#ifdef __DEBUG__
        for(unsigned long long i = 0; i < S.size(); ++i) 
            cerr << S[i] << "*X^" << i << " + ";
        cerr << endl;
#endif
        
        _Comm.stop(LVL_NORMAL,PARTIAL_RESULT) 
            << "Degree : " << v_degree(C)-v_val(C) 
            << " over GF(" << _domain.size() << "), 0:" << x << endl;
    }

public:
// ---------------------------------------------
// Massey
//
    void pseudo_rank(unsigned long long& rank) {
        PreferredPolynomial_t phi;
        massey(phi, 0);
        rank = v_degree(phi) - v_val(phi);
    };
 
    void valence(Type_t& valence, unsigned long long& rank) {
        _Comm.start("Valence",LVL_NORMAL,INTERNAL_DESCRIPTION) << endl;
        PreferredPolynomial_t phi;
        massey(phi, 1);
        rank = v_degree(phi) - v_val(phi);
        valence = phi[v_degree(phi)] ;

        _domain.write( _Comm.stop(LVL_NORMAL,PARTIAL_RESULT) , valence)
            << ", rank : " << rank
            << endl;
    }

    template<class Polynomial>
    void pseudo_minpoly(Polynomial& phi, unsigned long long& rank, bool full_poly = 1) {
        massey(phi, full_poly);
        rank = v_degree(phi) - v_val(phi);
        if (phi.size()) {
            long long dp = v_degree(phi);
            for(long long i = dp >> 1;i > 0; --i) {
                phi[0] = phi[i];
                phi[i] = phi[dp-i];
                phi[dp-i] = phi[0];
            }
            phi[0] = phi[dp];
            phi[dp] = _domain.one;
        }
    }

};


    
#endif
