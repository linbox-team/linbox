// =================================================================== //
// Givaro : Prime numbers
//              Factors,
// Needs list structures : stl ones for instance
// Time-stamp: <30 Aug 00 20:26:22 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef _GIVARO_FACTORISATION_INL_
#define _GIVARO_FACTORISATION_INL_

#include "givtimer.h"
#include "givinteger.h"
#include "givintprime.h"
#include "givintfactor.h"

template<class RandIter>
ostream& IntFactorDom<RandIter>::write(ostream& o, const Rep& n) const
{
//         // n = * Lf[i] ^ Lo[i]
//         // But Lf[i] might not be prime (cf. factor probability)
    Rep nn,g,r,u;
    nn = n;
    long c;
    int flag = 0;
    Integer tmp;
    
    while( isgt(nn,1) ) {
        factor(g,nn);
	if (flag)
		o << " * ";
	else
		flag = 1;
	IntegerDom::write(o,g);
        c=0;r=0;

	Rep::divexact(u, nn,g);
        while(iszero(r)) {
		nn.copy(u); 
		Rep::divmod( u, r, nn,g );
             	c++;
        }
	if (c>1) o << "^" << c;
    }
    return o;
}    

// =================================================================== //
// Set or Container of divisors, factors.
// =================================================================== //
template<class RandIter>
template< template<class> class Container> unsigned long 
IntFactorDom<RandIter>::set(Container<Rep>& Lf, Container<unsigned long>& Lo, const Rep& n)  const 
{
        // n = * Lf[i] ^ Lo[i]
        // But Lf[i] might not be prime (cf. factor probability)
    
    Rep nn,g,r,u;
    if (n<0) Rep::neg(nn,n); else nn=n;
    unsigned long c,nb=0;
    while(nn > 1) {        
        factor(g,nn);
        Lf.push_back(g);
        c=0;r=zero;
	Rep::divexact(u, nn,g);
        while(r == 0) {
        //	nn = nn / g;
        //	r = nn % g;
		nn.copy(u); 
		Rep::divmod( u, r, nn,g );
            	c++;
        }
        Lo.push_back( c );
        nb++;
    }
    return nb;
}

template<class RandIter>
template< template<class> class Container> 
void IntFactorDom<RandIter>::set(Container<Rep>& Lf, Container<long>& Lo, const Rep& n)  const 
{
        // n = * Lf[i] ^ Lo[i]
        // But Lf[i] might not be prime (cf. Pollard probability)
    Rep nn,g,r,u;
    nn = n;
    long c;
    while(nn > 1) {        
        factor(g,nn);
        Lf.push_back(g);
        c=0;r=0;
	Rep::divexact(u, nn,g);
        while(r == 0) {
        //	nn = nn / g;
        //	r = nn % g;
		nn.copy(u); 
		Rep::divmod( u, r, nn,g ); // divide(nn,g,u,r);
            c++;
        }
        Lo.push_back( c );
    }
}           

template<class RandIter>
template< template<class> class Container> 
void IntFactorDom<RandIter>::Erathostene(Container<Rep>& Lf, const Rep& p)  const 
{
        // Deterministic algorithm
        // Valid for p < BOUNDARY_factor
        // Lf is the Container of factors
    long n; access(n,p);
    if (! (n & 0x1)) {
        Lf.push_back(Rep(2));
        do
            n >>= 1;
        while (!(n & 0x1));
    }
    short * IP = new short[n+1];
    int i;
    for(int i=n+1;i--;)
	IP[i] = 0L;
    i=3;
    int j, ii, cofact;
    Rep sq;
    while (i<=sqrt(sq,Rep(n))) {
        ii= i << 1;
        j = i+ii;
        while (j<=n) {
            IP[j] = 1L;
            j+=ii;
        }
        if ((j-ii) == n) {
            Lf.push_back(Rep(i));
            do
                n /= i;
            while (!(n%i));
        }
        j = i+1;
        while (IP[++j]) { j++;}
        i = j;
    }
    if (!(IP[n]) && (n>1)) Lf.push_back(Rep(n));
    delete [] IP;
}
 

template<class RandIter>
template< template<class> class Container> 
void IntFactorDom<RandIter>::set( Container<Rep>& Lf,  const Rep& n)  const 
{
        // big_factor is executed until
        // a (sometimes probably) prime factor is found.
        // Lf is a Container of divisors
        // Lf is a Container of factors with probability -- TO EXPLICIT
        // something like : 1 - (big_factor to be composite)*(big_isprime)
    Rep nn,g,r,u;
    nn = n;
    while(nn > 1) {
        factor(g,nn);
        r=0;
	Rep::divexact(u, nn,g);
        while(r == 0) {
		nn.copy(u); 
		Rep::divmod( u, r, nn,g ); 
        }
//         if (isprime(g))
            Lf.push_back(g);
//         else
//             IntFactorDom::set(Lf, g);
    }
}


template<class RandIter>
template< template<class> class Container> Container< IntFactorDom<RandIter>::Rep >&  IntFactorDom<RandIter>::divisors( Container<Rep>& L, const Container<Rep>& Lf, const Container<unsigned long>& Le)  const 
{
    typename Container<Rep>::const_iterator li = Lf.begin();
    typename Container<unsigned long>::const_iterator lj = Le.begin();
    Container<Rep> Res(1,Rep(1));
    Container<Rep> Res2;
    typename Container<Rep>::iterator lr;
    Rep Itmp;
    for(;li!=Lf.end();++li,++lj) {
        for(lr = Res.begin();lr!=Res.end();++lr) {
		Itmp = *lr;
      		for(unsigned long i=*lj;i--;) {
			Itmp = Itmp * *li;
            		Res2.push_back( Itmp );
		}
	}
        Res.splice(Res.end(),Res2);
    }
    return L = Res;
}


template<class RandIter>
template< template<class> class Container> Container<IntFactorDom<RandIter>::Rep>& IntFactorDom<RandIter>::divisors( Container<Rep>& L, const Rep& n)  const 
{
    Container<Rep> Lf;
    Container<unsigned long> Le;
    IntFactorDom<RandIter>::set(Lf,Le,n);
    return divisors(L, Lf, Le);
}

template<class RandIter>
IntFactorDom<RandIter>::Rep& IntFactorDom<RandIter>::Pollard(RandIter& gen, Rep& g, const Rep& n, const unsigned long threshold = 0) const 
{
  // average number of iterations < 13/8*sqrt( Pi*n/2)
  // Sometimes the factor isn't prime -- TO EXPLICIT
    if (islt(n,3)) return g=n;
    if ( isprime(n) ) return g=n;
    g=1;
    Rep m(zero), x, y, p(one), t;
    random(gen, y, n);

    if (threshold) {
        unsigned long c = 0;
        while(isone(g) && (++c < threshold)) {
            if(  isequal(p, addin(m,one)) ) {
                x=y;
                mulin(p,2);
            }
            Pollard_fctin(y,n);
            gcd(g,sub(t,y,x),n);
        }
    } else {
        while(isone(g)) {
            if(  isequal(p, addin(m,one)) ) {
                x=y;
                mulin(p,2);
            }
            Pollard_fctin(y,n);
            gcd(g,sub(t,y,x),n);
        }
    }
    return g;
}


// =================================================================== //
// Elliptic curves routines
// =================================================================== //

inline void Add_Curve( const Integer& n, const Integer A, const Integer& ax, const Integer& az, Integer& cx, Integer& cz) {
    Integer t1,t2,t3,t4;
    t3 = ax+az; t4 = ax-az;
//     t1 = ((ax+az)*(ax+az))%n; 
//     t2 = ((ax-az)*(ax-az))%n;
    t1 = (t3*t3)%n; 
    t2 = (t4*t4)%n;
    cx = (t1*t2)%n; 
    cz = ((t1-t2)*((A*(t1-t2)+t2)%n))%n;
}

inline void one_Mul_Curve( const Integer& n, const Integer A, const Integer& mm, const Integer& nn, const Integer& px, const Integer& pz, Integer& aax, Integer& aaz) {
    Integer ax, az, bx, bz, cx, cz, tmpx, tmpz, d, e, t1, t2,t3,t4;
    cx = px;
    cz = pz;
    e = mm;
    d = nn - mm;
    if (e<d) {
        Add_Curve(n,A,px,pz,bx,bz);ax = px;az = pz;d = d - e;
    } else {
        Add_Curve(n,A,px,pz,ax,az);bx = px;bz = pz;e = e - d;
    }
    while (! iszero(e)) {
        if (e<d) {
            tmpx = bx; tmpz = bz; 
            t1 = ((ax-az)*(bx+bz))%n; 
            t2 = ((ax+az)*(bx-bz))%n; 
//             bx = (cz*(((t1+t2)*(t1+t2))%n))%n;
//             bz = (cx*(((t1-t2)*(t1-t2))%n))%n;
            t3 = t1+t2;
            t4 = t1-t2;
            bx = (cz*(((t3)*(t3))%n))%n;
            bz = (cx*(((t4)*(t4))%n))%n;
            d = d-e;
        } else {
            tmpx = ax; tmpz = az; 
            t1 = ((ax-az)*(bx+bz))%n; 
            t2 = ((ax+az)*(bx-bz))%n; 
            t3 = t1+t2;
            t4 = t1-t2;
//             ax = (cz*(((t1+t2)*(t1+t2))%n))%n;
//             az = (cx*(((t1-t2)*(t1-t2))%n))%n;
            ax = (cz*(((t3)*(t3))%n))%n;
            az = (cx*(((t4)*(t4))%n))%n;
            e = e-d;
        }
        cx = tmpx;
        cz = tmpz;
    }
    aax = ax;
    aaz = az;
}        
    

inline void Mul_Curve( const Integer& n, Integer& Ai, const Integer& mm, const Integer& nn, const Integer& B1, Integer& Xi, Integer& Zi) {
    Integer pow = nn;
    while (pow <= B1) {
        one_Mul_Curve(n,Ai,mm,nn,Xi,Zi,Xi,Zi);
        pow *= nn;
    }
}

// ======================================================================== //
// Lenstra algorithm for elliptic curves
// Returns -1 if failure to find factors
// heuristically exp( sqrt( (2+epsilon)(ln p)(ln ln p) ) ) multiplications
// to find a factor p of N.
// TODO : make it generic in regards to DOMAINLIKENESS
// ======================================================================== //
template<class RandIter>
IntFactorDom<RandIter>::Rep& IntFactorDom<RandIter>::Lenstra(RandIter& gen, Rep& g, const Rep& n, const Rep& B1 = 1000000, const unsigned long curves = 30) const 
{
    if (n<3) return g=n;
    if ( isprime(n) ) return g=n;
    if (iszero(n % 2)) return g=2;
    if (iszero(n % 3)) return g=3;

    Rep * A = new Rep[curves+1]
        , * X = new Rep[curves+1]
        , * Z = new Rep[curves+1];

    Rep r,a,asq,kg,kgg;
    for (long c=curves+1;c--;)
        Z[c] = one;
    
    Rep u,v, four;
    assign(four,4);
//     g = gcd(4,n,u,v);
    gcd(g,u,v,four,n);
    Rep inv4 = u, sixt;
    assign(sixt,16);
//     g = gcd(16,n,u,v);
    gcd(g,u,v,sixt,n);
    Rep inv16 = u,inva;

// Initialize # curves
    for (long i=1;i<=curves;++i) {
        a = 0, asq = 0;
        while ((( a*(asq-1)*(9*asq-1) ) % n) == 0 ) {
            random(gen,r,n);
//             kg = r*r + 6;
            mul(kg,r,r); addin(kg,6);
//             kgg = gcd(kg,n);
            gcd(kgg,kg,n);
            if (isone(kgg)) {
//                 g = gcd(kg,n,u,v); if (! isone(g)) { delete [] A, X, Z; return g; }
                gcd(g,u,v,kg,n); if (! isone(g)) { delete [] A, X, Z; return g; }
                a = (6*r*u)% n;
                asq = (a * a) % n ;
            } else
                return g=kgg;
        }
//         g = gcd(a,n,u,v); if (! isone(g)) { delete [] A, X, Z; return g; }
        gcd(g,u,v,a,n); if (! isone(g)) { delete [] A, X, Z; return g; }
        A[i] = ( (8-3*a + (1-6*a*a)*u*u*u )*inv16 ) % n;
        X[i] = (3*a*inv4)%n;
    }

        // .5*sqrt(5)-.5, 37 digits
//     Rep s("618033988749894848204586834370");
//     Rep si("1000000000000000000000000000000");
    Rep s("6180339887498948482045868343656381177");
    Rep si("10000000000000000000000000000000000000");

// Begins search with curves on primes up to B1
    Rep prime = 2, sp, f;
    while (prime <= B1) {
//         cerr << "p: " << prime << endl;
        sp = (prime*s)/si;
        Mul_Curve(n,A[1],sp,prime,B1,X[1],Z[1]);
        f = Z[1];
        for(long i=2;i<=curves;++i) {
            Mul_Curve(n,A[i],sp,prime,B1,X[i],Z[i]);
            f = (f*Z[i])%n;
        }
//         f = gcd(f,n);
        Rep ftm;
        gcd(ftm,f,n);
        f = ftm;
        if (isone(f)) {
            nextprime(ftm,prime);
            prime = ftm;
//             prime = nextprime(prime);
        } else {
            delete [] A, X, Z; 
            return g=f; 
        }
    }
    
    cerr << "*** Elliptic curves with " << curves << " curves, threshold " << B1 << " failed ***" << endl;
    delete [] A, X, Z; 
    return neg(g,one);
}

#endif
