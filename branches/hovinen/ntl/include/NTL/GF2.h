
#ifndef NTL_GF2__H
#define NTL_GF2__H

#include <NTL/ZZ.h>

class GF2 {

long rep;

public:

GF2() { rep = 0; }

GF2(const GF2& a) : rep(a.rep) { }

~GF2() { }

GF2& operator=(const GF2& a) { rep = a.rep; return *this; }

GF2& operator=(long a) { rep = a & 1; return *this; }

friend void conv(GF2& x, long a) { x.rep = a & 1; }
friend GF2 to_GF2(long a) 
   { GF2 x; conv(x, a); return x; }

friend void conv(GF2& x, const ZZ& a) { x.rep = IsOdd(a); }
friend GF2 to_GF2(const ZZ& a) 
   { GF2 x; conv(x, a); return x; }

friend long rep(GF2 a) { return a.rep; }

static long modulus() { return 2; }
static GF2 zero() { return GF2(); }

friend void clear(GF2& x) { x.rep = 0; }
friend void set(GF2& x) { x.rep = 1; }

friend void swap(GF2& x, GF2& y)
   { long t; t = x.rep; x.rep = y.rep; y.rep = t; }

friend void add(GF2& x, GF2 a, GF2 b)
   { x.rep = a.rep ^ b.rep; }

friend void sub(GF2& x, GF2 a, GF2 b)
   { x.rep = a.rep ^ b.rep; }

friend void negate(GF2& x, GF2 a)
   { x.rep = a.rep; }

friend void add(GF2& x, GF2 a, long b)
   { x.rep = a.rep ^ (b & 1); }

friend void add(GF2& x, long a, GF2 b)
   { x.rep = (a & 1) ^ b.rep; }

friend void sub(GF2& x, GF2 a, long b)
   { x.rep = a.rep ^ (b & 1); }

friend void sub(GF2& x, long a, GF2 b)
   { x.rep = (a & 1) ^ b.rep; }

friend GF2 operator+(GF2 a, GF2 b)
    { GF2 x; add(x, a, b); return x; }

friend GF2 operator+(GF2 a, long b)
    { GF2 x; add(x, a, b); return x; }

friend GF2 operator+(long a, GF2 b)
    { GF2 x; add(x, a, b); return x; }

friend GF2 operator-(GF2 a, GF2 b)
    { GF2 x; sub(x, a, b); return x; }

friend GF2 operator-(GF2 a, long b)
    { GF2 x; sub(x, a, b); return x; }

friend GF2 operator-(long a, GF2 b)
    { GF2 x; sub(x, a, b); return x; }

friend GF2 operator-(GF2 a)
   { GF2 x; negate(x, a); return x; }

friend GF2& operator+=(GF2& x, GF2 b)
   { add(x, x, b); return x; }

friend GF2& operator+=(GF2& x, long b)
   { add(x, x, b); return x; }

friend GF2& operator-=(GF2& x, GF2 b)
   { sub(x, x, b); return x; }

friend GF2& operator-=(GF2& x, long b)
   { sub(x, x, b); return x; }

friend GF2& operator++(GF2& x) { add(x, x, 1); return x; }
friend void operator++(GF2& x, int) { add(x, x, 1); }
friend GF2& operator--(GF2& x) { sub(x, x, 1); return x; }
friend void operator--(GF2& x, int) { sub(x, x, 1); }


friend void mul(GF2& x, GF2 a, GF2 b)
   { x.rep = a.rep & b.rep; }

friend void mul(GF2& x, GF2 a, long b)
   { x.rep = a.rep & b; }

friend void mul(GF2& x, long a, GF2 b)
   { x.rep = a & b.rep; }

friend void sqr(GF2& x, GF2 a)
   { x = a; }

friend GF2 sqr(GF2 a)
   { return a; }

friend GF2 operator*(GF2 a, GF2 b)
    { GF2 x; mul(x, a, b); return x; }

friend GF2 operator*(GF2 a, long b)
    { GF2 x; mul(x, a, b); return x; }

friend GF2 operator*(long a, GF2 b)
    { GF2 x; mul(x, a, b); return x; }


friend GF2& operator*=(GF2& x, GF2 b)
   { mul(x, x, b); return x; }

friend GF2& operator*=(GF2& x, long b)
   { mul(x, x, b); return x; }


friend void div(GF2& x, GF2 a, GF2 b);
friend void div(GF2& x, long a, GF2 b);
friend void div(GF2& x, GF2 a, long b);

friend void inv(GF2& x, GF2 a);

friend GF2 inv(GF2 a)
   { GF2 x; inv(x, a); return x; }

friend GF2 operator/(GF2 a, GF2 b)
    { GF2 x; div(x, a, b); return x; }

friend GF2 operator/(GF2 a, long b)
    { GF2 x; div(x, a, b); return x; }

friend GF2 operator/(long a, GF2 b)
    { GF2 x; div(x, a, b); return x; }


friend GF2& operator/=(GF2& x, GF2 b)
   { div(x, x, b); return x; }

friend GF2& operator/=(GF2& x, long b)
   { div(x, x, b); return x; }


friend void power(GF2& x, GF2 a, long e);
friend GF2 power(GF2 a, long e)
   { GF2 x; power(x, a, e); return x; }


friend long IsZero(GF2 a)
   { return a.rep == 0; }

friend long IsOne(GF2 a)
   { return a.rep == 1; }

friend long operator==(GF2 a, GF2 b)
   { return a.rep == b.rep; }

friend long operator!=(GF2 a, GF2 b)
   { return !(a == b); }

friend long operator==(GF2 a, long b) { return a.rep == (b & 1); }
friend long operator==(long a, GF2 b) { return (a & 1) == b.rep; }

friend long operator!=(GF2 a, long b) { return !(a == b); }
friend long operator!=(long a, GF2 b) { return !(a == b); }

friend void random(GF2& x)
   { x.rep = RandomBnd(2); }

friend GF2 random_GF2()
   { GF2 x; random(x); return x; }

friend ostream& operator<<(ostream& s, GF2 a);

friend istream& operator>>(istream& s, GF2& x);

};


#endif

