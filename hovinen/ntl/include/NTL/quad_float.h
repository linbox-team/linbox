
#ifndef NTL_quad_float__H
#define NTL_quad_float__H

#include <iostream.h>
#include <math.h>
#include <NTL/ZZ.h>


class quad_float {
public:
  double hi, lo;

  // Constructors
  quad_float() : hi(0), lo(0)  {}

  inline quad_float& operator=(double x);

  static long oprec;
  static void SetOutputPrecision(long p);
  static long OutputPrecision() { return oprec; }

  quad_float(double x, double y) : hi(x), lo(y) { } // internal use only

  ~quad_float() {}

};  // end class quad_float

#if (NTL_BITS_PER_INT < NTL_DOUBLE_PRECISION)
inline quad_float to_quad_float(int n) { return quad_float(n, 0); }
#else
quad_float to_quad_float(int n);
#endif

#if (NTL_BITS_PER_LONG < NTL_DOUBLE_PRECISION)
inline quad_float to_quad_float(long n) { return quad_float(n, 0); }
#else
quad_float to_quad_float(long n);
#endif


inline quad_float to_quad_float(double x) { return quad_float(x, 0); }
// On platforms with extended doubles, this may result in an
// improper quad_float object, but it should be converted to a proper
// one when passed by reference to any of the arithmetic routines,
// at which time x will be forced to memory.

inline quad_float to_quad_float(float x) 
   { return to_quad_float(double(x)); }

inline quad_float& quad_float::operator=(double x) 
   { *this = to_quad_float(x); return *this; }

quad_float operator+(const quad_float&, const quad_float& );

inline quad_float operator+(const quad_float& x, double y )
   { return x + to_quad_float(y); }

inline quad_float operator+(double x, const quad_float& y)
   { return to_quad_float(x) + y; }

quad_float operator-(const quad_float&, const quad_float& );

inline quad_float operator-(const quad_float& x, double y )
   { return x - to_quad_float(y); }

inline quad_float operator-(double x, const quad_float& y)
   { return to_quad_float(x) - y; }

quad_float operator*(const quad_float&, const quad_float& );

inline quad_float operator*(const quad_float& x, double y )
   { return x * to_quad_float(y); }

inline quad_float operator*(double x, const quad_float& y)
   { return to_quad_float(x) * y; }

quad_float operator/(const quad_float&, const quad_float& );

inline quad_float operator/(const quad_float& x, double y )
   { return x / to_quad_float(y); }

inline quad_float operator/(double x, const quad_float& y)
   { return to_quad_float(x) / y; }

quad_float operator-(const quad_float& x);

quad_float& operator+= (quad_float& x, const quad_float& y);
inline quad_float& operator += (quad_float& x, double y)
   { x += to_quad_float(y); return x; }

quad_float& operator-= (quad_float& x, const quad_float& y);
inline quad_float& operator-= (quad_float& x, double y)
   { x -= to_quad_float(y); return x; }

quad_float& operator*= (quad_float& x, const quad_float& y);
inline quad_float& operator*= (quad_float& x, double y)
   { x *= to_quad_float(y); return x; }

quad_float& operator/= (quad_float& x, const quad_float& y);
inline quad_float& operator/= (quad_float& x, double y)
   { x /= to_quad_float(y); return x; }

inline quad_float& operator++(quad_float& a) { a += 1.0; return a; }
inline void operator++(quad_float& a, int) { a += 1.0; }

inline quad_float& operator--(quad_float& a) { a -= 1.0; return a; }
inline void operator--(quad_float& a, int) { a -= 1.0; }


long operator> (const quad_float& x, const quad_float& y);
long operator>=(const quad_float& x, const quad_float& y);
long operator< (const quad_float& x, const quad_float& y);
long operator<=(const quad_float& x, const quad_float& y);
long operator==(const quad_float& x, const quad_float& y);
long operator!=(const quad_float& x, const quad_float& y);

inline long operator> (const quad_float& x, double y) 
   { return x > to_quad_float(y); }
inline long operator> (double x, const quad_float& y)
   { return to_quad_float(x) > y; }

inline long operator>=(const quad_float& x, double y) 
   { return x >= to_quad_float(y); }
inline long operator>=(double x, const quad_float& y)
   { return to_quad_float(x) >= y; }

inline long operator< (const quad_float& x, double y) 
   { return x < to_quad_float(y); }
inline long operator< (double x, const quad_float& y)
   { return to_quad_float(x) < y; }

inline long operator<=(const quad_float& x, double y) 
   { return x <= to_quad_float(y); }
inline long operator<=(double x, const quad_float& y)
   { return to_quad_float(x) <= y; }

inline long operator!=(const quad_float& x, double y) 
   { return x != to_quad_float(y); }
inline long operator!=(double x, const quad_float& y)
   { return to_quad_float(x) != y; }

inline long operator==(const quad_float& x, double y) 
   { return x == to_quad_float(y); }
inline long operator==(double x, const quad_float& y)
   { return to_quad_float(x) == y; }


inline long sign(const quad_float& x){
  if (x.hi>0.0) return 1; else if (x.hi<0.0) return -1; else return 0;
}

long compare(const quad_float&, const quad_float&);

inline long compare(const quad_float& x, double y)
   { return compare(x, to_quad_float(y)); }

inline long compare(double x, const quad_float& y)
   { return compare(to_quad_float(x), y); }



istream& operator >> (istream&, quad_float&);
ostream& operator << (ostream&, const quad_float&);


quad_float sqrt(const quad_float&);
quad_float floor(const quad_float&);
quad_float ceil(const quad_float&);
quad_float trunc(const quad_float&);
quad_float fabs(const quad_float&);

void power(quad_float&, const quad_float&, long);
inline quad_float power(const quad_float& x, long e)
   { quad_float z; power(z, x, e); return z; }

void power2(quad_float&, long);
inline quad_float power2_quad_float(long e)
   { quad_float z; power2(z, e); return z; }


long to_long(const quad_float&);
inline int to_int(const quad_float& x) { return int(to_long(x)); }

inline double to_double(const quad_float& x) { return x.hi; }

inline float to_float(const quad_float& x) { return float(x.hi); }


inline void conv(quad_float& x, int a) { x = to_quad_float(a); }
inline void conv(quad_float& x, long a) { x = to_quad_float(a); }
inline void conv(quad_float& x, float a) { x = to_quad_float(a); }
inline void conv(quad_float& x, double a) { x = to_quad_float(a); }


inline void conv(long& x, const quad_float& a) { x = to_long(a); }
inline void conv(int& x, const quad_float& a) { x = to_int(a); }
inline void conv(double& x, const quad_float& a) { x = to_double(a); }
inline void conv(float& x, const quad_float& a) { x = to_float(a); }

void conv(quad_float&, const ZZ&);
inline quad_float to_quad_float(const ZZ& x)
   { quad_float z; conv(z, x); return z; }

void conv(ZZ&, const quad_float&);
inline ZZ to_ZZ(const quad_float& a)
   { ZZ x; conv(x, a); NTL_OPT_RETURN(ZZ, x); }

inline void conv(quad_float& x, const quad_float& a)
   { x = a; }
inline quad_float to_quad_float(const quad_float& a)
   { return a; }

quad_float to_quad_float(const char *s);
inline void conv(quad_float& x, const char *s)
   { x = to_quad_float(s); }

long IsFinite(quad_float *x);

long PrecisionOK();

quad_float ldexp(const quad_float& x, long exp);

quad_float exp(const quad_float& x);
quad_float log(const quad_float& x);


#endif
