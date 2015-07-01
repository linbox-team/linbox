
#ifndef NTL_mat_GF2__H
#define NTL_mat_GF2__H


#include <NTL/vec_vec_GF2.h>


class mat_GF2 {  
private:  
  
   vec_vec_GF2 rep;  
   long numcols;  
  
public:  
  
   mat_GF2() { numcols = 0; }  
   mat_GF2(const mat_GF2& a);  
   mat_GF2& operator=(const mat_GF2& a);  
   ~mat_GF2() { }  
  
   mat_GF2(INIT_SIZE_TYPE, long n, long m);  
  
   void kill();  
  
   void SetDims(long n, long m);  
  
   long NumRows() const { return rep.length(); }  
   long NumCols() const { return numcols; }  
  
   vec_GF2& operator[](long i) { return rep[i]; }  
   const vec_GF2& operator[](long i) const { return rep[i]; }  
  
   vec_GF2& operator()(long i) { return rep[i-1]; }  
   const vec_GF2& operator()(long i) const { return rep[i-1]; }  

   GF2 get(long i, long j) const { return rep[i].get(j); }
   void put(long i, long j, GF2 a) { rep[i].put(j, a); }
   void put(long i, long j, long a) { rep[i].put(j, a); }

   subscript_GF2 operator()(long i, long j)
      { return subscript_GF2(rep[i-1], j-1); }

   const_subscript_GF2 operator()(long i, long j) const
      { return const_subscript_GF2(rep[i-1], j-1); }

   long position(const vec_GF2& a) const { return rep.position(a); }
  
   friend const vec_vec_GF2& rep(const mat_GF2& a)  
      { return a.rep; }  

  
   friend void swap(mat_GF2& X, mat_GF2& Y); 
  
   mat_GF2(mat_GF2& x, INIT_TRANS_TYPE)  :
      rep(x.rep, INIT_TRANS), numcols(x.numcols) { }
};  
  
void conv(mat_GF2& x, const vec_vec_GF2& a);  
inline mat_GF2 to_mat_GF2(const vec_vec_GF2& a)
   { mat_GF2 x; conv(x, a); NTL_OPT_RETURN(mat_GF2, x); }



long operator==(const mat_GF2& a, const mat_GF2& b); 
long operator!=(const mat_GF2& a, const mat_GF2& b); 


istream& operator>>(istream&, mat_GF2&); 
ostream& operator<<(ostream&, const mat_GF2&);  



void add(mat_GF2& X, const mat_GF2& A, const mat_GF2& B); 

inline void sub(mat_GF2& X, const mat_GF2& A, const mat_GF2& B)
   { add(X, A, B); }

inline void negate(mat_GF2& X, const mat_GF2& A)
   { X = A; }

void mul(mat_GF2& X, const mat_GF2& A, const mat_GF2& B); 
void mul(vec_GF2& x, const mat_GF2& A, const vec_GF2& b); 
void mul(vec_GF2& x, const vec_GF2& a, const mat_GF2& B); 

void mul(mat_GF2& X, const mat_GF2& A, GF2 b);
inline void mul(mat_GF2& X, GF2 a, const mat_GF2& B)
   { mul(X, B, a); }

inline void mul(mat_GF2& X, const mat_GF2& A, long b)
   { mul(X, A, to_GF2(b)); }
inline void mul(mat_GF2& X, long a, const mat_GF2& B)
   { mul(X, B, a); }

void ident(mat_GF2& X, long n); 
inline mat_GF2 ident_mat_GF2(long n)
   { mat_GF2 X; ident(X, n); NTL_OPT_RETURN(mat_GF2, X); }

long IsIdent(const mat_GF2& A, long n);
void transpose(mat_GF2& X, const mat_GF2& A);
void solve(GF2& d, vec_GF2& X, const mat_GF2& A, const vec_GF2& b);
void inv(GF2& d, mat_GF2& X, const mat_GF2& A);

inline void sqr(mat_GF2& X, const mat_GF2& A)
   { mul(X, A, A); }

inline mat_GF2 sqr(const mat_GF2& A)
   { mat_GF2 X; sqr(X, A); NTL_OPT_RETURN(mat_GF2, X); }

void inv(mat_GF2& X, const mat_GF2& A);

inline mat_GF2 inv(const mat_GF2& A)
   { mat_GF2 X; inv(X, A); NTL_OPT_RETURN(mat_GF2, X); }

void power(mat_GF2& X, const mat_GF2& A, const ZZ& e);
inline mat_GF2 power(const mat_GF2& A, const ZZ& e)
   { mat_GF2 X; power(X, A, e); NTL_OPT_RETURN(mat_GF2, X); }

inline void power(mat_GF2& X, const mat_GF2& A, long e)
   { power(X, A, ZZ_expo(e)); }
inline mat_GF2 power(const mat_GF2& A, long e)
   { mat_GF2 X; power(X, A, e); NTL_OPT_RETURN(mat_GF2, X); }


void diag(mat_GF2& X, long n, GF2 d);
inline mat_GF2 diag(long n, GF2 d)
   { mat_GF2 X; diag(X, n, d); NTL_OPT_RETURN(mat_GF2, X); }

long IsDiag(const mat_GF2& A, long n, GF2 d);


long gauss(mat_GF2& M);
long gauss(mat_GF2& M, long w);
void image(mat_GF2& X, const mat_GF2& A);
void kernel(mat_GF2& X, const mat_GF2& A);




void determinant(GF2& x, const mat_GF2& a);
inline GF2 determinant(const mat_GF2& a)
   { GF2 x; determinant(x, a); return x; }

inline mat_GF2 transpose(const mat_GF2 & a)
   { mat_GF2 x; transpose(x, a); NTL_OPT_RETURN(mat_GF2, x); }


void clear(mat_GF2& a);
// x = 0 (dimension unchanged)

long IsZero(const mat_GF2& a);
// test if a is the zero matrix (any dimension)


// operator notation:

mat_GF2 operator+(const mat_GF2& a, const mat_GF2& b);
mat_GF2 operator-(const mat_GF2& a, const mat_GF2& b);
mat_GF2 operator*(const mat_GF2& a, const mat_GF2& b);

inline mat_GF2 operator-(const mat_GF2& a)
   { return a; }


// matrix/scalar multiplication:

inline mat_GF2 operator*(const mat_GF2& a, GF2 b)
   { mat_GF2 x; mul(x, a, b); NTL_OPT_RETURN(mat_GF2, x); }

inline mat_GF2 operator*(const mat_GF2& a, long b)
   { mat_GF2 x; mul(x, a, b); NTL_OPT_RETURN(mat_GF2, x); }

inline mat_GF2 operator*(GF2 a, const mat_GF2& b)
   { mat_GF2 x; mul(x, a, b); NTL_OPT_RETURN(mat_GF2, x); }

inline mat_GF2 operator*(long a, const mat_GF2& b)
   { mat_GF2 x; mul(x, a, b); NTL_OPT_RETURN(mat_GF2, x); }




// matrix/vector multiplication:

vec_GF2 operator*(const mat_GF2& a, const vec_GF2& b);

vec_GF2 operator*(const vec_GF2& a, const mat_GF2& b);


// assignment operator notation:

inline mat_GF2& operator+=(mat_GF2& x, const mat_GF2& a)
{
   add(x, x, a);
   return x;
}   

inline mat_GF2& operator-=(mat_GF2& x, const mat_GF2& a)
{
   sub(x, x, a);
   return x;
}   


inline mat_GF2& operator*=(mat_GF2& x, const mat_GF2& a)
{
   mul(x, x, a);
   return x;
}   

inline mat_GF2& operator*=(mat_GF2& x, GF2 a)
{
   mul(x, x, a);
   return x;
}   

inline mat_GF2& operator*=(mat_GF2& x, long a)
{
   mul(x, x, a);
   return x;
}   
   

inline vec_GF2& operator*=(vec_GF2& x, const mat_GF2& a)
{
   mul(x, x, a);
   return x;
}   



#endif
