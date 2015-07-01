#ifndef NTL_matrix__H
#define NTL_matrix__H

#include <NTL/tools.h>
#include <NTL/vector.h>


// matrix templates


#define NTL_matrix_decl(T,vec_T,vec_vec_T,mat_T)  \
class mat_T {  \
private:  \
  \
   vec_vec_T l__rep;  \
   long l__numcols;  \
  \
public:  \
  \
   mat_T() { l__numcols = 0; }  \
   mat_T(const mat_T& l__a);  \
   mat_T& operator=(const mat_T& l__a);  \
   ~mat_T() { }  \
  \
   mat_T(INIT_SIZE_TYPE, long l__n, long l__m);  \
  \
   void kill();  \
  \
   void SetDims(long l__n, long l__m);  \
  \
   long NumRows() const { return l__rep.length(); }  \
   long NumCols() const { return l__numcols; }  \
  \
   vec_T& operator[](long l__i) { return l__rep[l__i]; }  \
   const vec_T& operator[](long l__i) const { return l__rep[l__i]; }  \
  \
   vec_T& operator()(long l__i) { return l__rep[l__i-1]; }  \
   const vec_T& operator()(long l__i) const { return l__rep[l__i-1]; }  \
  \
   T& operator()(long l__i, long l__j) { return l__rep[l__i-1][l__j-1]; }  \
   const T& operator()(long l__i, long l__j) const   \
      { return l__rep[l__i-1][l__j-1]; }  \
  \
   friend const vec_vec_T& rep(const mat_T& l__a)  \
      { return l__a.l__rep; }  \
  \
   friend void swap(mat_T& l__X, mat_T& l__Y); \
  \
   long position(const vec_T& l__a) const { return l__rep.position(l__a); } \
  mat_T(mat_T& l__x, INIT_TRANS_TYPE) :  \
    l__rep(l__x.l__rep, INIT_TRANS), l__numcols(l__x.l__numcols) { }  \
};  \
  \
void MakeMatrix(mat_T& l__x, const vec_vec_T& l__a);  \



#define NTL_eq_matrix_decl(T,vec_T,vec_vec_T,mat_T) \
long operator==(const mat_T& l__a, const mat_T& l__b); \
long operator!=(const mat_T& l__a, const mat_T& l__b); \



#define NTL_io_matrix_decl(T,vec_T,vec_vec_T,mat_T) \
istream& operator>>(istream&, mat_T&); \
ostream& operator<<(ostream&, const mat_T&);  \


#define NTL_matrix_impl(T,vec_T,vec_vec_T,mat_T)  \
mat_T::mat_T(const mat_T& l__a)  \
{  \
   l__numcols = 0;  \
   SetDims(l__a.NumRows(), l__a.NumCols());  \
   l__rep = l__a.l__rep;  \
}  \
  \
mat_T& mat_T::operator=(const mat_T& l__a)  \
{  \
   SetDims(l__a.NumRows(), l__a.NumCols());  \
   l__rep = l__a.l__rep;  \
   return *this;  \
}  \
  \
  \
mat_T::mat_T(INIT_SIZE_TYPE, long l__n, long l__m)  \
{  \
   l__numcols = 0;  \
   SetDims(l__n, l__m);  \
}  \
  \
void mat_T::kill()  \
{  \
   l__numcols = 0;  \
   l__rep.kill();  \
}  \
  \
void mat_T::SetDims(long l__n, long l__m)  \
{  \
   if (l__n < 0 || l__m < 0)  \
      Error("SetDims: bad args");  \
  \
   if (l__m != l__numcols) {  \
      l__rep.kill();  \
      l__numcols = l__m;  \
   }  \
        \
   long l__oldmax = l__rep.MaxLength();  \
   long l__i;  \
   l__rep.SetLength(l__n);  \
  \
   for (l__i = l__oldmax; l__i < l__n; l__i++)  \
      l__rep[l__i].FixLength(l__m);  \
}  \
     \
        \
void MakeMatrix(mat_T& l__x, const vec_vec_T& l__a)  \
{  \
   long l__n = l__a.length();  \
  \
   if (l__n == 0) {  \
      l__x.SetDims(0, 0);  \
      return;  \
   }  \
  \
   long l__m = l__a[0].length();  \
   long l__i;  \
  \
   for (l__i = 1; l__i < l__n; l__i++)  \
      if (l__a[l__i].length() != l__m)  \
         Error("nonrectangular matrix");  \
  \
   l__x.SetDims(l__n, l__m);  \
   for (l__i = 0; l__i < l__n; l__i++)  \
      l__x[l__i] = l__a[l__i];  \
}  \
  \
void swap(mat_T& l__X, mat_T& l__Y)  \
{  \
   swap(l__X.l__numcols, l__Y.l__numcols);  \
   swap(l__X.l__rep, l__Y.l__rep);  \
}  \
  \




   

#define NTL_eq_matrix_impl(T,vec_T,vec_vec_T,mat_T)  \
long operator==(const mat_T& l__a, const mat_T& l__b)  \
{  \
   if (l__a.NumCols() != l__b.NumCols())  \
      return 0;  \
  \
   if (l__a.NumRows() != l__b.NumRows())  \
      return 0;  \
  \
   long l__n = l__a.NumRows();  \
   long l__i;  \
  \
   for (l__i = 0; l__i < l__n; l__i++)  \
      if (l__a[l__i] != l__b[l__i])  \
         return 0;  \
  \
   return 1;  \
}  \
  \
  \
long operator!=(const mat_T& l__a, const mat_T& l__b)  \
{  \
   return !(l__a == l__b);  \
}  \




#define NTL_io_matrix_impl(T,vec_T,vec_vec_T,mat_T)  \
istream& operator>>(istream& l__s, mat_T& l__x)  \
{  \
   vec_vec_T l__buf;  \
   l__s >> l__buf;  \
   MakeMatrix(l__x, l__buf);  \
   return l__s;  \
}  \
  \
ostream& operator<<(ostream& l__s, const mat_T& l__a)  \
{  \
   long l__n = l__a.NumRows();  \
   long l__i;  \
   l__s << "[";  \
   for (l__i = 0; l__i < l__n; l__i++) {  \
      l__s << l__a[l__i]; \
      l__s << "\n"; \
   }  \
   l__s << "]";  \
   return l__s;  \
}  \




#endif
