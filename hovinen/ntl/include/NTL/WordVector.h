
#ifndef NTL_WordVector__H
#define NTL_WordVector__H

/**************************************************************

  A WordVector is essentially functionally similar to
  a  generic NTL vector of u_long.  

  Be careful! the MaxLength() function does not return 
    the max length ever set, but rather the max space allocated,
    which *may* be more.

  The FixLength() facility is not available.

  The reason for special-casing is efficiency (of course).

**************************************************************/



#include <NTL/tools.h>
#include <NTL/ZZ.h>

typedef unsigned long u_long;
typedef u_long *u_long_ptr;



#ifndef NTL_RANGE_CHECK
#define NTL_WV_RANGE_CHECK_CODE 
#else
#define NTL_WV_RANGE_CHECK_CODE if (i < 0 || !rep || i >= long(rep[-1])) RangeError(i);
#endif

// vectors are allocated in chunks of this size

#ifndef NTL_WordVectorMinAlloc
#define NTL_WordVectorMinAlloc (4)
#endif

// vectors are always expanded by at least this ratio

#ifndef NTL_WordVectorExpansionRatio
#define NTL_WordVectorExpansionRatio (1.2)
#endif

// controls initialization during input

#ifndef NTL_WordVectorInputBlock
#define NTL_WordVectorInputBlock 50
#endif


class WordVector {  
public:  
   u_long *rep;  
   void RangeError(long i) const;  

   WordVector(WordVector& x, INIT_TRANS_TYPE) { rep = x.rep; x.rep = 0; }


  
   WordVector() { rep = 0; }  
   WordVector(INIT_SIZE_TYPE, long n) { rep = 0; DoSetLength(n); }  
   WordVector(const WordVector& a) { rep = 0; *this = a; }     

   WordVector& operator=(const WordVector& a);  

   ~WordVector();  
   void kill(); 

   void DoSetLength(long n);
  
   void SetLength(long n)
   {
      u_long *x = rep;
      if (x && long(x[-2] >> 1) >= n && n >= 0)
         x[-1] = n;
      else
         DoSetLength(n);
   }

   void ZeroLength() { if (rep) rep[-1] = 0; }
         
   void SetMaxLength(long n); 
   void QuickSetLength(long n) { rep[-1] = u_long(n); } 
  
   long length() const { return (!rep) ?  0 : long(rep[-1]); }  
   long MaxLength() const 
   { return (!rep) ?  0 : long(rep[-2] >> 1); } 
  
   u_long& operator[](long i)   
   {  
      NTL_WV_RANGE_CHECK_CODE  
      return rep[i];  
   }  
  
   const u_long& operator[](long i) const 
   {  
      NTL_WV_RANGE_CHECK_CODE  
      return rep[i];  
   }  
  
   u_long& operator()(long i) { return (*this)[i-1]; }  
   const u_long& operator()(long i) const { return (*this)[i-1]; } 
   
   friend void swap(WordVector& x, WordVector& y);  
  
   const u_long* elts() const { return rep; }  
   u_long* elts() { return rep; }  
         
   friend void append(WordVector& v, u_long a); 
   friend void append(WordVector& v, const WordVector& w); 
}; 




istream& operator>>(istream&, WordVector&);  
ostream& operator<<(ostream&, const WordVector&);  


long operator==(const WordVector& a, const WordVector& b);  
long operator!=(const WordVector& a, const WordVector& b);


long InnerProduct(const WordVector& a, const WordVector& b);

void ShiftAdd(u_long *cp, const u_long* ap, long sa, long n);
// cp = cp + (a << n)



#endif
