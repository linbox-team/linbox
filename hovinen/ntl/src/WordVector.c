
#include <NTL/WordVector.h>



void WordVector::DoSetLength(long n)   
{   
   long m;  
  
   if (n < 0) {  
      Error("negative length in vector::SetLength");  
   }  

   if (n >= (1L << (NTL_BITS_PER_LONG-4))/NTL_BITS_PER_LONG) {
      Error("length too big in vector::SetLength");
   }
      
   if (n == 0) {  
      if (rep) rep[-1] = 0;  
      return;  
   }  
  
   if (!rep) {  
      m = ((n+NTL_WordVectorMinAlloc-1)/NTL_WordVectorMinAlloc) * NTL_WordVectorMinAlloc; 
      u_long *p = (u_long *) malloc(sizeof(u_long)*(m+2)); 
      if (!p) {  
	 Error("out of memory in SetLength()");  
      }  
      rep = p+2;

      rep[-1] = n;
      rep[-2] = m << 1;
 
      return;
   }  

   long max_length = (rep[-2] >> 1);


   if (n <= max_length) {  
      rep[-1] = n;  
      return;
   }  

   long frozen = (rep[-2] & 1);

   if (frozen) Error("Cannot grow this WordVector");
      
   if (n < NTL_WordVectorExpansionRatio*max_length)  
      m = (long) (NTL_WordVectorExpansionRatio*max_length);  
   else  
      m = n;  

   m = ((m+NTL_WordVectorMinAlloc-1)/NTL_WordVectorMinAlloc)*NTL_WordVectorMinAlloc; 
   u_long *p = rep - 2;
   p = (u_long *) realloc(p, sizeof(u_long)*(m+2)); 
   if (!p) {  
      Error("out of memory in SetLength()");  
   }  
   rep = p+2;

   rep[-1] = n;
   rep[-2] = m << 1;
}  
 
 
void WordVector::SetMaxLength(long n) 
{ 
   long OldLength = length(); 
   DoSetLength(n); 
   if (rep) rep[-1] = OldLength;
} 
 
  
WordVector& WordVector::operator=(const WordVector& a)  
{  
   long i, n;  
   u_long *p;  
   const u_long *ap;  
  
   if (this == &a) return *this;  
  
   n = a.length();  
   ap = a.elts();  
  
   SetLength(n);  
   p = elts();  
  
  
   for (i = 0; i < n; i++)  
      p[i] = ap[i];  

   return *this;
}  
       
  
WordVector::~WordVector()  
{  
   if (!rep) return;  
   if (rep[-2] & 1) Error("Cannot free this WordVector");
   free(rep-2);
}  
   
void WordVector::kill()  
{  
   if (!rep) return;  
   if (rep[-2] & 1) 
      Error("Cannot free this WordVector");
   free(rep-2);
   rep = 0; 
}  
  
void WordVector::RangeError(long i) const  
{  
   cerr << "index out of range in vector: ";  
   cerr << i;  
   if (!rep)  
      cerr << "(0)\n";  
   else  
      cerr << "(" << rep[-1] << ")\n";  
   abort();  
}  

void CopySwap(WordVector& x, WordVector& y)
{
   static WordVector t;
   t = x;
   x = y;
   y = t;
}
 
void swap(WordVector& x, WordVector& y)  
{  
   if ((x.rep && (x.rep[-2] & 1)) ||
       (y.rep && (y.rep[-2] & 1))) {
      CopySwap(x, y);
      return;
   }

   u_long* t;  
   t = x.rep;  
   x.rep = y.rep;  
   y.rep = t;  
} 
 
void append(WordVector& v, u_long a)  
{  
   long l = v.length();
   v.SetLength(l+1);  
   v[l] = a;  
}  
  
void append(WordVector& v, const WordVector& w)  
{  
   long l = v.length();  
   long m = w.length();  
   long i;  
   v.SetLength(l+m);  
   for (i = 0; i < m; i++)  
      v[l+i] = w[i];  
}


istream & operator>>(istream& s, WordVector& a)   
{   
   WordVector ibuf;  
   long c;   
   long n;   
   if (!s) Error("bad vector input"); 
   
   c = s.peek();  
   while (c == ' ' || c == '\n' || c == '\t') {  
      s.get();  
      c = s.peek();  
   }  
   if (c != '[') {  
      Error("bad vector input");  
   }  
 
   n = 0;   
   ibuf.SetLength(0);  
      
   s.get();  
   c = s.peek();  
   while (c == ' ' || c == '\n' || c == '\t') {  
      s.get();  
      c = s.peek();  
   }  
   while (c != ']' && c != EOF) {   
      if (n % NTL_WordVectorInputBlock == 0) ibuf.SetMaxLength(n + NTL_WordVectorInputBlock); 
      n++;   
      ibuf.SetLength(n);   
      if (!(s >> ibuf[n-1])) Error("bad vector input");   
      c = s.peek();  
      while (c == ' ' || c == '\n' || c == '\t') {  
         s.get();  
         c = s.peek();  
      }  
   }   
   if (c == EOF) Error("bad vector input");  
   s.get(); 
   
   a = ibuf; 
   return s;   
}    
   
   
ostream& operator<<(ostream& s, const WordVector& a)   
{   
   long i, n;   
  
   n = a.length();  
   
   s << '[';   
   
   for (i = 0; i < n; i++) {   
      s << a[i];   
      if (i < n-1) s << " ";   
   }   
   
   s << ']';   
      
   return s;   
}   

long operator==(const WordVector& a, const WordVector& b) 
{  
   long n = a.length();  
   if (b.length() != n) return 0;  
   const u_long* ap = a.elts(); 
   const u_long* bp = b.elts(); 
   long i;  
   for (i = 0; i < n; i++) if (ap[i] != bp[i]) return 0;  
   return 1;  
} 

long operator!=(const WordVector& a, const WordVector& b) 
{  return !(a == b); }

   



long InnerProduct(const WordVector& a, const WordVector& b)
{
   long n = min(a.length(), b.length());
   const u_long *ap = a.elts();
   const u_long *bp = b.elts();

   u_long acc;
   long i;

   acc = 0;
   for (i = 0; i < n; i++)
      acc ^= ap[i] & bp[i];

#if (NTL_BITS_PER_LONG == 32)
   acc ^= acc >> 16;
   acc ^= acc >> 8;
   acc ^= acc >> 4;
   acc ^= acc >> 2;
   acc ^= acc >> 1;
   acc &= 1;
#elif (NTL_BITS_PER_LONG == 64)
   acc ^= acc >> 32;
   acc ^= acc >> 16;
   acc ^= acc >> 8;
   acc ^= acc >> 4;
   acc ^= acc >> 2;
   acc ^= acc >> 1;
   acc &= 1;
#else
   u_long t = acc;
   while (t) {
      t = t >> 8;
      acc ^= t;
   }

   acc ^= acc >> 4;
   acc ^= acc >> 2;
   acc ^= acc >> 1;
   acc &= 1;
#endif

   return long(acc);
}


void ShiftAdd(u_long *cp, const u_long* ap, long sa, long n)
// c = c + (a << n)
{
   if (sa == 0) return;

   long i;

   long wn = n/NTL_BITS_PER_LONG;
   long bn = n - wn*NTL_BITS_PER_LONG;

   if (bn == 0) {
      for (i = sa+wn-1; i >= wn; i--)
         cp[i] ^= ap[i-wn];
   }
   else {
      u_long t = ap[sa-1] >> (NTL_BITS_PER_LONG-bn);
      if (t) cp[sa+wn] ^= t;
      for (i = sa+wn-1; i >= wn+1; i--)
         cp[i] ^= (ap[i-wn] << bn) | (ap[i-wn-1] >> (NTL_BITS_PER_LONG-bn));
      cp[wn] ^= ap[0] << bn;
   }
}
