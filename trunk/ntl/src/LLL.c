
#include <NTL/LLL.h>
#include <NTL/vec_long.h>
#include <NTL/tools.h>


static void ExactDiv(ZZ& qq, const ZZ& a, const ZZ& b)
{
   static ZZ q, r;

   DivRem(q, r, a, b);
   if (!IsZero(r)) {
      cerr << "a = " << a << "\n";
      cerr << "b = " << b << "\n";
      Error("ExactDiv: nonzero remainder");
   }
   qq = q;
}


static void BalDiv(ZZ& q, const ZZ& a, const ZZ& d)

//  rounds a/d to nearest integer, breaking ties
//    by rounding towards zero.  Assumes d > 0.

{
   static ZZ r;
   DivRem(q, r, a, d);


   add(r, r, r);

   long cmp = compare(r, d);
   if (cmp > 0 || (cmp == 0 && q < 0))
      add(q, q, 1);
}



static void MulAddDiv(ZZ& c, const ZZ& c1, const ZZ& c2, 
                      const ZZ& x, const ZZ& y, const ZZ& z)

// c = (x*c1 + y*c2)/z

{
   static ZZ t1, t2;

   mul(t1, x, c1);
   mul(t2, y, c2);
   add(t1, t1, t2);
   ExactDiv(c, t1, z);
}


static void MulSubDiv(ZZ& c, const ZZ& c1, const ZZ& c2, 
                      const ZZ& x, const ZZ& y, const ZZ& z)

// c = (x*c1 - y*c2)/z

{
   static ZZ t1, t2;

   mul(t1, x, c1);
   mul(t2, y, c2);
   sub(t1, t1, t2);
   ExactDiv(c, t1, z);
}
   




static void MulSubDiv(vec_ZZ& c, const vec_ZZ& c1, const vec_ZZ& c2,
                      const ZZ& x, const ZZ& y, const ZZ& z)

// c = (x*c1 + y*c2)/z

{
   long n = c1.length();
   if (c2.length() != n) Error("MulSubDiv: length mismatch");
   c.SetLength(n);

   long i;
   for (i = 1; i <= n; i++) 
      MulSubDiv(c(i), c1(i), c2(i), x, y, z);
}

static void RowTransform(vec_ZZ& c1, vec_ZZ& c2,
                         const ZZ& x, const ZZ& y, const ZZ& u, const ZZ& v)

// (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2)

{
   long n = c1.length();
   if (c2.length() != n) Error("MulSubDiv: length mismatch");
   static ZZ t1, t2, t3, t4;

   long i;
   for (i = 1; i <= n; i++) {
      mul(t1, x, c1(i));
      mul(t2, y, c2(i));
      add(t1, t1, t2);

      mul(t3, u, c1(i));
      mul(t4, v, c2(i));
      add(t3, t3, t4);

      c1(i) = t1;
      c2(i) = t3;
   }
}

static void RowTransform(ZZ& c1, ZZ& c2,
                         const ZZ& x, const ZZ& y, const ZZ& u, const ZZ& v)

// (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2)

{
   static ZZ t1, t2, t3, t4;

   mul(t1, x, c1);
   mul(t2, y, c2);
   add(t1, t1, t2);

   mul(t3, u, c1);
   mul(t4, v, c2);
   add(t3, t3, t4);

   c1 = t1;
   c2 = t3;
}



static void MulSub(ZZ& c, const ZZ& c1, const ZZ& c2, const ZZ& x)

// c = c1 - x*c2

{
   static ZZ t1;

   mul(t1, x, c2);
   sub(c, c1, t1);
}


static void MulSub(vec_ZZ& c, const vec_ZZ& c1, const vec_ZZ& c2,
                   const ZZ& x)

// c = c1 - x*c2

{
   long n = c1.length();
   if (c2.length() != n) Error("MulSub: length mismatch");
   c.SetLength(n);

   long i;
   for (i = 1; i <= n; i++)
      MulSub(c(i), c1(i), c2(i), x);
}


      
      
   
static long SwapTest(const ZZ& d0, const ZZ& d1, const ZZ& d2, const ZZ& lam,
                     long a, long b)

// test if a*d1^2 > b*(d0*d2 + lam^2)

{
   static ZZ t1, t2;

   mul(t1, d0, d2);
   sqr(t2, lam);
   add(t1, t1, t2);
   mul(t1, t1, b);

   sqr(t2, d1);
   mul(t2, t2, a);

   return t2 > t1;
}






static
void reduce(long k, long l, 
            mat_ZZ& B, vec_long& P, vec_ZZ& D, 
            vec_vec_ZZ& lam, mat_ZZ* U)
{
   static ZZ t1;
   static ZZ r;

   if (P(l) == 0) return;
   add(t1, lam(k)(P(l)), lam(k)(P(l)));
   abs(t1, t1);
   if (t1 <= D[P(l)]) return;

   long j;

   BalDiv(r, lam(k)(P(l)), D[P(l)]);
   MulSub(B(k), B(k), B(l), r);

   if (U) MulSub((*U)(k), (*U)(k), (*U)(l), r);

   for (j = 1; j <= l-1; j++)
      if (P(j) != 0)
         MulSub(lam(k)(P(j)), lam(k)(P(j)), lam(l)(P(j)), r);

   MulSub(lam(k)(P(l)), lam(k)(P(l)), D[P(l)], r);
}

static
void swap(long k, mat_ZZ& B, vec_long& P, vec_ZZ& D, 
          vec_vec_ZZ& lam, mat_ZZ* U, long m, long verbose)

// swaps vectors k-1 and k;  assumes P(k-1) != 0

{
   long i, j;
   static ZZ t1, t2, t3, e, x, y;


   if (P(k) != 0) {
      if (verbose) cerr << "swap case 1: " << k << "\n";

      swap(B(k-1), B(k));
      if (U) swap((*U)(k-1), (*U)(k));
   
      for (j = 1; j <= k-2; j++)
         if (P(j) != 0)
            swap(lam(k-1)(P(j)), lam(k)(P(j)));

      for (i = k+1; i <= m; i++) {
         MulAddDiv(t1, lam(i)(P(k)-1), lam(i)(P(k)), 
                   lam(k)(P(k)-1), D[P(k)-2], D[P(k)-1]); 
         MulSubDiv(t2, lam(i)(P(k)-1), lam(i)(P(k)), 
                   D[P(k)], lam(k)(P(k)-1), D[P(k)-1]);
         lam(i)(P(k)-1) = t1;
         lam(i)(P(k)) = t2;
      }

      MulAddDiv(D[P(k)-1], D[P(k)], lam(k)(P(k)-1),
                D[P(k)-2], lam(k)(P(k)-1), D[P(k)-1]);
   }
   else if (!IsZero(lam(k)(P(k-1)))) {
      if (verbose) cerr << "swap case 2: " << k << "\n";
      XGCD(e, x, y, lam(k)(P(k-1)), D[P(k-1)]);

      ExactDiv(t1, lam(k)(P(k-1)), e);
      ExactDiv(t2, D[P(k-1)], e);

      t3 = t2;
      negate(t2, t2);
      RowTransform(B(k-1), B(k), t1, t2, y, x);
      if (U) RowTransform((*U)(k-1), (*U)(k), t1, t2, y, x);
      for (j = 1; j <= k-2; j++)
         if (P(j) != 0)
            RowTransform(lam(k-1)(P(j)), lam(k)(P(j)), t1, t2, y, x);

      sqr(t2, t2);
      ExactDiv(D[P(k-1)], D[P(k-1)], t2);

      for (i = k+1; i <= m; i++)
         if (P(i) != 0) {
            ExactDiv(D[P(i)], D[P(i)], t2);
            for (j = i+1; j <= m; j++) {
               ExactDiv(lam(j)(P(i)), lam(j)(P(i)), t2);
            }
         }

      for (i = k+1; i <= m; i++) {
         ExactDiv(lam(i)(P(k-1)), lam(i)(P(k-1)), t3);
      }

      swap(P(k-1), P(k));
   }
   else {
      if (verbose) cerr << "swap case 3: " << k << "\n";

      swap(B(k-1), B(k));
      if (U) swap((*U)(k-1), (*U)(k));
   
      for (j = 1; j <= k-2; j++)
         if (P(j) != 0)
            swap(lam(k-1)(P(j)), lam(k)(P(j)));

      swap(P(k-1), P(k));
   }
}

   


static
void IncrementalGS(mat_ZZ& B, vec_long& P, vec_ZZ& D, vec_vec_ZZ& lam, 
                   long& s, long k)
{
   long n = B.NumCols();
   long m = B.NumRows();

   static ZZ u, t1, t2;

   long i, j;

   for (j = 1; j <= k-1; j++) {
      long posj = P(j);
      if (posj == 0) continue;

      InnerProduct(u, B(k), B(j));
      for (i = 1; i <= posj-1; i++) {
         mul(t1, D[i], u);
         mul(t2, lam(k)(i), lam(j)(i));
         sub(t1, t1, t2);
         div(t1, t1, D[i-1]);
         u = t1;
      }

      lam(k)(posj) = u;
   }

   InnerProduct(u, B(k), B(k));
   for (i = 1; i <= s; i++) {
      mul(t1, D[i], u);
      mul(t2, lam(k)(i), lam(k)(i));
      sub(t1, t1, t2);
      div(t1, t1, D[i-1]);
      u = t1;
   }

   if (u == 0) {
      P(k) = 0;
   }
   else {
      s++;
      P(k) = s;
      D[s] = u;
   }
}


static
long LLL(ZZ& det, mat_ZZ& B, mat_ZZ* U, long a, long b, long verbose)
{
   long m = B.NumRows();
   long n = B.NumCols();

   vec_long P;
   P.SetLength(m);

   vec_ZZ D;
   D.SetLength(m+1);
   D[0] = 1;

   vec_vec_ZZ lam;

   lam.SetLength(m);

   long j;
   for (j = 1; j <= m; j++)
      lam(j).SetLength(m);

   if (U) ident(*U, m);

   long s = 0;

   long k = 1;
   long max_k = 0;


   while (k <= m) {
      if (k > max_k) {
         IncrementalGS(B, P, D, lam, s, k);
         max_k = k;
      }

      if (k == 1) {
         k++;
         continue;
      }

      reduce(k, k-1, B, P, D, lam, U);

      if (P(k-1) != 0 && 
          (P(k) == 0 || 
           SwapTest(D[P(k)], D[P(k)-1], D[P(k)-2], lam(k)(P(k)-1), a, b))) {
         swap(k, B, P, D, lam, U, max_k, verbose);
         k--;
      }
      else {
         for (j = k-2; j >= 1; j--) 
            reduce(k, j, B, P, D, lam, U);
         k++;
      }
   }

   det = D[s];
   return s;
}

static
long image(ZZ& det, mat_ZZ& B, mat_ZZ* U, long verbose)
{
   long m = B.NumRows();
   long n = B.NumCols();

   vec_long P;
   P.SetLength(m);

   vec_ZZ D;
   D.SetLength(m+1);
   D[0] = 1;

   vec_vec_ZZ lam;

   lam.SetLength(m);

   long j;
   for (j = 1; j <= m; j++)
      lam(j).SetLength(m);

   if (U) ident(*U, m);

   long s = 0;

   long k = 1;
   long max_k = 0;


   while (k <= m) {
      if (k > max_k) {
         IncrementalGS(B, P, D, lam, s, k);
         max_k = k;
      }

      if (k == 1) {
         k++;
         continue;
      }

      reduce(k, k-1, B, P, D, lam, U);

      if (P(k-1) != 0 && P(k) == 0) {
         swap(k, B, P, D, lam, U, max_k, verbose);
         k--;
      }
      else {
         for (j = k-2; j >= 1; j--) 
            reduce(k, j, B, P, D, lam, U);
         k++;
      }
   }

   det = D[s];
   return s;
}

long LLL(ZZ& det, mat_ZZ& B, mat_ZZ& U, long verbose)
{
   return LLL(det, B, &U, 3, 4, verbose);
}

long LLL(ZZ& det, mat_ZZ& B, long verbose)
{
   return LLL(det, B, 0, 3, 4, verbose);
}

long LLL(ZZ& det, mat_ZZ& B, mat_ZZ& U, long a, long b, long verbose)
{
   if (a <= 0 || b <= 0 || a > b || b/4 >= a) Error("LLL: bad args");
   
   return LLL(det, B, &U, a, b, verbose);
}

long LLL(ZZ& det, mat_ZZ& B, long a, long b, long verbose)
{
   if (a <= 0 || b <= 0 || a > b || b/4 >= a) Error("LLL: bad args");

   return LLL(det, B, 0, a, b, verbose);
}


long image(ZZ& det, mat_ZZ& B, mat_ZZ& U, long verbose)
{
   return image(det, B, &U, verbose);
}

long image(ZZ& det, mat_ZZ& B, long verbose)
{
   return image(det, B, 0, verbose);
}

