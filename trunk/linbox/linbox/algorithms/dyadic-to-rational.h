#ifndef __DYADICTORATIONAL_H
#define __DYADICTORATIONAL_H
/* dyadic-to-rational.h
 *
 * dyadicToRational reconstructs a rational a/b from dyadic approximation n/2^k.
 *
 * It is used in rational-solver-sn
 *
 * "Rational reconstruction" starts from a p-adic approximation.
 * This is different though very similar
 *
 *  Evolved by bds from rational-reconstruction2.h by Z. Wan.
 *
 *  standard linbox copyright/license applies.  See COPYING.
 */

#include <stack>
#include <assert.h>
#include <vector>

//#include <linbox/integer.h>

namespace LinBox{

/** Rational reconstruction of a/b from n/d with denominator bound B.
 * We give a/b, the continued fraction approximant of n/d that
 * satisfies |a/b - n/d| < 1/2d (well approximated) and 0 < b <= B.
 * Return value is 0, if no such approximant exists.
 * Return value is 1, if either
 *   (i) a second well approximated rational with denominator bounded by B may exist, or
 *   (ii) the well approximated condition is not met for a/b.
 *   In these cases, a/b may be used speculatively.
 * Return value is 2, if the approximant is guaranteed (because bB <= d).

 * If no fraction is well approximated the last b <= B in the remainder sequence of n,d is given.
 *
 * If d = 2^k and n = sum_i=l to k n_i 2^i, then * n/d = sum_{i=l down to 0} n_i/2^{k-i}
 * is a {\em dyadic rational}.  Numbers of this form are produced for example by
 * numeric-symbolic iterations.
 *
 * If it is known that n/d is the most accurate approximation with denominator d
 * to a/b, and that the denominator b is bounded by B, i.e. b <= B, then such a/b is
 * uniquely determined, provided d >= bB.
 * ...in that case, such a/b is returned by dyadicToRational().
 * This follows from two facts:
 * First, by definition, n/d is an accurate approximation to a/b
 * with b <= d when |n/d - a/b| < 1/2d.
 * Otherwise (n-1)/d or (n+1)/d would be a better approximation.
 * Second, if a/b and a'/b' are distinct rationals, then |a/b - a'/b'| >= 1/bb'.
 * Thus if a'/b' is another rational accurately approximated by n/d,
 * we have 1/bb' <= |a/b - a'/b'| <= |a/b - n/d| + |n/d - a'/b'| <= 1/2d + 1/2d = 1/d.
 * So bb' > d >= bB, thus b' > B.
 *
 * In summary: If it exists, the unique a/b is given such that n/d approximates a/b
 * to within 1/2d and b <= B.  Otherwise a plausible a/b is given or failure is signaled.
 *
 * "Symbolic-Numeric Exact Rational Linear System Solver" by Saunders, Wood, Youse.
 * describes the construction.
 */
template<class Ring>
int dyadicToRational (
	const Ring& Z,
	typename Ring::Element& a, typename Ring::Element& b,
	const typename Ring::Element& n, const typename Ring::Element& d,
	const typename Ring::Element& B)
{
	typedef typename Ring::Element Int;
	Int e; Z.init(e);// error term

	Int an; Z.init(an); Z.abs(an, n/*q*/);

	// Partial_hegcd is defined below.
	bool found = partial_hegcd(Z, e, b, an, d/*q*/, B); // e = b*an - ?*d and |b| <= B
	Z.axmyin(e, b, an);
	Z.div(a, e, d); //a = (e - b*an)/d, div is exact.
//std::cout << "Z.axmyin(e, b, an); " << e << " " << b << " " << an << ", a = e/d exact " << a << " " << d << std::endl;
	// now a/b is solution but signs may be wrong
	Z.abs(a,a);
	Z.abs(b,b);
	Int zero; Z.init(zero, 0);
	if (Z.compare(n, zero) < 0)  Z.negin(a); // a = -a;

//std::cout << "DtR in n, d " << n << " "<< d << ", bound " << B << ", out a, b " << a << " " << b << std::endl;
	bool guarantee = b*B < d;
	if (found && guarantee) return 2;
	if (b == 0) return 0;
	return 1; //if ((!found && b > 0) || (found && ! guarantee)) return 1;
}

/** partial_hegcd() sets e, b from the remainder sequence of n,d.
 * It requires positive n and d.
 * It sets e to the first r_i (remainder) and
 * b to the corresponding q_i (coefficient of n)
 * such that 2r_i < |q_i| and |q_i| <= B (the given denominator bound).
 * True is returned iff such e, b exist.
 *
 * If not, b is the largest q_i such that |q_i| <= B,
 * and e is the corresponding remainder.  In this case b is the denominator
 * of a plausibly approximated but not well approximated rational.
 * It can be used speculatively.
 */
// later reintroduce swapping trick
template<class Ring>
bool partial_hegcd(Ring& Z, typename Ring::Element& e, typename Ring::Element& b, const typename Ring::Element& n, const typename Ring::Element& d, const typename Ring::Element& denBound){
	typedef typename Ring::Element Int;
	Int quo, r, tmp;  Z.init(quo); Z.init(r); Z.init(tmp);
	bool withinbound, wellapproximated;

	Int b0; Z.init(b0, 1); // and a0 = -0
	Int r0; Z.init(r0, n); // so that r0 = b0*n - a0*d
	Int b1; Z.init(b1, -0); // and a1 = 1
	Int r1; Z.init(r1, d); // so that r1 = b1*n - a1*d
//std::cout << "init 1 -0: " << b0 << " " << b1 << std::endl;

	do {
//std::cout << "quorem from " << r0 << " " << b0 << " and " << r1 << " " << b1 << std::endl;
		Z.quoRem(quo, e, r0, r1);
		//integer::divmod (quo, e, r0, r1);
		b = b0;
		Z.axmyin(b, quo, b1);
		Z.negin(b); // b = b0 - quo*b1;
//std::cout << "Z.axmyin(b, quo, b1);// b = b0 - quo*b1;" << b << " " << b0 << " " << quo << " " << b1 << std::endl;
		r0 = r1; b0 = b1;
		r1 = e; b1 = b;
		//assert(r1 >= 0);
		Z.abs(tmp, b);
		withinbound = (Z.compare(tmp, denBound) <= 0);
	    wellapproximated = (Z.compare(2*e , tmp) <= 0);
	}
	while ( ! wellapproximated && withinbound );
	if (! withinbound) {e = r0; b = b0;} // make available for speculation
//std::cout << "withinbound " << withinbound << " e b " << e << " " << b << ", n/d " << n << "/" << d << ", denBound " << denBound << std::endl;
	return withinbound;
	// returning with first well approximated (small remainder e) or last within bound denom b.

} // partial_hegcd

// vector rational reconstruction building num, den from numx, denx
template<class Ring>
int dyadicToRational(
	const Ring& Z,
	std::vector<typename Ring::Element>& num, typename Ring::Element& den,
	std::vector<typename Ring::Element>& numx, typename Ring::Element& denx,
	typename Ring::Element& denBound)
{
	typedef typename Ring::Element Int;
	Int q, rem, tmp_den, nx;
	Z.init(q); Z.init(rem); Z.init(tmp_den); Z.init(nx);
	Int one; Z.init(one, 1);
	Int two; Z.init(two, 2);
	Int denx2;
	Z.init(denx2, denx); Z.divin(denx2, two);// denx2 = denx/2, for balancing remainders.
	std::stack<std::pair<size_t, Int> > S;
	Int tmp; Z.init(tmp);

	Int den_lcm; Z.init(den_lcm, 1);
	den = den_lcm; // = 1.

	S.push(std::pair<int, Int>(0, 1));
	Int e; Z.init(e);// e for error
	int ret = 2; // 2 means guaranteed, 1 possible, 0 fail.
	for (size_t i = 0; i < num.size(); ++i) {
		Z.abs(nx, numx[i]);
		Z.mul(tmp, nx, den);
		Z.quoRem(num[i], e, tmp, denx); //nx*den - num[i]*denx = e, with num[i] and e nonneg.
		// we need |nx/denx - num[i]/den| == e/den*denx <= 1/2denx, so 2e <= den.
		// adjust to balanced remainder e.
		if (Z.compare(e, denx2) >= 0) {Z.subin(e, denx), Z.addin(num[i], one); }
		//nx*den = num[i]*denx + e , thus |nx/denx - num[i]/den| = e/denx*den

	// can try e < den && 2*e < den for speed
		Z.mul(tmp, two, Z.abs(e));
		if ( Z.compare(tmp, den) > 0)// 2|e| > den
		{   // we failed, so need another reconstruction
			int oneret = dyadicToRational (Z, tmp, tmp_den, nx, denx, denBound);
			if (oneret == 1) ret = 1;
			if ( oneret == 0 ) return oneret;
			//std::cerr << i << " tmp " << tmp << " num[i] " << num[i] << std::endl;
			num[i] = tmp;

			Z.lcm (den_lcm, tmp_den, den);
			Z.div( tmp, den_lcm, tmp_den ); // exact
			Z.mulin( num[i], tmp ); // num[i]/den_lcm = prev num[i]/tmp_den.

			Z.div(tmp, den_lcm, den); // exact
			// must multiply all prior num[i] by this tmp.
			S.push(std::pair<size_t, Int>(i, tmp));
			den = den_lcm;
			//assert(Z.compare(den, denBound)<=0);
			if (Z.compare(den, denBound)>0) return false; // den > denBound
		}

		Int zero; Z.init(zero);
		if (Z.compare(numx[i], zero) < 0) Z.negin(num[i]); // numx[i] < 0

	}
	// now fix shorties
	Int t; Z.init(t, 1);
	while ( S.size() > 1 ) {
		Z.mulin(t, S.top().second);
		int k = (int)S.top().first;
		S.pop();
		int j = (int)S.top().first;
		for (int i = k-1; i >= j; --i) {
			Z.mulin(num[i], t);
		}
	}
	S.pop();
	return ret;

} // vector dyadicToRational

#if 0
// vector rational reconstruction building num, den from numx, denx
// This one -- very inefficient -- just reconstructs each one, then goes thru to fix for lcm.
void rational_reconstruction(std::vector<integer>& num, integer& den, std::vector<integer>& numx, integer& denx, integer& denBound) {
	integer den_tmp, missing_factor;
	den = 1;
	for (size_t i = 0; i < numx.size(); ++i) {
		rational_reconstruction(num[i], den_tmp, numx[i], denx, denBound);
		lcm(missing_factor, den_tmp, den);
		den = missing_factor;
	}
	for (size_t i = 0; i < numx.size(); ++i) {
		rational_reconstruction(num[i], den_tmp, numx[i], denx, denBound);
		integer::divexact (missing_factor, den, den_tmp);
		num[i] *= missing_factor;
	}
} // vector rational_reconstruction
****
#endif

}// LinBox
#endif // __DYADICTORATIONAL_H


