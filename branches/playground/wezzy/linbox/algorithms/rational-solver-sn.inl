
/* Copyright (C) 2011 LinBox
 * Written Bryan Youse <>
 *
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

// rational-solver-sn.inl
// WARNING: this file is included from INSIDE the definition of class LinBox::RationalSolverSN

#ifdef DEBUGRC
#define debug(x) std::cerr << x << std::endl
#define debugneol(x) std::cerr << x
#define debug2(x, y) debug(x); y
#else
#define debug(x)
#define debug2(x, y)
#define debugneol(x)
#endif

#define log_2(x) (log(x)/M_LN2)

int zw_shift(NumericSolver & NS_S, size_t n, FVector &r, FVector &x)
{
	//  ZW method for calculating shift
	// compute ax
	FVector ax(n);
	NS_S.apply(ax, x);
	// compute ax = ax -r, the negative of residual
	for(size_t i=0; i<n; i++)
		field().sub(ax[i], ax[i], r[i]);
	// compute possible shift
	int zw_shift_loc;
	double normr1, normr2, normr3, shift1, shift2;
	normr1 = zw_dmax(n, &*r.begin(), 1);
	normr2 = zw_dmax(n, &*ax.begin(), 1);
	normr3 = zw_dmax(n, &*x.begin(), 1);

	//cerr << normr1 << " " << normr2 << " " << normr3 << endl;
	//try to find a good scalar
	if (normr2 <.0000000001)
		zw_shift_loc = 30;
	else {
		shift1 = floor(log_2 (normr1 / normr2)) - 2;
		zw_shift_loc = (int)(30 < shift1 ? 30 : shift1);
	}

	normr3 = normr3 > 2 ? normr3 : 2;
	shift2 = floor(53. / log_2 (normr3));
	zw_shift_loc = (int)(zw_shift_loc < shift2 ? zw_shift_loc : shift2);

	return zw_shift_loc;
}

template <class IMatrix>
int rat_sol(IVector& numx, Int& denx, FVector& xs_int, FVector& xs_frac, IVector &b, IVector &lastb, FVector& r, FVector& lastr, FVector& x, integer &loopBound, IMatrix &IM) {

	debug("Matrix norm: " << mnorm);
	int thresh_expt = 1;  // consider as low as 1 or 2 (at least when n large)
	double threshold = 1.0/(1 << thresh_expt);

	size_t n = r.size();
	// xs_int: integer part of scaled x.
	// xs_frac: fractional part of scaled x.
	// x * 2^shift = xs_int + xs_frac.
	FVector nextx(n), quo(n);

	integer denx_i;
	typename Field::Element one; field().init(one, 1);

	//  need to save original r for zw_shift calculation
	//  TODO: I took out the ZWSHIFT, still need last r??

	if(denx == 1){
		// compute first approximate solution x
		_numsolver.solve(x, r);

		//writeVec(x, "x from first solve");
		copy(r.begin(), r.end(), lastr.begin());
		copy(b.begin(), b.end(), lastb.begin());
		//-- compute xs_int, xs_frac, r (residue)
		update_xs(xs_int, xs_frac, x);
		if(exact_apply) update_r_exact(b, r, xs_int, IM);
		else update_r(r, xs_int);
	}

	//std::cerr << "while loopbound " << loopBound << std::endl;
	while (_ring.compare(denx, loopBound) < 0) {
		++iterations;
		//  x = DM^{-1}*r
		_numsolver.solve(nextx, r);

		for(size_t i=0; i<n; i++){
			// TODO - analyze logic here
			//  quo[i] = xs_frac[i] / nextx[i] - 1;
			//field().div(quo[i], xs_frac[i], nextx[i]);
			//field().subin(quo[i], one);
			field().sub(quo[i], xs_frac[i], nextx[i]);
		}

		double q = zw_dmax((int)n, &*quo.begin(), 1);

		//writeVec(nextx, "nextx from loop");  //  DEBUG PRINTOUTS
		//writeVec(xs_frac, "xs_frac from loop");
		//writeVec(quo, "quo");
		/* fails to write!!
		   _VDF.write(std::cerr << "nextx from loop: ", nextx) << std::endl;
		   _VDF.write(std::cerr << "xs_frac from loop: ", xs_frac) << std::endl;
		   _VDF.write(std::cerr << "quo minus 1: ", quo) << std::endl;
		   if ( 0 == iterations%500) std::cerr << iterations << " MAX DIFFERENCE: " << q << std::endl;
		   */

		//if (q == 0.0)   (QUIT HERE??)
		if (q < threshold) {
			HIT++;
			// update numx and denx
			denx <<= shift;
			update_num (numx, xs_int);

			if(_VDF.isZero(r))
				return 2;

			// consider increasing the shift for next iteration
			upshift();

			// make x = nextx, then compute new residual, xs_int, xs_frac
			swap(x, nextx);
			copy(r.begin(), r.end(), lastr.begin());
			copy(b.begin(), b.end(), lastb.begin());
			update_xs(xs_int, xs_frac, x);
			if(exact_apply) update_r_exact(b, r, xs_int, IM);
			else update_r(r, xs_int);
		}
		else {  //  q >= threshold, back-off
			MISS++;

			//  point of no return, quit
			if (shift < 2){
				/*
				   std::cerr << "rat_sol failure, no bits in x overlap in nextx." << std::endl;
				   std::cerr << "Iterations: " << iterations << std::endl;
				   writeVec(b, "b", 0, 5); writeVec(r, "r", 0, 5);
				   */
				return -1;
			}

			downshift();
			//-- compute xs_int, xs_frac, r (residue)
			// but use last r as input
			copy(lastr.begin(), lastr.end(), r.begin());
			copy(lastb.begin(), lastb.end(), b.begin());
			update_xs(xs_int, xs_frac, x);
			if(exact_apply) update_r_exact(b, r, xs_int, IM);
			else update_r(r, xs_int);
		}
	}
	return 0;
}// rat_sol

inline void upshift()
{
	switch(sstatus){
		//  exponential increase
	case SHIFT_GROW:
		shift_prev = shift;
		debugneol("G");
		shift_max = shift<<=1;
		break;
	case SHIFT_SEARCH:
		shift_prev = shift;
		debugneol("S");
		shift = (shift + shift_max)>>1;
		searchPeak = true;
		break;
	case SHIFT_PEAK:
		debugneol("P");
		//  maybe increase if we have been successful for a while
		break;
	case SHIFT_MAX:
		debugneol("M");
		//  machine precision-- can go no higher
		break;
	case SHIFT_SHRINK:
		debugneol("H");
		break;
	}
	if(shift > SHIFT_BOUND){
		shift_max = shift = SHIFT_BOUND;
		sstatus = SHIFT_MAX;
	}
	debug("^shift: "  << shift << "  max: " << shift_max << "  prev: " << shift_prev);
}

inline void downshift()
{
	/*
	   shift -= 2;
	   sstatus = SHIFT_PEAK;
	   cerr << "peaked at " << shift << endl;
	   */
	//  back up
	switch(sstatus){
	case SHIFT_GROW:
		debugneol("G");
	case SHIFT_MAX:
		debugneol("M");
	case SHIFT_SEARCH:
		//  TODO - previous shift could fail
		debugneol("S");
		shift_max = shift;
		shift = (shift_prev + shift)>>1;
		sstatus = SHIFT_SEARCH;
		//  searchPeak true means we were going up but got knocked back
		if(shift == shift_prev || searchPeak)
			sstatus = SHIFT_PEAK;
		break;
	case SHIFT_SHRINK:
		debugneol("H");
		shift >>= 1;
		break;
	case SHIFT_PEAK:
		debugneol("P");
		shift -= 1;
		break;
	default:
		break;
	}

	debug("vshift: "  << shift << "  max: " << shift_max << "  prev: " << shift_prev);
} // downshift

inline void update_xs(FVector& xs_int, FVector& xs_frac, FVector& x)
{
	Float scalar, tmp;
	int64_t shifted = ((int64_t)1 << shift);
	field().init(scalar, (double) shifted);

	//  make xs_int and xs_frac such that x*scalar = xs_int + xs_frac.
	for(size_t i = 0; i < xs_int.size(); ++i){
		// TODO: tmp can overflow a double
		tmp = x[i]*scalar;
		xs_int[i] = floor(tmp + 0.5);  // TODO: ceiling?
		// TODO: TRY THIS xs_int[i] = ceiling(tmp);
		xs_frac[i] = tmp - xs_int[i];
	}
	return;
}

inline void update_r(FVector& r, FVector& xs_int)
{
	Float scalar;
	size_t n = r.size();
	int64_t shifted = ((int64_t)1 << shift);
	field().init(scalar, (double)shifted);
	FVector y(n);

	//update r = r * 2^shift - Mat*xs_int
	_VDF.mulin(r, scalar);
	_numsolver.apply(y, xs_int);
	_VDF.subin(r, y);

	return;
} // update_r

template <class IMatrix>
inline void update_r_exact(IVector& r_exact, FVector& r, FVector& xs_int, IMatrix &IM){
	size_t n = r.size();

	IVector x_i(n), y_i(n);

	typename Ring::Element scalar = ((int64_t)1 << shift);

	// update r = r * 2^shift - Mat*xs_int
	//  r *= 2 ^shift
	_VDR.mulin(r_exact, scalar);

	// determine if exact apply is needed
	double vnorm = zw_dOOnorm(&*xs_int.begin(), (int)n, 1);

	int64_t th = ((int64_t)1 << 52);  // double mantissa
	Float thresh;
	field().init(thresh, (double)th);

	debugneol("vnorm " << vnorm);

	//  r -= Mat * xs_int
	if(field().mulin(vnorm, mnorm) < thresh){
		debugneol("Numeric ");
		FVector y(n);
		_numsolver.apply(y, xs_int);
		for(size_t i = 0; i < n; ++i)
			_ring.init(y_i[i], y[i]);
	}
	else{
		SHIFT_BOUND--; // to less this possibility
		debugneol("Exact ");
		for(size_t i = 0; i < n; ++i)
			_ring.init(x_i[i], xs_int[i]);
		IM.apply(y_i, x_i);
	}

	_VDR.subin(r_exact, y_i);

	//  convert exactly computed residue back to double (for next solve)
	typename FVector::iterator rp = r.begin();
	typename IVector::iterator rep = r_exact.begin();
	for(; rp!= r.end(); ++rp, ++rep)
		field().init(*rp, *rep);

	return;
} // update_r_exact

// no longer called...
inline int HadamardBound(integer& B, FMatrix& DM)
{
	size_t n = DM.rowdim();
	zw_hbound (B, n, n, &*DM.Begin()); // compute the Hadamard bound
	B = B * B;
	double mnorm_loc = zw_dOOnorm(&*DM.Begin(), n, n);

	// [don't know what this comment is about] should be a check for 2 * mnorm + zw_dmax (n, b, 1);
	// TODO what is "b"? from copied code it is the RHS array of doubles
	// zw_max just seems to get abs(max value of b)
	// next line false, just to compile
	double *b;
	B *= 2 * mnorm_loc + zw_dmax (n, b, 1); // [don't know what this factor is about]
	B <<= 1; // [extra factor of 2 for some reason... ]
	return B;
}

//update num, *num <- *num * 2^shift + d
inline IVector& update_num (IVector& num, const FVector& d)
{
	size_t n = d.size();
	IVector d_i(n);
	for (size_t i = 0; i < n; ++i) {
		_ring.init(d_i[i], d[i]);
	}
	Int scalar; _ring.init(scalar, 1UL << shift);
	//  TODO - analyze GMP shifting capability
	_VDR.mulin(num, scalar);
	_VDR.addin(num, d_i);
	return num;
}

//update r = r * shift - M d
inline static int update_r_ll (double* r, int n, const double* M, const double* d, int shift)
{
	long long int tmp;
	double* p1;
	const double* p2;
	const double* pd;
	for (p1 = r, p2 = M; p1 != r + n; ++ p1) {
		tmp = (long long int) *p1;
		tmp <<= shift;
		for (pd = d; pd != d + n; ++ pd, ++ p2) {
			tmp -= (long long int)*pd * (long long int) *p2;
		}
		*p1 = tmp;
	}
	return 0;
}

inline static size_t nextPower2(size_t n)
{
	size_t p = 1;
	while(p < n) p <<= 1;
	return p;
}

inline static double highAbs(FMatrix M)
{
	double max = 0;
	typename FMatrix::Iterator ri = M.Begin();
	for(; ri != M.End(); ++ri){
		double tmp = fabs(*ri);
		if(max < tmp) max = tmp;
	}
	return max;
}

inline static double zw_dOOnorm(const double* M, int m, int n)
{
	double norm = 0;
	double old = 0;
	const double* p;
	for (p = M; p != M + (m * n); ) {
		old = norm;
		norm = cblas_dasum (n, p ,1);
		if (norm < old) norm = old;
		p += n;
	}
	return norm;
}

inline static double zw_dmax (const int N, const double* a, const int inc)
{
	return fabs(a[cblas_idamax (N, a, inc)]);
}

inline static int zw_hbound (integer& b, int m, int n, const double* M)
{
	double norm = 0;
	const  double* p;
	integer tmp;
	b = 1;
	for (p = M; p != M + (m * n); ) {
		norm = cblas_dnrm2 (n, p ,1);
		tmp =  norm;
		integer::mulin (b, tmp);
		p += n;
	}

	return 0;
}

/*  out:  vector to print
 *  tag:  prepend this to output
 *  bound:  upper bound on printed value (0 - ignore bound) */
template<class Vec>
std::ostream& writeVec(Vec& out, const char *tag="", integer bound=/*40000000*/0,
		       size_t numEntries=5, std::ostream &os=std::cerr)
{
	os << tag << ": [";

	size_t n = (out.size() < numEntries ? out.size() : numEntries);
	for( size_t i=0; i<n; ++i){
		if(bound && ((integer)(out[i]) > bound || (integer)(out[i]) < -bound)){
			os << " entry over bound ]"  << std::endl;
			return os;
		}
		else
			os << " " << out[i];
	}
	if(out.size() > numEntries)
		os << " ...";
	os << " ]" << std::endl;
	return os;
}

template<class Vec>
void writeVecFile(Vec& out, const char* file)
{
	std::ofstream os;
	os.open(file, std::ios::out);

	typename Vec::const_iterator vi = out.begin();

	for( ; vi != out.end(); ++vi){
		os << *vi << std::endl;
	}
	os.close();
}

// to write diagnostic info to files for further testing
template <class Matrix>
void dumpData(const Matrix &M, const IVector &b, IVector &numx, integer &denx, integer &denBound)
{
#ifdef SPITOUT
	std::ofstream matout;
	matout.open("debug.mat", std::ios::out);
	M.write(matout);
	matout.close();
	writeVecFile(b, "debug.rhs");
	writeVecFile(numx, "debug.num");
	std::ofstream dout;
	dout.open("debug.den", std::ios::out);
	dout << denx << std::endl;
	dout << denBound << std::endl;
	dout << numx.size() << std::endl;
	dout.close();
#endif
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

