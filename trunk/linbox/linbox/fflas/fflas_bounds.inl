/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* fflas/fflas_bounds.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifdef _LINBOX_CONFIG_H
#define FFLAS_INT_TYPE Integer
#else
#define FFLAS_INT_TYPE long unsigned int
#endif

//---------------------------------------------------------------------
// DotProdBound :
// Computes the maximal dimension k so that the product A*B+beta.C over Z,
// where A is m*k and B is k*n can be performed correctly with w Winograd
// recursion levels on the 53 bits of double mantissa
//---------------------------------------------------------------------
template  <class Field> 
inline size_t DotProdBoundCompute (const Field& F, const size_t w, 
				   const typename Field::Element& beta)
{
	typename Field::Element mone;
	static FFLAS_INT_TYPE p;
	F.characteristic(p);
	F.init (mone, -1.0);
	size_t kmax;
	if (p == 0)
		kmax = 2;
	else
		if (w > 0) {
			size_t ex=1;
			for (size_t i=0; i < w; ++i) 	ex *= 3;
			//FFLAS_INT_TYPE c = (p-1)*(ex)/2; //bound for a centered representation
			FFLAS_INT_TYPE c = (p-1)*(1+ex)/2;
			kmax =  lround(( double(1ULL << DOUBLE_MANTISSA) /double(c*c) + 1)*(1 << w));
			if (kmax ==  ( 1ULL << w))
				kmax = 2;
		}
		else{
			FFLAS_INT_TYPE c = p-1;
			FFLAS_INT_TYPE cplt=0;
			if (!F.isZero (beta))
				if (F.isOne (beta) || F.areEqual (beta, mone))
					cplt = c;
				else cplt = c*c;
			kmax =  lround(( double((1ULL << DOUBLE_MANTISSA) - cplt)) /double(c*c));
			if (kmax  < 2)
				kmax = 2;
		}
	return  MIN(kmax,1ULL<<31);
}
template  <> 
inline size_t DotProdBoundCompute (const Modular<double>& F, const size_t w, 
				   const double& beta)
{
	double mone;
	static FFLAS_INT_TYPE p;
	F.characteristic(p);
	F.init (mone, -1.0);
	size_t  kmax;
	if (p == 0)
		kmax = 2;
	else
		if (w > 0) {
			size_t ex=1;
			for (size_t i=0; i < w; ++i) 	ex *= 3;
			long long c;
#ifndef _LINBOX_CONFIG_H
			if (F.balanced)
 				c = (p-1)*(ex)/2; // balanced representation
 			else
#endif
			c = (p-1)*(1+ex)/2; // positive representation
			kmax =  lround(double(1ULL << DOUBLE_MANTISSA) /(double(c*c) + 1)*(1ULL << w));
			if (kmax ==  ( 1ULL << w))
				kmax = 2;
		}
		else{
			FFLAS_INT_TYPE c = p-1;
			FFLAS_INT_TYPE cplt=0;
			if (!F.isZero (beta))
				if (F.isOne (beta) || F.areEqual (beta, mone))
					cplt = c;
				else cplt = c*c;
			kmax =  lround( double((1ULL << DOUBLE_MANTISSA) - cplt) /(double(c*c)));
			if (kmax  < 2)
				kmax = 2;
		}
	return (size_t) MIN(kmax,1ULL<<31);
}

template  < class Field > 
inline size_t FFLAS::DotProdBound (const Field& F, const size_t w, 
				   const typename Field::Element& beta)
{
	static Field G = F;
	static FFLAS_INT_TYPE pig;
	G.characteristic(pig);
	FFLAS_INT_TYPE pif;
	F.characteristic(pif);
	static typename Field::Element b = beta;
	static size_t w2 = w;
	static size_t kmax = DotProdBoundCompute (F, w, beta);
     	if ((b != beta) || (pif != pig) ||  (w2 != w)) {
		G = F;
		w2 = w;
		b = beta;
		kmax =  DotProdBoundCompute (F, w, beta);
	}	
	return kmax;
}

//---------------------------------------------------------------------
// TRSMBound
// Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^53
//---------------------------------------------------------------------
size_t bound_compute(const long long pi) {
	
	long long p=pi,p1=1,p2=1;
	size_t nmax=1;
	double max = ( (  1ULL<<(DOUBLE_MANTISSA+1) )/(p-1));
	while ( (p1 + p2) < max ){
		p1*=p;
		p2*=p-2;
		nmax++;
	}
	nmax--;
	return nmax;
}
size_t bound_compute_balanced(const long long pi) {
	
	long long p=(pi+1)/2,p1=1;
	size_t nmax=0;
	double max = ( (  1ULL<<(DOUBLE_MANTISSA))/(p-1));
	while ( (p1) < max ){
		p1*=p;
		nmax++;
	}
	return nmax;
}

template<class Field>
size_t FFLAS::TRSMBound (const Field& F) {
	FFLAS_INT_TYPE pi;
	F.characteristic(pi);
	static long unsigned int p=pi;
	static size_t nmax=bound_compute(p);
	if (p == pi) 
		return nmax;
	else 
		return nmax=bound_compute(p=pi);
}

template<>
size_t FFLAS::TRSMBound (const Modular<double>& F) {
	FFLAS_INT_TYPE pi;
	F.characteristic(pi);
	static FFLAS_INT_TYPE p=pi;
#ifdef _LINBOX_CONFIG_H
	static size_t nmax = bound_compute(pi);
#else
	static size_t nmax = (F.balanced) ? bound_compute_balanced(pi) : bound_compute(pi);
#endif	
	if (p == pi) 
		return nmax;
	else
#ifdef _LINBOX_CONFIG_H
		return nmax= bound_compute (p=pi); //(F.balanced) ? bound_compute_balanced(p=pi) : bound_compute(p=pi);
#else
        	return (F.balanced) ? bound_compute_balanced(p=pi) : bound_compute(p=pi);
#endif
}

