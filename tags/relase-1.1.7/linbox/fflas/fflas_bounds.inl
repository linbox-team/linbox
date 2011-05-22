
/* fflas/fflas_bounds.inl
 * Copyright (C) 2008 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifdef _LINBOX_LINBOX_CONFIG_H
#define FFLAS_INT_TYPE Integer
#else
#define FFLAS_INT_TYPE long unsigned int
#endif

/**
 * MatMulParameters
 *
 * \brief Computes the threshold parameters for the cascade
 *        Matmul algorithm
 *
 * 
 * \param F Finite Field/Ring of the computation.
 * \param k Common dimension of A and B, in the product A x B
 * \param bet Computing AB + beta C
 * \param delayedDim Returns the size of blocks that can be multiplied
 *                   over Z with no overflow
 * \param base Returns the type of BLAS representation to use
 * \param winoRecLevel Returns the number of recursion levels of
 *                     Strassen-Winograd's algorithm to perform
 * \param winoLevelProvided tells whether the user forced the number of
 *                          recursive level of Winograd's algorithm
 */
template <class Field>
inline void FFLAS::MatMulParameters (const Field& F,
				     const size_t k,
				     const typename Field::Element& beta,
				     size_t& delayedDim,
				     FFLAS_BASE& base,
				     size_t& winoRecLevel,
				     bool winoLevelProvided) {

	// Strategy : determine Winograd's recursion first, then choose appropriate
	// floating point representation, and finally the blocking dimension.
	// Can be improved for some cases.

	if (!winoLevelProvided)
		winoRecLevel = WinoSteps (k);
	base = BaseCompute (F, winoRecLevel);
	delayedDim = DotProdBound (F, winoRecLevel, beta, base);

	size_t n = k;
	size_t winoDel = winoRecLevel;

	// Computes the delayedDim, only depending on the recursive levels
	// that must be performed over Z
	while (winoDel > 0 && delayedDim < n) {
		winoDel--;
		delayedDim = DotProdBound (F, winoDel, beta, base);
		n >>= 1;
	}
	delayedDim = MIN (n, delayedDim);
}

/**
 * DotProdBound
 *
 * \brief  computes the maximal size for delaying the modular reduction
 *         in a dotproduct
 *
 * This is the default version assuming a conversion to a positive modular representation
 * 
 * \param F Finite Field/Ring of the computation
 * \param winoRecLevel Number of recusrive Strassen-Winograd levels (if any, 0 otherwise)
 * \param beta Computing AB + beta C
 * \param base Type of floating point representation for delayed modular computations
 * 
 */
template <class Field>
inline size_t FFLAS::DotProdBound (const Field& F,
				   const size_t w, 
				   const typename Field::Element& beta,
				   const FFLAS_BASE base) {
	
	FFLAS_INT_TYPE p;
	F.characteristic(p);
	typename Field::Element mone;
	F.init (mone, -1.0);

	unsigned long mantissa =
		(base == FflasDouble) ? DOUBLE_MANTISSA : FLOAT_MANTISSA;

	if (p == 0)
		return 1;
	
	double kmax;
	if (w > 0) {
		double c = computeFactor (F,w);
		double d = (double (1ULL << mantissa) /(c*c) + 1);
		if (d < 2)
			return 1;
		kmax = floor (d * (1ULL << w));
	} else {
		////// A fixer: (p-1)/2 si balanced

		double c = p-1;
		double cplt=0;
		if (!F.isZero (beta)){
			if (F.isOne (beta) || F.areEqual (beta, mone)) cplt = c;
			else cplt = c*c;
		}
		kmax = floor ( (double ((1ULL << mantissa) - cplt)) / (c*c));
		if (kmax  <= 1)
			return 1;
		}
	//kmax--; // we computed a strict upper bound
	return  (size_t) MIN (kmax, 1ULL << 31);
}

/**
 * Internal function for the bound computation
 * Generic implementation for positive representations
 */
template <class Field>
inline double FFLAS::computeFactor (const Field& F, const size_t w){
	FFLAS_INT_TYPE p;
	F.characteristic(p);
	size_t ex=1;
	for (size_t i=0; i < w; ++i) 	ex *= 3;
	return double(p - 1) * (1 + ex) / 2;
}

/**
 * WinoSteps
 *
 * \brief Computes the number of recursive levels to perform
 *
 * \param m the common dimension in the product AxB
 * 
 */
inline size_t FFLAS::WinoSteps (const size_t m) {
	size_t w = 0;
	size_t mt = m;
	while (mt >= WINOTHRESHOLD) {w++; mt >>= 1;}
	return w;
}

/**
 * BaseCompute
 *
 * \brief Determines the type of floating point representation to convert to,
 *        for BLAS computations
 * \param F Finite Field/Ring of the computation
 * \param w Number of recursive levels in Winograd's algorithm
 * 
 */
template <class Field>
inline FFLAS::FFLAS_BASE FFLAS::BaseCompute (const Field& F, const size_t w){
	
	FFLAS_INT_TYPE pi;
	F.characteristic(pi);
	FFLAS_BASE base;
	switch (w) {
	case 0: base = (pi < FLOAT_DOUBLE_THRESHOLD_0)? FflasFloat : FflasDouble;
		break;
	case 1:  base = (pi < FLOAT_DOUBLE_THRESHOLD_1)? FflasFloat : FflasDouble;
		break;
	case 2:  base = (pi < FLOAT_DOUBLE_THRESHOLD_2)? FflasFloat : FflasDouble;
		break;
	default: base = FflasDouble;
		break;
	}
	return base;
}


/*************************************************************************************
 * Specializations for ModularPositive and ModularBalanced over double and float
 *************************************************************************************/

template <class Element>
inline double computeFactor (const ModularBalanced<Element>& F, const size_t w){
	FFLAS_INT_TYPE p;
	F.characteristic(p);
	size_t ex=1;
	for (size_t i=0; i < w; ++i) 	ex *= 3;
	return  (p - 1) * ex / 2; 
}

template <>
inline FFLAS::FFLAS_BASE FFLAS::BaseCompute (const Modular<double>& ,
					     const size_t ){
	return FflasDouble;
}

template <>
inline FFLAS::FFLAS_BASE FFLAS::BaseCompute (const Modular<float>& ,
					     const size_t ){
	return FflasFloat;
}

template <>
inline FFLAS::FFLAS_BASE FFLAS::BaseCompute (const ModularBalanced<double>& ,
					     const size_t ){
	return FflasDouble;
}

template <>
inline FFLAS::FFLAS_BASE FFLAS::BaseCompute (const ModularBalanced<float>& ,
					     const size_t ){
	return FflasFloat;
}

/**
 * TRSMBound
 *
 * \brief  computes the maximal size for delaying the modular reduction
 *         in a triangular system resolution
 *
 * This is the default version over an arbitrary field.
 * It is currently never used (the recursive algorithm is run until n=1 in this case)
 * 
 * \param F Finite Field/Ring of the computation
 * 
 */
template <class Field>
inline size_t FFLAS::TRSMBound (const Field& F) {
	return 1;	
}

/**
 * Specialization for positive modular representation over double
 * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^53
 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
 */
template<>
inline size_t FFLAS::TRSMBound (const Modular<double>& F){

	FFLAS_INT_TYPE pi;
	F.characteristic(pi);
	unsigned long p(pi);
        unsigned long long p1(1UL), p2(1UL);
	size_t nmax = 0;
	unsigned long long max = ( (1ULL << (DOUBLE_MANTISSA + 1) ) / ((unsigned long long)(p - 1)));
	while ( (p1 + p2) < max ){
		p1*=p;
		p2*=p-2;
		nmax++;
	}
	return nmax;
}


/**
 * Specialization for positive modular representation over float
 * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^24
 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
 */
template<>
inline size_t FFLAS::TRSMBound (const Modular<float>& F){

	FFLAS_INT_TYPE pi;
	F.characteristic(pi);
	unsigned long p(pi);
        unsigned long long p1(1UL), p2(1UL);
	size_t nmax = 0;
	unsigned long long max = ( (1ULL << (FLOAT_MANTISSA + 1) ) / ((unsigned long long)(p - 1)));
	while ( (p1 + p2) < max ){
		p1*=p;
		p2*=p-2;
		nmax++;
	}
	return nmax;
}

/**
 * Specialization for balanced modular representation over double
 * Computes nmax s.t. (p-1)/2*(((p+1)/2)^{nmax-1}) < 2^53
 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
 */
template<>
inline size_t FFLAS::TRSMBound (const ModularBalanced<double>& F){

	FFLAS_INT_TYPE pi;
	F.characteristic (pi);
	unsigned long p= (pi + 1) / 2;
        unsigned long long p1(1UL);
	size_t nmax = 0;
	unsigned long long max = ((1ULL << (DOUBLE_MANTISSA + 1)) / ((unsigned long long)(p - 1)));
	while (p1 < max){
		p1 *= p;
		nmax++;
	}
	return nmax;
}

/**
 * Specialization for balanced modular representation over float
 * Computes nmax s.t. (p-1)/2*(((p+1)/2)^{nmax-1}) < 2^24
 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
 */
template<>
inline size_t FFLAS::TRSMBound (const ModularBalanced<float>& F){

	FFLAS_INT_TYPE pi;
	F.characteristic (pi);
	unsigned long p = (pi + 1) / 2;
        unsigned long long p1(1UL);
	size_t nmax = 0;
	unsigned long long max = ((1ULL << (FLOAT_MANTISSA + 1)) / ((unsigned long long) (p - 1)));
	while (p1 < max){
		p1 *= p;
		nmax++;
	}
	return nmax;

}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
