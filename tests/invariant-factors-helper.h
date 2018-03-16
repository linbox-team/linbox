class TestInvariantFactorsHelper {
public:
	const size_t _p;
	const Field F;
	const PolynomialRing R;
	const MatrixDomain<Field> MD;
	
	TestInvariantFactorsHelper(size_t p) : _p(p), F(p), R(p), MD(F) {};
	
	std::ostream & writeInvariantFactors(std::ostream &os, std::vector<Polynomial> &factors, bool full=true) const {
		for (size_t i = 0; i < factors.size(); i++) {
			Polynomial f;
			R.monic(f, factors[i]);
			if (full)
				R.write(os, f) << std::endl;
			else
				if (not R.isOne(f))
				R.getCoeffField().write(os << R.deg(f) << ", const ", f[0]) << std::endl;
		}
		return os;
	}

	std::ostream & writeInvariantFactorDegrees(std::ostream &os, std::vector<Polynomial> &factors, bool full=true) const {
		os << R.deg(factors[0]); 
		for (size_t i = 0; i < factors.size(); i++) 
			os << ", " << R.deg(factors[i]);
		return os;
	}

	// count nontrivial (not 1 or x) invariant factors.
	size_t countInvariantFactors(std::vector<Polynomial> &factors) const {
		Polynomial x; 
		std::vector<long> coeffs; coeffs.push_back(0); coeffs.push_back(1);
		R.init(x, coeffs);

		size_t count = 0;
		for (size_t i = 0; i < factors.size(); i++) {
			Polynomial & p = factors[i];
			if (not (R.deg(p) < 1) and 
				not (R.deg(p) == 1 and R.getCoeffField().isZero(p[0])))
				count++;
		}
	//	if (count >= 5) //R.write(std::cout << "long ", factors[13]) << ", " << count << std::endl;
		//writeInvariantFactorDegrees(std::cout << "long ", factors) << std::endl;
		return count;
	}
	
	// total degree of invariant factors.
	size_t degInvariantFactors(std::vector<Polynomial> &factors) const {
		size_t d = 0;
		for (size_t i = 0; i < factors.size(); i++) 
			d += R.deg(factors[i]);
		return d;
	}
	
	double computeMinpoly(std::vector<size_t> &degree, std::vector<Matrix> &minpoly, const ProductMat &M, size_t b) const {
		size_t n = M.rowdim();
		
		RandIter RI(F);
		RandomMatrix RM(F, RI);
		
		Matrix U(F, b, n);
		Matrix V(F, n, b);
		
		RM.random(U);
		RM.random(V);
		
		// Construct block sequence to input to BM
		Sequence seq(&M, F, U, V);
		
		// Compute minimal generating polynomial matrix
		MasseyDom BMD(&seq); // pascal
		CoppersmithDom BCD(MD, &seq, 10); // george
				
		TW.clear();
		TW.start();
		
		//BMD.left_minpoly_rec(minpoly, degree);
		degree = BCD.right_minpoly(minpoly);
		
		TW.stop();
		double bm_time = TW.usertime();
		std::cout << "seq->gen " << bm_time << " " << std::flush;
		
		return bm_time;
	}
	
	double convertMinPolyToPolyMatrix(PolyMatrix &G, const std::vector<Matrix> &minpoly) {
		size_t b = G.rowdim();
		
		TW.clear();
		TW.start();
		
		for (size_t i = 0; i < b; i++) {
			for (size_t j = 0; j < b; j++) {
				std::vector<long> coeffs;
				for (size_t k = 0; k < minpoly.size(); k++) {
					long coeff;
					F.convert(coeff, minpoly[k].getEntry(i, j));
					coeffs.push_back(coeff);
				}
				Polynomial tmp;
				R.init(tmp, coeffs);
				
				G.setEntry(i, j, tmp);
			}
		}
		
		TW.stop();
		double cv_time = TW.usertime();
		// std::cout << cv_time << " " << std::flush;
		
		return cv_time;
	}
	
	void computeDet(Polynomial &det, std::vector<Polynomial> &factors) {
		R.assign(det, R.one);
		for (size_t i = 0; i < factors.size(); i++) {
			R.mulin(det, factors[i]);
		}
		
		Polynomial monic_det;
		R.monic(monic_det, det);
		R.assign(det, monic_det);
	}
	
	void timeTextbook(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveTextBook(result, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << ", text_time " << sf_time << " " << std::flush;
	}
	
	double timeKannanBachem(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		bool mem_error = false;
		try {
			SFD.solve(result, G);
		} catch (const std::bad_alloc &e) {
			mem_error = true;
		}
		
		TW.stop();
		double sf_time = TW.usertime();
		
		if (!mem_error) {
			std::cout << ", kb_time " << sf_time << " " << std::flush;
		} else {
			std::cout << "MEM " << std::flush;
		}
		
		return mem_error ? -1 : sf_time;
	}
	
	void timeHybrid(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveAdaptive(result, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << ", hybrid " << sf_time << " " << std::flush;
	}
	
	void timeHybrid2(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveAdaptive2(result, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << ", hybrid2 " << sf_time << " " << std::flush;
	}
	
	double timeIliopoulos(std::vector<Polynomial> &result, const PolyMatrix &M, const Polynomial &det) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveIliopoulos(result, G, det);
		
		TW.stop();
		double sf_time = TW.usertime();
		//std::cout << sf_time << " " << std::flush;
		
		return sf_time;
	}
	
	void local(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &f,
		long multiplicity) const {
	
		SmithFormLocalDom SFD;
		
		Polynomial modulus;
		R.pow(modulus, f, multiplicity);
		
		QuotientRing QR(_p, modulus);
		
		QuotMatrix QM(M, QR);
		
		std::list<QPolynomial> L;
		SFD(L, QM, QR);
		
		Hom<PolynomialRing, QuotientRing> hom(R, QR);
		
		size_t j = 0;
		std::list<QPolynomial>::const_iterator it;
		for (it = L.begin(); it != L.end(); it++) {
			Polynomial tmp;
			hom.preimage(tmp, *it);
								
			if (R.isOne(tmp)) {
				// noop
			} else if (R.isZero(tmp)) {
				R.mulin(result[j], modulus);
			} else {
				R.mulin(result[j], tmp);
			}
			j++;
		}
	}
	
	void localFactored(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &sf_factor,
		long multiplicity) const {
	
		std::vector<std::pair<Polynomial, long>> factors;
		R.factor(factors, sf_factor);
		
		for (size_t i = 0; i < factors.size(); i++) {
			local(result, M, factors[i].first, factors[i].second * multiplicity);
		}
	}
	
	double timeFactoredLocal(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) {	
	
		std::vector<std::pair<Polynomial, long>> factors;
		
		TW.clear();
		TW.start();
		
		R.squareFree(factors, det);
		
		result.clear();
		for (size_t i = 0; i < M.rowdim(); i++) {
			result.push_back(R.one);
		}
		
		for (size_t i = 0; i < factors.size(); i++) {
			if (factors[i].second == 1) {
				R.mulin(result[result.size() - 1], factors[i].first);
			} else {
				localFactored(result, M, factors[i].first, factors[i].second);
			}
		}
		
		TW.stop();
		double fp_time = TW.usertime();
		std::cout << fp_time << " " << std::flush;
		
		return fp_time;
	}
	
	size_t detLimit(const PolyMatrix &M, size_t dim) {
		size_t limit1 = 0;
		for (size_t i = 0; i < M.rowdim(); i++) {
			size_t max_degree = 0;
			for (size_t j = 0; j < M.coldim(); j++) {
				size_t deg = R.deg(M.getEntry(i, j));
				if (deg > max_degree) {
					max_degree = deg;
				}
			}
			
			limit1 += max_degree;
		}
		
		size_t limit2 = 0;
		for (size_t i = 0; i < M.coldim(); i++) {
			size_t max_degree = 0;
			for (size_t j = 0; j < M.rowdim(); j++) {
				size_t deg = R.deg(M.getEntry(j, i));
				if (deg > max_degree) {
					max_degree = deg;
				}
			}
			
			limit2 += max_degree;
		}
				
		return std::min(std::min(limit1, limit2), dim) + 1;
	}
	
	double timeLocalX(Polynomial &det, const PolyMatrix &M, size_t exponent) {
		LocalSmithFormDom SFD(R, exponent);
		LocalRing L(R, exponent);
		
		std::vector<Polynomial> result;
		result.clear();
		
		LocalMatrix G(M, L);
		
		TW.clear();
		TW.start();
		
		SFD.solveDet(det, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		//std::cout << ", local " << sf_time << " " << std::flush;
		
		NTL::MakeMonic(det);
		
		return sf_time;
	}
}; // End of TestInvariantFactorsHelper

// preconditioners

void identity(SparseMat & P)
{
	for (size_t i=0; i < P.rowdim(); ++i)
		P.setEntry(i,i,P.field().one);
	P.finalize();
}

void identityVec(std::vector<size_t> & V, size_t n) 
{
	for (size_t i=0; i < n; ++i)
		V[i]=i;
}

void randomVec(std::vector<size_t> & V, size_t n) 
{
	V.resize(n);
	for (size_t i = 0; i < n; ++i) V[i] = i;
	// Knuth construction
	for (size_t i = 0; i < n-1; ++i) {
		size_t j = i + rand()%(n-i);
		std::swap(V[i], V[j]);
	}
}

void randomCycleVec (std::vector<size_t> & V, size_t n) 
{
	std::vector<size_t> Perm;
	randomVec(Perm, n);
//	for (size_t i = 0; i < Perm.size(); ++i) std::cout << Perm[i] << " ";
//	std::cout << std::endl;
	size_t first = Perm.back(); Perm.pop_back();
	size_t i = first;
	while (Perm.size() > 0) {
		size_t j = Perm.back(); Perm.pop_back();
		V[i] = j;
		i = j;
	}
	V[i] = first;
} 

void randomCycle(SparseMat & P)
{
	size_t n = P.rowdim();
	std::vector<size_t> Perm;
	randomVec(Perm, n);
//	for (size_t i = 0; i < Perm.size(); ++i) std::cout << Perm[i] << " ";
//	std::cout << std::endl;
	size_t first = Perm.back(); Perm.pop_back();
	size_t i = first;
	while (Perm.size() > 0) {
		size_t j = Perm.back(); Perm.pop_back();
		P.setEntry(i,j, P.field().one);
		i = j;
	}
	P.setEntry(i,first, P.field().one);
	P.finalize();
} 

void randomPerm(SparseMat & P)
{
	size_t n = P.rowdim();
	std::vector<size_t> Perm;
	randomVec(Perm, n);
	for (size_t i = 0; i < n; ++i) P.setEntry(i, Perm[i], P.field().one);
	P.finalize();
}

void fillBlock(SparseMat & P, size_t s, size_t b)
{
	const Field& F = P.field();
	BlasMatrix<Field> B(F,b,b);
	Field::Element x; F.init(x);
//	SparseMat B(P.field(),b,b);

	Field::Element d; F.init(d);
	do { // hack because BlasMatrix::random() is broken 
		for (size_t i = 0; i < b; ++i)
			for (size_t j = 0; j < b; ++j) {
				F.init(x,rand());
				B.setEntry(i,j,x);
			}
	} while (F.isZero(det(d,B))); 

	for (size_t i = 0; i < b; ++i)
		for (size_t j = 0; j < b; ++j)
			P.setEntry(i+s,j+s,B.getEntry(x,i,j));
}

void blockDiag(SparseMat & P, size_t b=5)
{
	BlasMatrix<Field> B(P.field(), b, b);
	size_t n = P.rowdim(); size_t q = n/b; size_t r = n - q*b;
	for (size_t k = 0; k < q-1; ++k)
		fillBlock(P,b*k,b);
	fillBlock(P,b*(q-1),b+r);
	P.finalize();
}

template<class Sp>
void permuteRows(Sp & MP, std::vector<size_t> & P, Sp& M) {
	assert(MP.rowdim() == M.rowdim());
	assert(P.size() == M.rowdim());

	size_t n = M.rowdim();
	for (size_t i = 0; i< n; ++i) 
		if (P[i] >= n) std::cout << "whoa" << std::endl;
	
	typename Sp::IndexedIterator p;
	for (p = M.IndexedBegin(); p != M.IndexedEnd(); ++p)
			MP.setEntry(P[p.rowIndex()],p.colIndex(),p.value());
	/*
	Element x; M.field().init(x);
	for (size_t i = 0; i < M.rowdim(); ++i)
	for (size_t j = 0; j < M.coldim(); ++j)
		if (not M.field().isZero(M.getEntry(x,i,j)))
		*/
	MP.finalize();
}

template<class Sp>
void permutePlusId(Sp & MP, std::vector<size_t> & P, Sp& M) {
	assert(MP.rowdim() == M.rowdim());
	assert(P.size() == M.rowdim());

	/*
	size_t n = M.rowdim();
	for (size_t i = 0; i< n; ++i) 
		if (P[i] >= n) std::cout << "whoa" << std::endl;
	*/
	
	Element x; M.field().init(x);
	Element y; M.field().init(y);
	typename Sp::IndexedIterator p;
	for (p = M.IndexedBegin(); p != M.IndexedEnd(); ++p) {
		size_t i = p.rowIndex();
		size_t j = p.colIndex();
		size_t k = P[i];
		// identity apply
		MP.getEntry(x,i,j);
		M.field().addin(x,p.value());
		MP.setEntry(i,j,x);  
		// permutation apply
		MP.getEntry(x,k,j);
		M.field().addin(x,p.value());
		MP.setEntry(k,j,x);  
	}
	/*
	Element x; M.field().init(x);
	for (size_t i = 0; i < M.rowdim(); ++i)
	for (size_t j = 0; j < M.coldim(); ++j)
		if (not M.field().isZero(M.getEntry(x,i,j)))
		*/
	MP.finalize();
}
#include "linbox/vector/vector-domain.h"

template<class Sp>
void mul(Sp & C, Sp & A, Sp& B) {
	const typename Sp::Field& F  = A.field();
	LinBox::VectorDomain<typename Sp::Field> VD(A.field());
	size_t n = A.rowdim();
	std::vector<typename Sp::Field::Element> x, y, z;
	typename Sp::Field::Element e;
	x.resize(n, F.zero); y.resize(n, F.zero); z.resize(n, F.zero);
	for (size_t i = 0; i < n; ++i) {
		F.assign(x[i], F.one );
		B.apply(y,x);
		A.apply(z,y);
		for (size_t j = 0; j < n; ++j) {
			if (not F.isZero(z[j])) 
				C.getEntry(e, i, j);
				F.addin(e, z[j]);
				C.setEntry(i,j,e);
		}
		F.assign(x[i], F.zero );
	}
	C.finalize();
}
