#include <linbox/linbox-config.h>

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

// NTL for random irreducible polynomial generation
#include "linbox/ring/ntl.h"

// Matrix Domains
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

#ifndef __LINBOX_SPARSE_MATRIX_GENERATOR_H
#define __LINBOX_SPARSE_MATRIX_GENERATOR_H

namespace LinBox
{
	template<class Field, class PolynomialRing>
	class SparseMatrixGenerator
	{
	public:
		typedef typename Field::Element Element;
		
		typedef typename PolynomialRing::Element Polynomial;
		typedef typename PolynomialRing::RandIter RandIter;
		typedef typename PolynomialRing::CoeffField CoeffField;
		typedef typename PolynomialRing::Coeff Coeff;
		typedef typename CoeffField::RandIter CoeffRandIter;
		
		typedef MatrixDomain<Field> MatrixDom;
		typedef typename MatrixDom::Matrix SubMatrix;
		
	private:
		Field _F;
		
		PolynomialRing _R;
		RandIter _RI;
		CoeffRandIter _CRI;
		
	public:
		SparseMatrixGenerator(const Field &F, const PolynomialRing &R): _F(F), _R(R), _RI(R), _CRI(R.getCoeffField()) {}
		
		void linearPolynomial(Polynomial &p, const Coeff& a) const {
			const CoeffField& F=_R.getCoeffField();
			Coeff ma; F.init(ma); 
			F.neg(ma, a);
			_R.init(p);
			_R.setCoeff(p, 1, F.one);
			_R.setCoeff(p, 0, ma);
		}

		void randomPolynomial(Polynomial &p, size_t d) const {
			_RI.random(p, d);
		}
		
		void randomIrreducible(Polynomial &p, size_t d) const {
			_RI.randomIrreducible(p, d);
		}
		
		void randomTrinomial(Polynomial &p, size_t d) const {
			if (d < 2) {
				std::cout << "Error: trinomial must have degree > 1" << std::endl;
				exit(1);
			}
			
			_R.init(p);
			_R.setCoeff(p, d, _R.getCoeffField().one);
			
			Coeff c;
			_CRI.random(c);
			_R.setCoeff(p, 0, c);
			
			size_t i = rand() % d;
			_CRI.random(c);
			_R.setCoeff(p, i, c);
		}
		
		void randomIrreducibleTrinomial(Polynomial &p, size_t d) const {
			do {
				randomTrinomial(p, d);
			} while (!_R.isIrreducible(p));
		}
		
		void readDivisor(std::ifstream &in, Polynomial &divisor) const {
			std::string type;
			in >> type;
						
			if (type == "c") {
				_R.read(in, divisor);
			} else if (type == "d") {
				size_t d;
				in >> d;
				
				randomPolynomial(divisor, d);
			} else if (type == "i") {
				size_t d;
				in >> d;
				
				randomIrreducible(divisor, d);
			} else if (type == "t") {
				size_t d;
				in >> d;
				
				randomTrinomial(divisor, d);
			} else if (type == "r") {
				size_t d;
				in >> d;
				
				randomIrreducibleTrinomial(divisor, d);
			} else {
				std::cout << "Error: unknown divisor type (" << type << ")" << std::endl;
				exit(1);
			}
		}
		
		/**
		 * Reads a file of format:
		 * <number of bumps>
		 * <multiplicity 1> <bump type 1> <bump 1>
		 * <multiplicity 2> <bump type 2> <bump 2>
		 * ...
		 * <multiplicity n> <bump type n> <bump n>
		 * 
		 * divisor types:
		 * c - const (format of polynomial ring read)
		 * d - random polynomial of degree d (use: d 9)
		 * i - random irreducible of degree d (use: i 4)
		 * t - random trinomial of degree d (use: t 5)
		 * r - random irreducible trinomial of degree d (use: r 3)
		 */
		void readFile(std::vector<Polynomial> &fs, const std::string &filename) const {			
			std::ifstream in(filename);
			size_t nbumps;
			in >> nbumps;
			
			Polynomial f;
			_R.assign(f, _R.one);
			
			fs.resize(0);
			
			for (size_t i = 0; i < nbumps; i++) {
				size_t multiplicity;
				in >> multiplicity;
				
				Polynomial bump;
				readDivisor(in, bump);
				_R.write(std::cout << "divisor (" << multiplicity << " times): ", bump) << std::endl;
				
				for (size_t j = 0; j < multiplicity; j++) {
					//Polynomial tmp;
					//_R.assign(tmp, bump);
					fs.push_back(bump);
				}
			}
		}
		
		template<class Matrix>
		bool makeCompanion(Matrix &M, const Polynomial &f, size_t offset) const {
			size_t dim = _R.deg(f);
			
			for (size_t i = 1; i < dim; i++) {
				M.setEntry(offset + i, offset + i - 1, _F.one);
			}
			
			Polynomial monic_f;
			_R.monic(monic_f, f);
			// _R.write(std::cout << "monic f: ", monic_f) << std::endl;
			
			std::vector<integer> coeffs;
			_R.convert(coeffs, monic_f);
			
			for (size_t i = 0; i < coeffs.size() - 1; i++) {
				Element tmp;
				_F.init(tmp, coeffs[i]);
				_F.negin(tmp);
				
				if (!_F.isZero(tmp)) {
					M.setEntry(offset + i, offset + dim - 1, tmp);
				}
			}
			
			return true;
		}
		
		/**
		 * Populates matrix with frobenius form with blocks corresponding to
		 * F_k = C_{f_1 * f_2 * ... * f_k}
		 * return: True if matrix successfully generated.
		 */
		template<class Matrix>
		bool build(Matrix &M, const std::vector<Polynomial> &fs) const {
			size_t min_dim = 0;
			for (size_t i = 0; i < fs.size(); i++) {
				min_dim += _R.deg(fs[i]);
			}
			
			if (M.rowdim() < min_dim || M.coldim() < min_dim) {
				M.resize(min_dim, min_dim);
			}
			
			size_t offset = 0;
			Polynomial f;
			for (size_t i = 0; i < fs.size(); i++) {				
				_R.assign(f, fs[i]);
				size_t d = _R.deg(f);
				
				makeCompanion(M, f, offset);
				
				offset += d;
			}
			
			M.finalize();
			
			return true;
		}
		
		// specialization for format SMM (sparse-map-map)	
		void fillIn(SparseMatrix<Field, SparseMatrixFormat::SMM> &M, 
					double targetSparsity) const {
			M.randomSim(int(M.rowdim()*M.coldim()*targetSparsity));
		}
		
		template<class Matrix1, class Matrix2>
		void copy(Matrix1 &M, Matrix2 &T) const {
			if (M.rowdim() != T.rowdim() || M.coldim() != T.coldim()) {
				M.resize(T.rowdim(), T.coldim());
			}
			
			for (size_t i = 0; i < M.rowdim(); i++) {
				for (size_t j = 0; j < M.coldim(); j++) {
					Element tmp;
					T.getEntry(tmp, i, j);
					
					if (_F.isZero(tmp)) {
						continue;
					}
					
					M.setEntry(i, j, tmp);
				}
			}
			
			M.finalize();
		}
		
	public:
		template<class Matrix>
		double nnz(const Matrix &M) const {
			double nnz = 0;
			
			for (size_t i = 0; i < M.rowdim(); i++) {
				for (size_t j = 0; j < M.coldim(); j++) {
					if (_F.isZero(M.getEntry(i,j))) {
						continue;
					}
					
					nnz++;
				}
			}
			
			return nnz;
		}
		
		template<class Matrix>
		double sparsity(const Matrix &M) const {
			return nnz(M) / (M.rowdim() * M.coldim());
		}
		
		template<class Matrix>
		void generate(Matrix &M, Polynomial &det, const std::vector<Polynomial> & fs, double sparsity) const {

			_R.assign(det, _R.one);
			for (size_t i = 0; i < fs.size(); i++) {
				_R.mulin(det, fs[i]);
			}
			if (_R.deg(det) < M.rowdim()) {
				Polynomial x, xm;
				std::vector<long> coeffs = {0, 1};
				_R.init(x, coeffs);
				_R.pow(xm, x, M.rowdim() - _R.deg(det));
				
				_R.mulin(det, xm);
			}
			
			// Use SMM format to generate the matrix
			// fill in is much faster that other formats
			SparseMatrix<Field, SparseMatrixFormat::SMM> T(_F, M.rowdim(), M.coldim());
			build(T, fs);
			fillIn(T, sparsity);
			
			copy(M, T);
		} // generate from vec

		template<class Matrix>
		void generate(Matrix &M, Polynomial &det, const std::string &filename, double sparsity) const {
			std::vector<Polynomial> fs;
			readFile(fs, filename);
			generate(M, det, fs, sparsity);
		} // generate from file

		std::vector<Polynomial>& augment(std::vector<Polynomial>& fs, size_t k, Polynomial& p) {
			for (size_t i = 0; i < k; ++i) fs.push_back(p);
			return fs;
		}

		std::vector<Polynomial>& invariants(std::vector<Polynomial>& fs, size_t n, int fsnum) {
			//choose within a small set of fs build schemes
			Polynomial p; _R.init(p);
			if (fsnum >= 0) {
				Polynomial xm1, xm1s, x2m1;
				_R.init(xm1); _R.init(xm1s); _R.init(x2m1); _R.init(p);
				linearPolynomial(xm1,_R.getCoeffField().one);
				linearPolynomial(x2m1,_R.getCoeffField().mOne); _R.mulin(x2m1, xm1);
				_R.mul(xm1s, xm1, xm1);
				switch (fsnum) {
				case 0: // but don't use
				case 1: //identity n (x-1)'s 
					augment(fs, n, xm1);
					break;
				case 2: // n/3 x-1, n/3 (x-1)^2
					augment(fs, n-2*(n/3), xm1);
					augment(fs, n/3, xm1s);
					break;
				case 3: // n/3 x-1, n/3 (x^2-1) -- distinct roots except /F2.
					augment(fs, n-2*(n/3), xm1);
					augment(fs, n/3, x2m1);
					break;
				case 4: // tight stack
					// k(k+1)/2 <= n < (k+1)(k+2)/2
					size_t k = 1; 
					while ((k*(k+1))/2 <= n) ++k; 
					--k;
					size_t l = n-(k*(k+1))/2;
					_R.assign(p, xm1);
					for (size_t i = 1; i <= k; ++i) {
						fs.push_back(p);
						if (i == l) fs.push_back(p);
						_R.mulin(p, xm1);
					}
					break;
				//case 5: // spread stack
				//	break;
				}// switch
			} else { // fsnum is neg
				invariants(fs, n/2, -fsnum);
				randomPolynomial(p,n-n/2);
				augment(fs, 1, p);
			}
			return fs;
		}
			
	}; // SparseMatrixGenerator
} // linbox

#endif // __LINBOX_SPARSE_MATRIX_GENERATOR_H
