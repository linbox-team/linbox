#include <linbox/linbox-config.h>

#include <iostream>
#include <vector>

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
		typedef MatrixDomain<Field> MatrixDom;
		typedef typename MatrixDom::Matrix SubMatrix;
		
	private:
		Field _F;
		PolynomialRing _R;
		
	public:
		SparseMatrixGenerator(const Field &F, const PolynomialRing &R): _F(F), _R(R) {}
		
		template<class Matrix1>
		void printMatrix(const Matrix1 &A) const {
			std::cout << "matrix(K, " << A.rowdim() << ", " << A.coldim() << ", [" << std::endl;
			for (size_t i = 0; i < A.rowdim(); i++) {
				std::cout << "\t[";
				for (size_t j = 0; j < A.coldim(); j++) {
					_F.write(std::cout, A.getEntry(i, j));
					if (j < A.coldim() - 1) {
						std::cout << ", ";
					}
				}
				std::cout << "]";
				if (i < A.rowdim() - 1) {
					std::cout << ",";
				}
				std::cout << std::endl;
			}
			std::cout << "])" << std::endl;
		}
		
		/**
		 * Reads a file of format:
		 * <number of bumps>
		 * <multiplicity 1> <bump 1>
		 * <multiplicity 2> <bump 2>
		 * ...
		 * <multiplicity n> <bump n>
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
				_R.read(in, bump);
				_R.write(std::cout << "bump: ", bump) << std::endl;
				
				for (size_t j = 0; j < multiplicity; j++) {
					Polynomial tmp;
					_R.assign(tmp, bump);
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
			_R.write(std::cout << "monic f: ", monic_f) << std::endl;
			
			std::vector<integer> coeffs;
			_R.convert(coeffs, monic_f);
			
			for (size_t i = 0; i < coeffs.size() - 1; i++) {
				Element tmp;
				_F.init(tmp, coeffs[i]);
				_F.negin(tmp);
				
				M.setEntry(offset + i, offset + dim - 1, tmp);
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
				std::cout << "Matrix too small (min dim: " << min_dim << ")" << std::endl;
				return false;
			}
			
			for (size_t i = 0; i < M.rowdim(); i++) {
				for (size_t j = 0; j < M.coldim(); j++) {
					M.setEntry(i, j, _F.zero);
				}
			}
			
			size_t offset = 0;
			Polynomial f;
			for (size_t i = 0; i < fs.size(); i++) {				
				_R.assign(f, fs[i]);
				size_t d = _R.deg(f);
				
				makeCompanion(M, f, offset);
				
				offset += d;
			}
			
			return true;
		}
		
		template<class Matrix>
		double sparsity(const Matrix &M) {
			double nnz = 0;
			
			for (size_t i = 0; i < M.rowdim(); i++) {
				for (size_t j = 0; j < M.coldim(); j++) {
					if (_F.isZero(M.getEntry(i,j))) {
						continue;
					}
					
					nnz++;
				}
			}
			
			return nnz / (M.rowdim() * M.coldim());
		}
		
		// Next in in [0, limit), not equal to except
		size_t nextInt(size_t limit, size_t except) {
			size_t rv = rand() % limit;
			
			while(rv == except) {
				rv = rand() % limit;
			}
			
			return rv;
		}
		
		template<class Matrix>
		void addRow(Matrix &M, size_t row1, size_t row2, Element &z) {
			for (size_t col = 0; col < M.coldim(); col++) {
				Element tmp, a, b;
				
				M.getEntry(a, row1, col);
				M.getEntry(b, row2, col);
				
				_F.mul(tmp, b, z);
				_F.addin(tmp, a);
				
				M.setEntry(row1, col, tmp);
				
				Element tmp2;
				M.getEntry(tmp2, row1, col);
				
				_F.write(std::cout, tmp) << " = ";
				_F.write(std::cout, a) << " + ";
				_F.write(std::cout, b) << " * ";
				_F.write(std::cout, z) << " = ";
				_F.write(std::cout, tmp2) << std::endl;
				
				assert(_F.areEqual(tmp, tmp2));
			}
		}
		
		template<class Matrix>
		void addCol(Matrix &M, size_t col1, size_t col2, Element &z) {
			for (size_t row = 0; row < M.rowdim(); row++) {
				Element tmp, a, b;
				
				M.getEntry(a, row, col1);
				M.getEntry(b, row, col2);
				
				_F.mul(tmp, b, z);
				_F.addin(tmp, a);
				
				M.setEntry(row, col1, tmp);
				
				Element tmp2;
				M.getEntry(tmp2, row, col1);
				
				_F.write(std::cout, tmp) << " = ";
				_F.write(std::cout, a) << " + ";
				_F.write(std::cout, b) << " * ";
				_F.write(std::cout, z) << " = ";
				_F.write(std::cout, tmp2) << std::endl;
				
				assert(_F.areEqual(tmp, tmp2));
			}
		}
		
		/**
		 *
		 */
		template<class Matrix>
		void fillIn(Matrix &M, double targetSparsity) {
			size_t dim = M.rowdim();
			
			printMatrix(M);
			
			while (sparsity(M) < targetSparsity) {
				Element z;
				do {
					_F.init(z, rand());
				} while(_F.isZero(z));
				_F.write(std::cout << "Scale Row: ", z) << std::endl;
				
				size_t a = nextInt(dim, -1);
				size_t b = nextInt(dim, a);
				
				std::cout << "rows: " << a << ", " << b << std::endl;
				
				addRow(M, a, b, z);
				printMatrix(M);
				
				_F.negin(z);
				_F.write(std::cout << "Scale Col: ", z) << std::endl;
				std::cout << "cols: " << b << ", " << a << std::endl;
				addCol(M, b, a, z);
				printMatrix(M);
				std::cout << std::endl;
				// return;
			}
			
			printMatrix(M);
		}
	};
}

#endif // __LINBOX_SPARSE_MATRIX_GENERATOR_H