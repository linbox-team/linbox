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
		
	private:
		Field _F;
		PolynomialRing _R;
		MatrixDom _MD;
		
	public:
		SparseMatrixGenerator(const Field &F, const PolynomialRing &R): _F(F), _R(R), _MD(F) {}
		
		/**
		 * Reads a file of format:
		 * <number of bumps>
		 * <index of bump 1> <bump 1>
		 * <index of bump 2> <bump 2>
		 * ...
		 * <index of end of frobenius form> 0
		 */
		void readFile(std::vector<Polynomial> &fs, const std::string filename) const {			
			std::ifstream in(filename);
			size_t nbumps;
			in >> nbumps;
			
			Polynomial f;
			_R.assign(f, _R.one);
			_R.write(std::cout << "initial f: ", f) << std::endl;
			
			fs.resize(0);
			
			for (size_t i = 0; i < nbumps; i++) {
				size_t idx;
				in >> idx;
				
				std::cout << "index: " << idx << std::endl;
				std::cout << "size: " << fs.size() << std::endl;
				_R.write(std::cout << "f: ", f) << std::endl;
				
				while (fs.size() < idx) {
					Polynomial tmp;
					_R.assign(tmp, f);
					fs.push_back(tmp);
				}
				
				Polynomial bump;
				_R.read(in, bump);
				
				if (_R.isZero(bump)) {
					return;
				}
				
				_R.mulin(f, bump);
			}
		}
		
		/**
		 * Populates matrix with frobenius form with blocks corresponding to
		 * F_k = C_{f_1 * f_2 * ... * f_k}
		 * return: True if matrix successfully generated.
		 */
		template<class Matrix>
		bool build(Matrix &M, std::vector<Polynomial> fs) {
			size_t min_dim = 0;
			for (size_t i = 0; i < fs.size(); i++) {
				for (size_t j = i; j < fs.size(); j++) {
					min_dim += _R.deg(fs[j]);
				}
			}
			
			if (M.rowdim() < min_dim || M.coldim() < min_dim) {
				std::cout << "Matrix too small (min dim: " << min_dim << ")" << std::endl;
				return false;
			}
			
			size_t offset = 0;
			Polynomial f;
			_R.assign(f, _R.one);
			for (size_t i = 0; i < fs.size(); i++) {
				_R.mulin(f, fs[i]);
				
				
			}
			
			return true;
		}
	};
}

#endif // __LINBOX_SPARSE_MATRIX_GENERATOR_H