#include <linbox/linbox-config.h>
#include <vector>

// Polynomial Matrix / Order Basis
#include "linbox/ring/modular.h"
#include "linbox/ring/givaro-poly.h"
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"

//#define LINBOX_USES_OMP 1
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

#ifndef __TEST_POLY_SMITH_FORM_H
#define __TEST_POLY_SMITH_FORM_H

namespace LinBox
{
	template<class Field>
	class TestPolySmithFormUtil
	{
	public:
		typedef typename Field::Element Element;
		
		typedef MatrixDomain<Field> MatDom;
		typedef typename MatDom::Matrix SubMatrix;
		
		typedef BlasMatrix<Field> Matrix;
		
	private:
		Field _F;
		MatDom _MD;
		
	public:
		TestPolySmithFormUtil(const Field &F): _F(F), _MD(F) {}
		
		template<class Matrix1>
		void printMatrix(const Matrix1 &A) const {
			std::cout << "[" << std::endl;
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
			std::cout << "]" << std::endl;
		}
		
		std::vector<Element> &makeDiag(std::vector<Element> &diag, 
			std::vector<size_t> idx, std::vector<Element> bumps, size_t n) const {
			
			Element elm;
			_F.assign(elm, _F.one);
			
			size_t j = 0;
			for (size_t i = 0; i < n; i++) {
				if (i == idx[j]) {
					_F.mulin(elm, bumps[j]);
					j++;
				}
				
				diag.push_back(elm);
			}
			
			return diag;
		}
		
		void makeExample(Matrix &A, const std::vector<Element> &diag, 
			const std::vector<Element> &lumps) 
		{
			A.zero();
			for (size_t i = 0; i < diag.size(); i++) {
				A.setEntry(i, i, diag[i]);
			}
						
			Matrix L(_F, A.coldim(), A.coldim());
			for (size_t i = 0; i < A.coldim(); i++) {
				L.setEntry(i, i, _F.one);
				
				for (size_t j = 0; j < i; j++) {
					L.setEntry(i, j, lumps[rand() % lumps.size()]);
				}
			}
						
			Matrix U(_F, A.rowdim(), A.rowdim());
			for (size_t i = 0; i < A.rowdim(); i++) {
				U.setEntry(i, i, _F.one);
				
				for (size_t j = i+1; j < A.rowdim(); j++) {
					U.setEntry(i, j, lumps[rand() % lumps.size()]);
				}
			}
						
			_MD.mulin(A, L);
			_MD.rightMulin(A, U);
						
			size_t row_perms = rand() % (2 * A.rowdim());
			for (size_t i = 0; i < row_perms; i++) {
				size_t idx1 = rand() % A.rowdim();
				size_t idx2 = rand() % A.rowdim();
				
				SubMatrix row1(A, idx1, 0, 1, A.coldim());
				SubMatrix row2(A, idx2, 0, 1, A.coldim());
				
				row1.swap(row2);
			}
						
			size_t col_perms = rand() % (2 * A.coldim());
			for (size_t i = 0; i < col_perms; i++) {
				size_t idx1 = rand() % A.coldim();
				size_t idx2 = rand() % A.coldim();
				
				SubMatrix col1(A, 0, idx1, A.rowdim(), 1);
				SubMatrix col2(A, 0, idx2, A.rowdim(), 1);
				
				col1.swap(col2);
			}
			
			size_t row_adds = 1 + (rand() % (2 * A.rowdim()));
			for (size_t i = 0; i < row_adds; i++) {
				for (size_t idx1 = 0; idx1 < A.rowdim(); idx1++) {
					size_t idx2 = rand() % A.rowdim();
					while (idx2 == idx1) {
						idx2 = rand() % A.rowdim();
					}
										
					SubMatrix row1(A, idx1, 0, 1, A.coldim());
					SubMatrix row2(A, idx2, 0, 1, A.coldim());
					_MD.addin(row1, row2);
				}
			}
			
			size_t col_adds = 1 + (rand() % (2 * A.coldim()));
			for (size_t i = 0; i < col_adds; i++) {
				for (size_t idx1 = 0; idx1 < A.coldim(); idx1++) {
					size_t idx2 = rand() % A.coldim();
					while (idx2 == idx1) {
						idx2 = rand() % A.coldim();
					}
					
					SubMatrix col1(A, 0, idx1, A.rowdim(), 1);
					SubMatrix col2(A, 0, idx2, A.rowdim(), 1);
					_MD.addin(col1, col2);
				}
			}
		}
	};
}

#endif // __TEST_POLY_SMITH_FORM_H