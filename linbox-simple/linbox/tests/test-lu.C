#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <linbox/field/modular-int.h>
#include <linbox/solutions/lu.h>
#include <linbox/blackbox/dense.h>
#include <linbox/util/commentator.h>
#include <vector>
#include <iomanip>

typedef LinBox::Modular<int> Field;

namespace LinBox
{
	/*- test if M==L*U;
	 * where L is lower triangular matrix, whose diagonal entries are 1.
	 * U is an upper triangular matrix.
	 */
	template<class Field>
	bool LU_MUL_TEST(const DenseMatrix<Field>& M, const DenseMatrix<Field>& L,const DenseMatrix<Field>& U) {

		linbox_check( (M.rowdim() == L. rowdim()) && (M. coldim() == U. coldim())&&
			      (L.rowdim() == L. coldim()) && (U. rowdim() == U. coldim())&&
			      (L.coldim() == U. rowdim()));
      
		typename DenseMatrix<Field>::ConstColIterator ccolp;

		typename DenseMatrix<Field>::Element e;

		typename DenseMatrix<Field>::ConstRowIterator crowp1;   

		typename DenseMatrix<Field>::ConstColIterator ccolp2;       

		typename DenseMatrix<Field>::ConstCol::const_iterator cep;

		typename DenseMatrix<Field>::ConstRow::const_iterator cep1,cep1e;

		typename DenseMatrix<Field>::ConstCol::const_iterator cep2;

		int i, j;

		i=0;

		for (ccolp = M. colBegin(), ccolp2 = U. colBegin();
		     ccolp != M. colEnd();
		     ++ ccolp, ++ ccolp2, ++ i) {

			j=0;

			for (cep = ccolp -> begin(), crowp1 = L. rowBegin();
			     cep != ccolp -> end();
			     ++ cep, ++ crowp1, ++j) {

				if ( j <= i)

					cep1e = crowp1 -> begin() + j;

				else

					cep1e = crowp1 -> begin() + (i + 1);
              
				M. field(). init (e, 0);

				for (cep1 = crowp1 -> begin(), cep2 = ccolp2 -> begin();
				     cep1 != cep1e;
				     ++ cep1, ++ cep2)

					M. field(). axpyin (e, *cep1, *cep2);
              
				if (j <= i)

					M. field(). addin (e, *cep2);

				if (! M. field(). areEqual (e, *cep))

					return false;
			}     
		}
		return true;
	}
}

bool test (int _SIZE) {

	Field field (1073741789UL);

	Field::Element e;

	typedef LinBox::DenseMatrix<Field>  Matrix;

	Matrix M (field, _SIZE, _SIZE);


	for (int i = 0; i < (int)M. rowdim(); ++i)

		for (int j = 0; j < (int)M. coldim(); ++j) {

			field. init (e,  random() % 100 + 1);

			M. setEntry (i, j, e); 
		}
 
	LinBox::DenseMatrix<Field> M_C(M);

	LinBox::LU (M);

	return LinBox::LU_MUL_TEST(M_C,M,M);
}
  
int main() {

	srand(time(0));

	int iteration = 120;

	LinBox::commentator.start("Test LU over random matrix in modular int field","",iteration-10);

	for(int i = 10;i < iteration; ++i) {

		LinBox::commentator.progress();

		if (!test (i)) {

			LinBox::commentator.stop("Failed","");

			return -1;
		}
	}

	LinBox::commentator.stop("passed","");

	return 0;
}
