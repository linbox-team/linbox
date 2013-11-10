/* linbox/blackbox/rational-matrix-factory.h
 * Copyright (C) 2009 Anna Marszalek
 *
 * Written by Anna Marszalek <aniau@astronet.pl>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */
#ifndef __LINBOX_rational_dense_factory_H
#define __LINBOX_rational_dense_factory_H

#include "linbox/blackbox/factory.h"
//#include "linbox/field/gmp-rational.h"
#include "linbox/field/PID-integer.h"
#include "linbox/integer.h"

#include <vector>

namespace LinBox
{

	/*
	 * aniau@astronet.pl 06/2009
	 * Given rational matrix _matA, computes parameters needed to found best (usually integer) representation:
	 * See: denominator, rationalNorm, normAtilde, normAprim, getOmega
	 * See others: maxNorm, hadamard
	 * Computes the representation
	 * See: makeAtilde, makeAprim
	 * Works with dense matrices, needs sparse spcialization (due to use of row::iterator)
	 */


	template<class Integers, class Rationals, class QMatrix>
	class RationalMatrixFactory /*: public MatrixFactory<Integers,typename Rationals::Element >*/ {
		//typedef GMPRationalField Rationals;
		typedef typename Rationals::Element Quotient;
		//typedef PID_integer Integers;
		//typedef typename PID_integer::Element Integer;
		//

	private:
		const QMatrix*  _matA;
		const Rationals _ratField;

		mutable size_t  _ratOmega;
		mutable size_t  _omega;
		mutable Integer _denA; //once computed becomes constant
		mutable std::vector<Integer> _denAi; //once computed becomes constant

	public:

		RationalMatrixFactory(const QMatrix* A) :
			_ratField(), _denAi(A->rowdim(),1)
		{_ratOmega = 0; _omega = 0; _denA = 1; _matA = A;}

		size_t rowdim() { return _matA->rowdim(); }
		size_t coldim() { return _matA->coldim(); }

		typedef typename QMatrix::ConstIterator ConstIterator;
		typedef typename QMatrix::ConstRowIterator RowIterator;

		double maxNorm(double& res) const {
			typename QMatrix::ConstIterator i;
			res = 0.0;

			for( i = _matA->Begin(); i != _matA->End(); ++i ) {
				double tmp;
				Integer d,n;
				_ratField.get_den(d,*i);
				_ratField.get_num(n,*i);
				tmp = abs( (double)(n)/(double)(d) );
				if( res < tmp ) res = tmp;
			}
			return res;
		}

		double hadamardBound(double& res) const {
			typename QMatrix::ConstRowIterator r;
			typename QMatrix::ConstRow::const_iterator c;

			res = 1.0;

			for( r = _matA->rowBegin(); r != _matA->rowEnd(); ++r ) {
			double temp;
				temp = 0.0;
				for( c = r->begin(); c != r->end(); ++c ) {
					Integer d,n;
					_ratField.get_den(d,*c);
					_ratField.get_num(n,*c);
					double a = (double) n / (double) d;
					temp = temp + (a* a);
				}
				res *= temp;
			}
			res = sqrt(res);
			return res;
		}

		//rational specialization

		/*
		 * Returns the common denominator _denA of _matA
		 */
		Integer& denominator(Integer& da) const {
			if (_denA ==1 ) {
				ConstIterator i;
				//_denA = 1L;
				for (i = _matA->Begin (); i != _matA->End (); ++i) {
					Integer d;
					_ratField.get_den(d,*i);
					lcm(_denA,_denA,d);
				}
			}
			return da=_denA;;
		}

		/*
		 * returns common denominator _denAi[(size_t)i] of i-th row
		 */
		Integer& denominator(Integer& di, const int i) const {
			if (_denAi[(size_t)i]==1) {
				for (size_t j=0; j < _matA->coldim(); ++j) {
					Integer d;
					_ratField.get_den(d,_matA->getEntry((size_t)i,(size_t)j));
					lcm(_denAi[(size_t)i],_denAi[(size_t)i],d);
				}
			}
			return di=_denAi[(size_t)i];
		}

		//returns max of abs(numerators) and denominators of _matA
		Integer& rationalNorm(Integer& res) const {
			ConstIterator i;
			res = 0L;
			Integer tmp;
			for (i = _matA->Begin (); i != _matA->End (); ++i) {
				Integer n ;
				_ratField.get_num(n,*i);
				Integer d ;
				_ratField.get_den(d,*i);
				tmp = abs (n);
				if (tmp < d) tmp = d;
				if (res < tmp) res = tmp;
			}
			return res;
		}

		//returns norm of A'= _denA * _matA
		Integer& normAprim(Integer& res) const {

			Integer DA;
			denominator(DA);

			double norm; maxNorm(norm);
			res = (Integer) ( (double) DA * norm ) ;

			return res;
		}

		//returns norm of tilde{A} = diag(_denAi[(size_t)i])_matA
		Integer& normAtilde(Integer& res) const {
			res = 0L;
			double dres = 0;
			//int i=0;
			Integer tmp;
			double dtmp;
			for (int i=0; i < _matA->rowdim(); ++i) {
				Integer di; denominator(di,i);
				for (int j=0; j < _matA->coldim(); ++j ) {
					Integer n ;
					_ratField.get_num(n,_matA->getEntry(i,j));
					Integer d ;
					_ratField.get_den(d,_matA->getEntry(i,j));
					dtmp = (double)di/double(d)*double(n);
					//tmp = di/d;
					//tmp*=abs(n);
					if (dtmp > dres) dres = dtmp;
				}
			}
			res = (Integer)dres;

			return res;
		}

		/*
		 * optimization: computes normAprim, normAprim, rationalNormat the same time
		 */
		Integer getNorms(Integer& ratnorm, Integer& normaprim, Integer& normatilde) const {
			ratnorm = 0L; normaprim=0L; normatilde= 0L;
			Integer da=1L;
			std::vector<integer> di(_matA->rowdim(),1L);

			for (size_t i=0; i < _matA->rowdim(); ++i)  {
				if (_denAi[(size_t)i]==1) {
					for (size_t j=0; j < _matA->coldim(); ++j ) {
						Integer d ;
						_ratField.get_den(d,_matA->getEntry(i,j));
						lcm(_denAi[(size_t)i],_denAi[(size_t)i],d);
					}
				}
				di[(size_t)i] = _denAi[(size_t)i];
				lcm(_denA,_denA,di[(size_t)i]);
			}
			da = _denA;

			for (size_t i=0; i < _matA->rowdim(); ++i )  {
				for (size_t j=0; j < _matA->coldim(); ++j ) {
					Integer n ;
					_ratField.get_num(n,_matA->getEntry(i,j));
					Integer d ;
					_ratField.get_den(d,_matA->getEntry(i,j));

					Integer tmp = abs(n);
					if (tmp > ratnorm) ratnorm = tmp;
					if (d > ratnorm) ratnorm = d;

					Integer tmp2 = (di[(size_t)i]) / d;
					tmp2 *=tmp;
					if (tmp2 > normatilde) normatilde = tmp2;

					tmp2 = da/d;
					tmp2 *= tmp;
					if (tmp2 > normaprim) normaprim = tmp2;

				}
			}
			integer minnorm = (ratnorm > normatilde) ? normatilde : ratnorm;
			return minnorm;
		}

		/*
		 * Creates Aprim = _denA * _matA
		 */
		template <class Matrix>
		Matrix& makeAprim(Matrix& Aprim) const {
			Integer da; denominator(da);
			Aprim.resize(_matA->rowdim(),_matA->coldim());

			for( size_t i=0; i < _matA->rowdim(); ++i) {
				for (size_t j=0; j < _matA->coldim(); ++j) {
					Quotient  Aij =  _matA->getEntry(i,j);
					Integer n ;
					_ratField.get_num(n,Aij);
					Integer d ;
					_ratField.get_den(d,Aij);
					Integer tmp = da/d;
					tmp *=n;
					typename Matrix::Field F=Aprim.field();
					typename Matrix::Field::Element ftmp; F.init(ftmp,tmp);
					Aprim.setEntry(i,j, ftmp);
				}
			}
			return Aprim;
		}

		/*
		 * Creates Atilde = diag(_denAi[(size_t)i]) * _matA
		 */
		template <class Matrix>
		Matrix& makeAtilde(Matrix& Atilde) const {

			Atilde.resize(_matA->rowdim(),_matA->coldim());
			std::vector<integer> di(_matA->rowdim());
			for (size_t i=0; i < (size_t)_matA->rowdim(); ++i)
				denominator(di[(size_t)i],(int)i);

			for( size_t i=0; i < _matA->rowdim(); ++i) {
				for (size_t j=0; j < _matA->coldim(); ++j) {
					Quotient  Aij = _matA->getEntry(i,j);
					Integer n ;
					_ratField.get_num(n,Aij);
					Integer d ;
					_ratField.get_den(d,Aij);
					Integer tmp = di[(size_t)i]/d;
					tmp *=n;
					typename Matrix::Field F=Atilde.field();
					typename Matrix::Field::Element ftmp; F.init(ftmp,tmp);
					Atilde.setEntry(i,j, ftmp);
				}
			}
			return Atilde;
		}

		/*
		 * Counts the number of non-zero and non-integer elements
		 */
		size_t getOmega(size_t& o, size_t& ro) const {
			if (_omega ==0 ) {
				for( size_t i=0; i < _matA->rowdim(); ++i) {
					for (size_t j=0; j < _matA->coldim(); ++j) {
						Quotient  Aij =  _matA->getEntry(i,j);
						Integer n ;
						_ratField.get_num(n,Aij);
						Integer d ;
						_ratField.get_den(d,Aij);
						if (n!=0) {
							++_omega;
							if (d!=1) ++_ratOmega;
						}
					}
				}
			}
			ro = _ratOmega;
			return o=_omega;
		}
	};

} // namespace LinBox


#endif  //__LINBOX_rational_dense_factory_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

