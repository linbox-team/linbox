/* linbox/blackbox/rational-matrix-factory.h
 * Copyright (C) 2009 Anna Marszalek
 * 
 * Written by Anna Marszalek <aniau@astronet.pl> 
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#ifndef __LINBOX_rational_dense_factory_H
#define __LINBOX_rational_dense_factory_H

#include <linbox/blackbox/factory.h>
#include <linbox/blackbox/dense.h>
//#include <linbox/matrix/dense.h>
//#include <linbox/field/gmp-rational.h>
#include <linbox/field/PID-integer.h>
#include <linbox/integer.h>

#include <vector>

namespace LinBox 
{

/*
 * aniau@astronet.pl 06/2009
 * Given rational matrix _A, computes parameters needed to found best (usually integer) representation:
 * See: denominator, rationalNorm, normAtilde, normAprim, getOmega
 * See others: maxNorm, hadamard 
 * Computes the representation
 * See: makeAtilde, makeAprim
 * Works with dense matrices, needs sparse spcialization (due to use of row::iterator)  
 */


template<class Integers, class Rationals, class QMatrix>
class RationalMatrixFactory //: public DenseMatrixFactory<Integers,typename Rationals::Element >
{
//typedef GMPRationalField Rationals;
typedef typename Rationals::Element Quotient;
//typedef PID_integer Integers;
//typedef typename PID_integer::Element Integer;
//

    private:
    	const QMatrix* _A;
	const Rationals Q;

	mutable size_t rat_omega;
	mutable size_t omega;
	mutable Integer denA; //once computed becomes constant
	mutable std::vector<Integer> denAi; //once computed becomes constant

    public:

    	RationalMatrixFactory(const QMatrix* A):  Q(), denAi(A->rowdim(),1) {rat_omega = 0; omega = 0; denA = 1; _A = A;}

	size_t rowdim() { return _A->rowdim(); }
        size_t coldim() { return _A->coldim(); }

	typedef typename QMatrix::ConstRawIterator ConstRawIterator;
        typedef typename QMatrix::ConstRowIterator RowIterator;

	double maxNorm(double& res) const {
		typename QMatrix::ConstRawIterator i;
	        res = 0.0;
	        double tmp;

		for( i = _A->rawBegin(); i != _A->rawEnd(); ++i ) {
			Integer d,n; Q.get_den(d,*i);Q.get_num(n,*i);
			tmp = abs( (double)(n)/(double)(d) );
			if( res < tmp ) res = tmp;
		}
	        return res;
	}

	double hadamardBound(double& res) const {
	        typename QMatrix::ConstRowIterator r;
	        typename QMatrix::ConstRow::const_iterator c;

	        res = 1.0;
	        double temp;

		for( r = _A->rowBegin(); r != _A->rowEnd(); ++r ) {
			temp = 0.0;
			for( c = r->begin(); c != r->end(); ++c ) {
				Integer d,n; Q.get_den(d,*c);Q.get_num(n,*c);
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
	 * Returns the common denominator denA of _A
	 */
        Integer& denominator(Integer& da) const {
	    if (denA ==1 ) {
	        ConstRawIterator i;
	        //denA = 1L;
	        for (i = _A->rawBegin (); i != _A->rawEnd (); ++i) {
	 	        Integer d; Q.get_den(d,*i);
		        lcm(denA,denA,d);
	        }
	    }	
	    return da=denA;;
	}

	/*
	 * returns common denominator denAi[i] of i-th row
	 */
	Integer& denominator(Integer& di, const int i) const {
	  if (denAi[i]==1) {   
		typedef typename QMatrix::ConstRow::const_iterator EltIterator;
		for (size_t j=0; j < _A->coldim(); ++j) {
			Integer d; Q.get_den(d,_A->getEntry(i,j));
			lcm(denAi[i],denAi[i],d);
		}
	  }
	  return di=denAi[i];
	}

//returns max of abs(numerators) and denominators of _A 
	Integer& rationalNorm(Integer& res) const {
	        ConstRawIterator i;
	        res = 0L;
	        Integer tmp;
	        for (i = _A->rawBegin (); i != _A->rawEnd (); ++i) {
		        Integer n ;Q.get_num(n,*i);
		        Integer d ;Q.get_den(d,*i);
			tmp = abs (n);
			if (tmp < d) tmp = d;
			if (res < tmp) res = tmp;
		}
		return res;
	}

//returns norm of A'= denA * _A
	Integer& normAprim(Integer& res) const {
	
                Integer DA;
                denominator(DA);
		
		double norm; maxNorm(norm);
		res = (Integer) ( (double) DA * norm ) ;

                return res;
        }

//returns norm of tilde{A} = diag(denAi[i])_A
        Integer& normAtilde(Integer& res) const {
	        //typedef typename DenseMatrixBase<Quotient>::ConstRow::const_iterator EltIterator;
                res = 0L;
		double dres = 0;
                //int i=0;
                Integer tmp;
		double dtmp;
                for (int i=0; i < _A->rowdim(); ++i) {
	                Integer di; denominator(di,i);
 	                for (int j=0; j < _A->coldim(); ++j ) {
	                        Integer n ; Q.get_num(n,_A->getEntry(i,j));
                                Integer d ; Q.get_den(d,_A->getEntry(i,j));
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
		typedef typename QMatrix::ConstRow::const_iterator EltIterator;
		ratnorm = 0L; normaprim=0L; normatilde= 0L;
		Integer da=1L;
		std::vector<integer> di(_A->rowdim(),1L);

		for (size_t i=0; i < _A->rowdim(); ++i)  {
		   if (denAi[i]==1) { 	
			for (size_t j=0; j < _A->coldim(); ++j ) {
		        	Integer d ; Q.get_den(d,_A->getEntry(i,j));
				lcm(denAi[i],denAi[i],d);
			}
		   } 
		   di[i] = denAi[i];
		   lcm(denA,denA,di[i]);
		}
		da = denA;

		for (size_t i=0; i < _A->rowdim(); ++i )  {
			for (size_t j=0; j < _A->coldim(); ++j ) {
				Integer n ; Q.get_num(n,_A->getEntry(i,j));
				Integer d ; Q.get_den(d,_A->getEntry(i,j));
				
				Integer tmp = abs(n);
				if (tmp > ratnorm) ratnorm = tmp;
				if (d > ratnorm) ratnorm = d;

				Integer tmp2 = (di[i]) / d;
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
	DenseMatrix< Rationals >& makeA(DenseMatrix<Rationals >& A) {
	        //DenseMatrix<Rationals > local(Q,_A);
		//return A = local;
		
		for( int i=0; i < _A->rowdim(); ++i) {
	                for (int j=0; j < _A->coldim(); ++j) {
	                        A.setEntry(i,j, _A->getEntry(i,j));
			}
		}
	}
*/
	/*
	 * Creates Aprim = denA * _A
	 */
	template <class Matrix>
	Matrix& makeAprim(Matrix& Aprim) const {
		Integer da; denominator(da);
		Aprim.resize(_A->rowdim(),_A->coldim());

                for( size_t i=0; i < _A->rowdim(); ++i) {
			for (size_t j=0; j < _A->coldim(); ++j) {
		                Quotient  q =  _A->getEntry(i,j);
			        Integer n ; Q.get_num(n,q);
			        Integer d ; Q.get_den(d,q);
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
	 * Creates Atilde = diag(denAi[i]) * _A
	 */
	template <class Matrix>
        Matrix& makeAtilde(Matrix& Atilde) const {

		Atilde.resize(_A->rowdim(),_A->coldim());
	        std::vector<integer> di(_A->rowdim());
	        for (size_t i=0; i < _A->rowdim(); ++i) denominator(di[i],i);
		           
		for( size_t i=0; i < _A->rowdim(); ++i) {
                        for (size_t j=0; j < _A->coldim(); ++j) {
                        	Quotient  q = _A->getEntry(i,j);
	                        Integer n ; Q.get_num(n,q);
        	                Integer d ; Q.get_den(d,q);
			        Integer tmp = di[i]/d;
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
		if (omega ==0 ) {
		    for( size_t i=0; i < _A->rowdim(); ++i) {
                        for (size_t j=0; j < _A->coldim(); ++j) {
				Quotient  q =  _A->getEntry(i,j);
				Integer n ; Q.get_num(n,q);
 		                Integer d ; Q.get_den(d,q);
				if (n!=0) {
					++omega;
					if (d!=1) ++rat_omega;
				}
			}
		    }
		}
		ro = rat_omega;
		return o=omega;
	}
};

} // namespace LinBox

#include "dense.inl"

#endif  //__LINBOX_rational_dense_factory_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
