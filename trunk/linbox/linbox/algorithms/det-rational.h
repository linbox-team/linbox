/* linbox/blackbox/rational-reconstruction-base.h
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

#ifndef __LINBOX_det_rational_H
#define __LINBOX_det_rational_H

#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"
#include "linbox/field/modular-double.h"

//#include "linbox/field/gmp-rational.h"
#include "linbox/field/PID-integer.h"
#include "linbox/blackbox/rational-matrix-factory.h"
#include "linbox/blackbox/dense.h"
#include "linbox/algorithms/varprec-cra-early-single.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/rational-reconstruction-base.h"
#include "linbox/algorithms/classic-rational-reconstruction.h"
#include "linbox/solutions/det.h"
//#include "linbox/blackbox/apply.h"
#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/algorithms/rational-solver.h"

namespace LinBox
{
/*
 * Computes the determinant of a rational dense matrix
 */

/* looks for the same integer value by 2 methods */
template<class T1, class T2>
struct MyModularDet{
	T1* t1;
	T2* t2;

	int switcher;

	MyModularDet(T1* s1, T2* s2, int s = 1)  {t1=s1; t2=s2;switcher = s;}

	int setSwitcher(int s) {return switcher = s;}

	template<typename Int, typename Field>
	Int& operator()(Int& P, const Field& F) const {
		if (switcher ==1) {
			t1->operator()(P,F);
		} else {
			t2->operator()(P,F);
		}
		return P;
	}
};

/* PrecDet variant of algorithm 
 * mul is D, and div is d>1 if corrections are run
 * */
template <class Blackbox, class MyMethod>
struct MyRationalModularDet {
	const Blackbox &A;
	const MyMethod &M;
	const Integer &mul;//multiplicative prec;
	const Integer &div;

	MyRationalModularDet(const Blackbox& b, const MyMethod& n, 
				  const Integer & p1, const Integer & p2): A(b), M(n), mul(p1), div(p2) {}
	MyRationalModularDet(MyRationalModularDet& C) :
		MyRationalModularDet(C.A,C.M,C.mul)
	{}
	
	void setDiv (const Integer& d) {div = d;}
	void detMul (const Integer& m) {mul = m;}

	template<typename Int, typename Field>
	Int& operator()(Int& P, const Field& F) const {
		typedef typename Blackbox::template rebind<Field>::other FBlackbox;
		FBlackbox Ap(A, F);
		det (P, Ap, typename FieldTraits<Field>::categoryTag(), M);
		typename Field::Element e;
		F.init(e, mul);
		F.mulin(P,e);
		F.init(e, div);
		return F.divin(P,e);
	}
};

/* PrecMat variant of algorithm 
 * no multiple. no divisor
 */
template <class Blackbox, class MyMethod>
struct MyIntegerModularDet {
	const Blackbox &A;
	const MyMethod &M;

	MyIntegerModularDet(const Blackbox& b, const MyMethod& n)
	: A(b), M(n) {}

	MyIntegerModularDet(MyIntegerModularDet& C) : MyIntegerModularDet(C.A,C.M){}

	template<typename Int, typename Field>
	Int& operator()(Int& P, const Field& F) const {
		typedef typename Blackbox::template rebind<Field>::other FBlackbox;
		FBlackbox Ap(A, F);
		return det( P, Ap, typename FieldTraits<Field>::categoryTag(), M);
	}
};

/*
 * Corrects the matrix by removing gcd of columns
 */

template <class Blackbox>
typename Blackbox::Field::Element& corrections(Blackbox& IntBB, typename Blackbox::Field::Element& f) {
	f = 1;
	for (size_t j=0; j < IntBB.coldim(); ++j) {
		Integer g = 0;
		for (size_t i=0; i < IntBB.rowdim(); ++i) {
			gcd(g,g,IntBB.getEntry(i,j));
		}
		f*=g;
		for (size_t i=0; i < IntBB.rowdim(); ++i) {
			Integer x = IntBB.getEntry(i,j);
		        IntBB.setEntry(i,j,x/g);
		}
	}
	return f;
}

template <class Rationals, class MyMethod >
typename Rationals::Element& rational_det (typename Rationals::Element &d,
	                                                      const DenseMatrix<Rationals > &A,
							      const MyMethod &Met=  Method::Hybrid())
{

typedef Modular<double> myModular;
typedef typename Rationals::Element Quotient;

commentator.start ("Rational Det", "Rdeterminant");

RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));

Integer F = 1; 
Integer M = 1;

//DenseMatrixBase<Quotient> ABase(A);
RationalMatrixFactory<PID_integer,Rationals, DenseMatrix<Rationals > > FA(&A);
Integer di=1;

for (int i=A.rowdim()-1; i >= 0 ; --i) {
	FA.denominator(di,i);
	M *=di;
}

PID_integer Z;
DenseMatrix<PID_integer> Atilde(Z,A.rowdim(), A.coldim());
FA.makeAtilde(Atilde);

UserTimer t0, t1,t2;bool term = false;
t0.clear();
t0.start();

corrections(Atilde,F);

ChineseRemainder< VarPrecEarlySingleCRA<Modular<double> > > cra(3UL);
MyRationalModularDet<DenseMatrix<Rationals > , MyMethod> iteration1(A, Met, M, F);
MyIntegerModularDet<DenseMatrix<PID_integer>, MyMethod> iteration2(Atilde, Met);
MyModularDet<MyRationalModularDet<DenseMatrix<Rationals > , MyMethod>, 
		  MyIntegerModularDet<DenseMatrix<PID_integer>, MyMethod> >  iteration(&iteration1,&iteration2);

RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > RR;

Integer dd; // use of integer due to non genericity of cra. PG 2005-08-04
size_t k = 4;
t1.clear();
t2.clear();
t1.start();
	term = cra(k,dd,iteration1,genprime);
t1.stop();
t2.start();
	term = cra(k,dd,iteration2,genprime);
t2.stop();

double s1 = t1.time(), s2 = t2.time();

if (t1.time() < t2.time()) {
	//cout << "ratim " << flush;
	iteration.setSwitcher(1);
} else {
	//cout << "intim " << flush;
	iteration.setSwitcher(2);
	s1 = s2;
}

//switch: LIF
vector<typename PID_integer::Element> v(A.rowdim()), r(A.rowdim());
++genprime;
for(size_t i=0; i < v.size(); ++i) {
	v[i] = rand() % (*genprime) ;
}
t2.clear();
t2.start();
for (size_t i=0; i < k; ++i)
	Atilde.apply(r,v);
t2.stop();
s2 = t2.time();

Integer lif = 1;
if ((s1 > 4*s2) && (!term)){
	//cout << "lif " << flush;
	RationalSolver < PID_integer , Modular<double>, RandomPrimeIterator, DixonTraits > RSolver;
	LastInvariantFactor < PID_integer ,RationalSolver < PID_integer, Modular<double>, RandomPrimeIterator, DixonTraits > >  LIF(RSolver);
	vector<Integer> r_num2 (Atilde. coldim());
t1.clear();
t1.start();	
	if (LIF.lastInvariantFactor1(lif, r_num2, Atilde)==0) {
		cerr << "error in LIF\n" << flush;
	        dd = 0;
	}
t1.stop();
	cra.changePreconditioner(lif,1);
} 

k *= 2;
t1.clear();t2.clear();
while (! cra(k,dd, iteration, genprime)) {
	k *=2;
	Integer m,res,res1; //Integer r; Integer a,b;
	cra.getModulus(m);
	cra.getResidue(res);//no need to divide
	Integer D_1; 
	inv(D_1,M,m);
	res = (res*D_1) % m;
	res1 = (res*lif) % m;

	Integer den;
	Integer num;

	bool rrterm = false;

t1.start();
 if ((rrterm = RR.reconstructRational(num,den,res1,m))) {
		t1.stop();
		cra.changePreconditioner(num*M,den*F);
	        if (cra(5,dd,iteration,genprime) ) {
		        cra.getResidue(dd);
			num *=dd;
		 	//cout << "without lif " << flush;
		} else rrterm = false;
	} else 	{//we use num approx by lif
		t1.stop();
		t2.start();
		if ((rrterm = RR.reconstructRational(num,den,res,m))) {
			t2.stop();
			cra.changePreconditioner(num*M,den*F*lif);
			if (cra(5,dd,iteration,genprime) ) {
				cra.getResidue(dd);
				num *=dd;
				den *=lif;
				//cout << "with lif " << flush;
			} else rrterm = false;
		} else t2.stop();
	}	
	if (rrterm) {	
			integer t,tt;
			Rationals Q;
	        	Q.init(d, num, den);
			Q.get_den(t,d);Q.get_num(tt,d);
			//cout << "terminated by RR at step "<< (int)(log(k)/log(2)) << " in " << flush;
			t0.stop();
			return d;	
	} else {
		//cout << "rollback" << flush;
		cra.changePreconditioner(lif,1);
	}
}
cra.result(dd);

integer t; 

Rationals Q;
dd *=F;
Q.init(d, dd,M);
Q.get_den(t, d);
//err = M/t;
//cout << "terminated by ET at step "<< (int)(log(k)/log(2)) << "in " << flush;

commentator.stop ("done", NULL, "Iminpoly");

t0.stop();
return d;

}

}

#endif //__LINBOX_det_rational_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
