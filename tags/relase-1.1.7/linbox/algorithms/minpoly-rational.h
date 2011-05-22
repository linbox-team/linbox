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

#ifndef __LINBOX_rat_minpoly_H
#define __LINBOX_rat_minpoly_H

#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"
#include "linbox/field/modular-double.h"

//#include "linbox/field/gmp-rational.h"
#include "linbox/field/PID-integer.h"
#include "linbox/blackbox/rational-matrix-factory.h"
#include "linbox/blackbox/dense.h"
#include "linbox/algorithms/cra-early-multip.h"
#include "linbox/algorithms/cra-domain.h"
//#include "linbox/algorithms/rational-cra.h"
#include "linbox/algorithms/rational-reconstruction-base.h"
#include "linbox/algorithms/classic-rational-reconstruction.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/diagonal.h"

namespace LinBox
{

/*
 * Computes the minimla polynomial of a rational dense matrix
 */

template<class T1, class T2>
struct MyModularMinpoly{
	T1* t1;
	T2* t2;

	int switcher;

	MyModularMinpoly(T1* s1, T2* s2, int s = 1)  {t1=s1; t2=s2;switcher = s;}

	int setSwitcher(int s) {return switcher = s;}

	template<typename Polynomial, typename Field>
	Polynomial& operator()(Polynomial& P, const Field& F) const {
		if (switcher ==1) {
			t1->operator()(P,F);
		} else {
			t2->operator()(P,F);
		}
		return P;
	}
};

template <class Blackbox, class MyMethod>
struct MyRationalModularMinpoly {
	const Blackbox &A;
	const MyMethod &M;
	const std::vector<Integer> &mul;//multiplicative prec;

	MyRationalModularMinpoly(const Blackbox& b, const MyMethod& n, 
				  const std::vector<Integer >& p): A(b), M(n), mul(p) {}
	MyRationalModularMinpoly(MyRationalModularMinpoly& C) :
		MyRationalModularMinpoly(C.A,C.M,C.mul)
	{}

	template<typename Polynomial, typename Field>
	Polynomial& operator()(Polynomial& P, const Field& F) const {
		typedef typename Blackbox::template rebind<Field>::other FBlackbox;
		FBlackbox * Ap;
		MatrixHom::map(Ap, A, F);
		minpoly( P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
		int s = A.rowdim() +1 - P.size();//shift in preconditioning
		typename std::vector<Integer >::const_iterator it = mul.begin() + s;
		typename Polynomial::iterator it_p = P.begin();
		for (;it_p !=P.end(); ++it, ++it_p) {
			typename Field::Element e;
			F.init(e, *it);
			F.mulin(*it_p,e);
		}

		delete Ap;
		return P;
	}
};

template <class Blackbox, class MyMethod>
struct MyIntegerModularMinpoly {
	const Blackbox &A;
	const MyMethod &M;
	const std::vector<typename Blackbox::Field::Element> &vD;//diagonal div. prec;
	const std::vector<typename Blackbox::Field::Element > &mul;//multiplicative prec;

	MyIntegerModularMinpoly(const Blackbox& b, const MyMethod& n, const std::vector<typename Blackbox::Field::Element>& ve,
					const std::vector<typename Blackbox::Field::Element >& p	)
	: A(b), M(n), vD(ve), mul(p) {}

	MyIntegerModularMinpoly(MyIntegerModularMinpoly& C) : MyIntegerModularMinpoly(C.A,C.M,C.vD,C.mul){}

	template<typename Polynomial, typename Field>
	Polynomial& operator()(Polynomial& P, const Field& F) const {
		typedef typename Blackbox::template rebind<Field>::other FBlackbox;
		FBlackbox * Ap;
		MatrixHom::map(Ap, A, F);

		typename std::vector<typename Blackbox::Field::Element>::const_iterator it;
		
		int i=0;
		for (it = vD.begin(); it != vD.end(); ++it,++i) {
			typename Field::Element t,tt;
			F.init(t,*it);
			F.invin(t);
			for (int j=0; j < A.coldim(); ++j) {
				F.mulin(Ap->refEntry(i,j),t);
			}
		}
			
		minpoly( P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
		int s = A.rowdim() +1 - P.size();//shift in preconditioning
		typename std::vector<typename Blackbox::Field::Element >::const_iterator it2 = mul.begin() + s;
                typename Polynomial::iterator it_p = P.begin();
                for (;it_p !=P.end(); ++it2, ++it_p) {
                        typename Field::Element e;
                        F.init(e, *it2);
                        F.mulin(*it_p,e);
                }

		delete Ap;
		return P;
	}
};

template <class Rationals, template <class> class Vector, class MyMethod >
Vector<typename Rationals::Element>& rational_minpoly (Vector<typename Rationals::Element> &p,
	                                                      const DenseMatrix<Rationals > &A,
							      const MyMethod &Met=  Method::Hybrid())
{

typedef Modular<double> myModular;
typedef typename Rationals::Element Quotient;

commentator.start ("Rational Minpoly", "Rminpoly");

RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));

std::vector<Integer> F(A.rowdim()+1,1);
std::vector<Integer> M(A.rowdim()+1,1);
std::vector<Integer> Di(A.rowdim());

RationalMatrixFactory<PID_integer,Rationals, DenseMatrix<Rationals > > FA(&A);
Integer da=1, di=1; Integer D=1;
FA.denominator(da);

for (int i=M.size()-2; i >= 0 ; --i) {
	//c[m]=1, c[0]=det(A);
	FA.denominator(di,i);
	D *=di;
	Di[i]=di;
	M[i] = M[i+1]*da;
}
for (int i=0; i < M.size() ; ++i ) {
	gcd(M[i],M[i],D);
}

PID_integer Z;
DenseMatrix<PID_integer> Atilde(Z,A.rowdim(), A.coldim());
FA.makeAtilde(Atilde);

ChineseRemainder< EarlyMultipCRA<Modular<double> > > cra(4UL);
MyRationalModularMinpoly<DenseMatrix<Rationals > , MyMethod> iteration1(A, Met, M);
MyIntegerModularMinpoly<DenseMatrix<PID_integer>, MyMethod> iteration2(Atilde, Met, Di, M);
MyModularMinpoly<MyRationalModularMinpoly<DenseMatrix<Rationals > , MyMethod>, 
		  MyIntegerModularMinpoly<DenseMatrix<PID_integer>, MyMethod> >  iteration(&iteration1,&iteration2);

RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > RR;

std::vector<Integer> PP,PP2; // use of integer due to non genericity of cra. PG 2005-08-04
UserTimer t1,t2;
t1.clear();
t2.clear();
t1.start();
	cra(2,PP,iteration1,genprime);
t1.stop();
t2.start();
	cra(2,PP,iteration2,genprime);
t2.stop();

if (t1.time() < t2.time()) {
	//cout << "ratim";
	iteration.setSwitcher(1);
} else {
	//cout << "intim";
	iteration.setSwitcher(2);
}

int k=2;
while (! cra(k,PP, iteration, genprime)) {
	k *=2;
	Integer m; //Integer r; Integer a,b;
	cra.getModulus(m);
	cra.result(PP);//need to divide
	for (int i=0; i < PP.size(); ++i) {
		Integer D_1;
		inv(D_1,M[i],m);
		PP[i] = (PP[i]*D_1) % m;
	}
	Integer den,den1;
	std::vector<Integer> num(A.rowdim()+1);
	std::vector<Integer> num1(A.rowdim()+1);
	if (RR.reconstructRational(num,den,PP,m,-1)) {//performs reconstruction strating form c[m], use c[i] as prec for c[i-1]
		cra(1,PP,iteration,genprime);
		cra.getModulus(m);
		for (int i=0; i < PP.size(); ++i) {
                        Integer D_1;
                        inv(D_1,M[i],m);
                	PP[i] = (PP[i]*D_1) % m;
                }
		bool terminated = true;
		if (RR.reconstructRational(num1,den1,PP,m,-1)) {
			if (den==den1) {
				for (int i=0; i < num.size(); ++i) {
					if (num[i] != num1[i]) {
						terminated =false;
						break;
					}
				}
			} else {
				terminated = false;
			}
			//set p
			if (terminated) {
				size_t i =0;
				integer t,tt,ttt;
				integer err;
				size_t max_err = 0;
				Quotient qerr;
				p.resize(PP.size());
					typename Vector <typename Rationals::Element>::iterator it;
					Rationals Q;
				for (it= p.begin(); it != p.end(); ++it, ++i) {
			        	A.field().init(*it, num[i],den);
					Q.get_den(t,*it);
					if (it != p.begin()) Q.get_den(tt,*(it-1));
					else tt = 1;
					Q.init(qerr,t,tt);

			        }
				return p;	
				break;
			}
		}
	}
}

cra.result(PP);

size_t i =0;
integer t,tt; 
integer err; 
size_t max_res=0;int max_i; double rel;
size_t max_resu=0; int max_iu;
size_t max_err = 0;
Quotient qerr;
p.resize(PP.size());

typename Vector <typename Rationals::Element>::iterator it;

Rationals Q;
for (it= p.begin(); it != p.end(); ++it, ++i) {
	A.field().init(*it, PP[i],M[i]);
	Q.get_den(t, *it);
	Q.get_num(tt,*it);
	err = M[i]/t;
	size_t resi = err.bitsize() + tt.bitsize() -1;
	size_t resu = t.bitsize() + tt.bitsize() -1;
	if (resi > max_res) {max_res = resi; max_i=i;}
	if (resu > max_resu) {max_resu = resu; max_iu =i;} 
}

max_res=0;
for (it= p.begin()+1; it != p.end(); ++it) {
        Q.get_den(t, *it);
	Q.get_den(tt, *(it-1));
        Q.init(qerr,t,tt);
	Q.get_num(tt, *it);
	size_t resi = Q.bitsize(t,qerr) + tt.bitsize() -2;
	if (resi > max_res) {max_res = resi; max_i=i;}
        }

commentator.stop ("done", NULL, "Iminpoly");

return p;

}

}

#endif //__LINBOX_rat_minpoly_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
