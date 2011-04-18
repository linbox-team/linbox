/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/lifting-container-base.h
 * Copyright (C) 2005  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_ternary_lattice_H
#define __LINBOX_ternary_lattice_H

#include <iostream>
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/field/PID-integer.h>
#include <linbox/util/timer.h>
#include <linbox/integer.h>
#include <algorithm>
#include <vector>

#include <dpe.h>
#include <math.h>



int int_div;
int int_gauss;
int int_comb;


struct StopReduce{};

namespace LinBox
{

	int large_double_division(integer &x, const integer &y, const integer &z)
	{
		double x_m, y_m, z_m;
		long   x_e, y_e, z_e;
		Timer t;

		y_m = mpz_get_d_2exp(&y_e, y.get_mpz());
		z_m = mpz_get_d_2exp(&z_e, z.get_mpz());
		x_e = y_e - z_e;


		if (x_e <53) {
			x_m = y_m / z_m;

			if (x_m == 0.){
				x=0;
				return 0;
			}
			else {
				int e;
				x_m = frexp(x_m,&e);
				x_e +=e;
				x= round(ldexp(x_m,x_e));
				return 0;
				//return x= ldexp(x_m,x_e);
			}
		}
		else {
			int_div++;std::cout<<"Exact Division\n";
			x=y/z;
			return 1;
		}
	}


	class TernaryLattice {
	public:


	protected:

		integer b1[3], b2[3], b3[3];
		integer lb1, lb2, lb3;
		integer b1b2, b1b3, b2b3 ;

		Timer Tgauss, Tcoeff, Tgetcoeff, Tcombine;

		integer y1, y2, TMP;
		integer B1B2LB1, B1B2LB2, B2B3LB2, B1B3LB1;
		integer x10,x11,x12,x20,x21,x22;

		inline void innerProduct(integer &z, const integer x[3] , const integer  y[3]){
			integer::mul(z,x[0],y[0]);
			integer::axpyin(z,x[1],y[1]);
			integer::axpyin(z,x[2],y[2]);
		}

		inline void SquareEuclideanLength(integer& l, const integer y[3]){
			innerProduct(l, y, y);
		}


		inline void EuclideanLength(integer& l, const integer y[3] )
		{
			innerProduct(l, y, y);
			::sqrt(l,l);
		}

		inline void swap(integer x[3], integer y[3] )
		{
			integer tmp;
			tmp=x[0];
			x[0]=y[0];
			y[0]=tmp;
			tmp=x[1];
			x[1]=y[1];
			y[1]=tmp;
			tmp=x[2];
			x[2]=y[2];
			y[2]=tmp;
		}

		inline void maxpyin(integer r[3], const integer &a, const integer y[3])
		{
			if ((a.bitsize()<32) && (a >= 0)){
				unsigned long aa= a;
				mpz_submul_ui(r[0].get_mpz(), y[0].get_mpz(), aa);
				mpz_submul_ui(r[1].get_mpz(), y[1].get_mpz(), aa);
				mpz_submul_ui(r[2].get_mpz(), y[2].get_mpz(), aa);
			}
			else {
				integer::maxpyin(r[0],a,y[0]);
				integer::maxpyin(r[1],a,y[1]);
				integer::maxpyin(r[2],a,y[2]);
			}
		}

		inline void axpyin(integer r[3], const integer &a, const integer y[3])
		{
			if ((a.bitsize()<32) && (a >= 0)){
				unsigned long aa= a;
				mpz_addmul_ui(r[0].get_mpz(), y[0].get_mpz(), aa);
				mpz_addmul_ui(r[1].get_mpz(), y[1].get_mpz(), aa);
				mpz_addmul_ui(r[2].get_mpz(), y[2].get_mpz(), aa);
			}
			else {
				integer::axpyin(r[0],a,y[0]);
				integer::axpyin(r[1],a,y[1]);
				integer::axpyin(r[2],a,y[2]);
			}
		}

		inline void sort()
		{
			if (lb1 > lb2){
				swap(b1,b2);
				std::swap(lb1,lb2);
				std::swap(b1b3,b2b3);
			}
			if (lb1 > lb3){
				swap(b1,b3);
				std::swap(lb1,lb3);
				std::swap(b1b2,b2b3);
			}
			if (lb2 > lb3){
				swap(b2,b3);
				std::swap(lb2,lb3);
				std::swap(b1b2, b1b3);
			}
		}

		inline void binaryGaussReduce()
		{
			int_gauss++;
			integer  r;
			r=b1b2/lb1;
			integer a[3], la;
			assign(a,b2);
			maxpyin(a,r,b1);
			la= (-r*(b1b2))<<1;
			integer::axpyin(la, r*r, lb1);
			la+=lb2;

			if (la < lb1){
				assign(b2,b1);lb2=lb1;
				assign(b1,a);lb1=la;
				std::swap(b2b3,b1b3);
				integer::maxpyin(b1b2,r,lb2);
				integer::maxpyin(b1b3,r, b2b3);
				binaryGaussReduce();
			}else{
				if (la < lb2){
					assign(b2,a);lb2=la;
					integer::maxpyin(b2b3,r, b1b3);
					integer::maxpyin(b1b2,r,lb1);
				}
			}
		}

		inline void assign(integer x[3], const integer y[3])
		{
			x[0]=y[0];
			x[1]=y[1];
			x[2]=y[2];
		}


		inline integer SEL(integer &l, const integer &x1, const integer &y1, const integer &y2, const integer &y3)
		{
			Timer t;
			t.start();
			integer tmp, tmp_min, y_min;
			integer b2b3_2   = b2b3<<1;
			integer rr = b2b3_2;
			integer::axpyin(rr,x1<<1,b1b2);
			tmp= b1b3<<1;
			l=lb3;
			integer::axpyin(tmp,x1,lb1);
			integer::axpyin(l, tmp, x1);

			tmp = rr;
			integer::axpyin(tmp, y1, lb2);
			tmp*=y1;
			tmp_min=tmp;
			y_min  =y1;

			tmp = rr;
			integer::axpyin(tmp, y2, lb2);
			tmp*=y2;
			if (tmp < tmp_min){
				tmp_min=tmp;
				y_min=y2;
			}

			tmp = rr;
			integer::axpyin(tmp, y3, lb2);
			tmp*=y3;
			if (tmp < tmp_min){
				tmp_min=tmp;
				y_min=y3;
			}


			l+=tmp_min;
			t.stop();
			Tgetcoeff+=t;
			return y_min;
		}


	public:

		TernaryLattice(const std::vector<integer> &L)
		{
			linbox_check(L.size()==9);
			b1[0]=L[0];
			b1[1]=L[1];
			b1[2]=L[2];
			b2[0]=L[3];
			b2[1]=L[4];
			b2[2]=L[5];
			b3[0]=L[6];
			b3[1]=L[7];
			b3[2]=L[8];
			//print();
			SquareEuclideanLength(lb1, b1);
			SquareEuclideanLength(lb2, b2);
			SquareEuclideanLength(lb3, b3);
			sort();
			innerProduct(b1b2, b1, b2);
			innerProduct(b1b3, b1, b3);
			innerProduct(b2b3, b2, b3);
		}


		template<class Blackbox>
		TernaryLattice(const Blackbox &L)
		{
			int_div=0;int_gauss=0;int_comb=0;

			Tgauss.clear();
			Tcoeff.clear();
			Tgetcoeff.clear();
			Tcombine.clear();
			L.field().convert(b1[0], L.getEntry(0,0));
			L.field().convert(b1[1], L.getEntry(0,1));
			L.field().convert(b1[2], L.getEntry(0,2));
			L.field().convert(b2[0], L.getEntry(1,0));
			L.field().convert(b2[1], L.getEntry(1,1));
			L.field().convert(b2[2], L.getEntry(1,2));
			L.field().convert(b3[0], L.getEntry(2,0));
			L.field().convert(b3[1], L.getEntry(2,1));
			L.field().convert(b3[2], L.getEntry(2,2));
			SquareEuclideanLength(lb1, b1);
			SquareEuclideanLength(lb2, b2);
			SquareEuclideanLength(lb3, b3);
			sort();
			innerProduct(b1b2, b1, b2);
			innerProduct(b1b3, b1, b3);
			innerProduct(b2b3, b2, b3);
		}


		void reduce()
		{

			integer lmin_a, x1, x2;
			{
				Timer t;

				/*
				   std::cout<<"Calling reduce (";
				   std::cout<<lb1.bitsize()<<",";
				   std::cout<<lb2.bitsize()<<",";
				   std::cout<<lb3.bitsize()<<")";;
				   */
				//std::cout<<lb1<<",";
				//std::cout<<lb2<<",";
				//std::cout<<lb3;
				//std::cout<<")\n";
				//print();

				t.start();
				binaryGaussReduce();//std::cout<<"Gauss reduce done\n";
				t.stop();
				Tgauss+=t;
				t.clear();
				t.start();


				integer::div(B1B2LB1, b1b2, lb1);
				integer::div(B1B2LB2, b1b2, lb2);
				integer::div(B2B3LB2, b2b3, lb2);
				integer::div(B1B3LB1, b1b3, lb1);

				TMP= integer(integer(1));
				integer::maxpyin(TMP, B1B2LB1, B1B2LB2);


				y2= B2B3LB2;
				integer::maxpyin(y2, B1B2LB2, B1B3LB1);
				integer::divin(y2,TMP);
				integer::negin(y2);

				y1= B1B3LB1;
				integer::maxpyin(y1, B1B2LB1, B2B3LB2);
				integer::divin(y1,TMP);
				integer::negin(y1);

				x10= y1;
				x11= 1+y1;
				x12= -1+y1;
				x20= y2;
				x21= 1+y2;
				x22= -1+y2;

				t.stop();
				Tcoeff+=t;

				integer la, x_tmp;


				x2   = SEL(la, x10, x20, x21, x22); lmin_a=la;x1=x10;
				x_tmp = SEL(la, x11, x20, x21, x22); if (la < lmin_a) {lmin_a=la;x1=x11;x2=x_tmp;}
				x_tmp = SEL(la, x12, x20, x21, x22); if (la < lmin_a) {lmin_a=la;x1=x12;x2=x_tmp;}


				//if (::abs(x1)-::abs(y1)  != 1) std::cout<<"|x1|-|y1|: "<< ::abs(x1)-::abs(y1)<<"\n";
				//if (::abs(x2)-::abs(y2)  != 1) std::cout<<"|x2|-|y2|: "<< ::abs(x2)-::abs(y2)<<"\n";

			}


			if (lmin_a < lb3){
				int_comb++;
				Timer tt;
				tt.clear();
				tt.start();
				axpyin(b3, x1, b1);
				axpyin(b3, x2, b2);
				lb3=lmin_a;
				integer::axpyin(b1b3,x1,lb1);
				integer::axpyin(b1b3,x2,b1b2);
				integer::axpyin(b2b3,x1,b1b2);
				integer::axpyin(b2b3,x2,lb2);

				if (lb3 < lb1){
					swap(b1, b3);std::swap(lb1,lb3);std::swap(b1b2, b2b3);
				}
				if (lb3 < lb2){
					swap(b2, b3);std::swap(lb2,lb3);std::swap(b1b3, b1b2);
				}

				//sort();
				tt.stop();
				Tcombine+=tt;
				reduce();
			}
		}

		integer* operator[](size_t i)
		{
			if (i==0) return b1;
			if (i==1) return b2;
			if (i==2) return b3;
		}


		void timing()
		{
			std::cout<<"Gauss reduce     : "<<Tgauss<<" with "<<int_gauss<<" calls \n";
			std::cout<<"Coeff computation: "<<Tcoeff<<" with "<<int_div<<" exact division \n";
			std::cout<<"       get coeff : "<<Tgetcoeff<<"\n";
			std::cout<<"Combine row      : "<<Tcombine<<" with "<<int_comb<<" combinaison \n";
		}

		void print()
		{
			PID_integer Z;
			BlasBlackbox<PID_integer> M(Z,3,3);
			M.setEntry(0,0,b1[0]);
			M.setEntry(0,1,b1[1]);
			M.setEntry(0,2,b1[2]);
			M.setEntry(1,0,b2[0]);
			M.setEntry(1,1,b2[1]);
			M.setEntry(1,2,b2[2]);
			M.setEntry(2,0,b3[0]);
			M.setEntry(2,1,b3[1]);
			M.setEntry(2,2,b3[2]);
			//M.write(std::cout);
		}

		template<class Blackbox>
		void getLattice(Blackbox &M)
		{
			M.setEntry(0,0,b1[0]);
			M.setEntry(0,1,b1[1]);
			M.setEntry(0,2,b1[2]);
			M.setEntry(1,0,b2[0]);
			M.setEntry(1,1,b2[1]);
			M.setEntry(1,2,b2[2]);
			M.setEntry(2,0,b3[0]);
			M.setEntry(2,1,b3[1]);
			M.setEntry(2,2,b3[2]);
		}

	};

	class LargeDouble{
	protected:
		double _m;
		long   _e;
	public:
		LargeDouble(const integer &x)
		{
			_m=mpz_get_d_2exp(&_e, x.get_mpz());
		}

		LargeDouble() {}

		LargeDouble(const LargeDouble &x) :
			_m(x._m), _e(x._e)
		{}

		LargeDouble& operator= (const LargeDouble &x)
		{
			_m = x._m;
			_e = x._e;
			return *this;
		}

		integer& convert(integer &x)
		{
			if (_e < 53 ){
				if (_m == 0.)
					return x=integer(0);
				else
					return x=integer(round(ldexp(_m,_e)));
			}
			else{
				x = (_m* 0x1P53);//9007199254740992;
				x = x<<(_e-53);
				return x;
			}
		}

		static LargeDouble& div(LargeDouble &x, const LargeDouble &y, const LargeDouble &z)
		{
			x._e = y._e - z._e;
			x._m = y._m / y._m;
			int e;
			x._m = frexp(x._m, &e);
			x._e+=e;
			//if (x._e >= 53) throw StopReduce();
			return x;
		}

		static LargeDouble& divin(LargeDouble &x, const LargeDouble &y)
		{
			x._e = x._e - y._e;
			x._m = x._m / y._m;
			int e;
			x._m = frexp(x._m, &e);
			x._e+=e;
			//if (x._e >= 53) throw StopReduce();
			return x;
		}

		static LargeDouble& axpyin(LargeDouble &x, const LargeDouble &a, const LargeDouble &y)
		{
			long e = a._e + y._e;
			if (x._e > e + 53)
				return x;
			else {
				double m = a._m * y._m;
				if (e > x._e + 53){
					x._e = e;
					x._m = m;
					return x;
				}
				else {
					long ee = x._e - e;
					if (ee >=0){
						x._m += ldexp( m, -ee);
					}
					else {
						x._m = m + ldexp( x._m, ee);
						x._e = e;
					}
					int t;
					x._m = frexp(x._m, &t);
					x._e+=t;
					return x;
				}
			}
		}

		static LargeDouble& maxpyin(LargeDouble &x, const LargeDouble &a, const LargeDouble &y)
		{

			long e = a._e + y._e;
			if (x._e > e + 53)
				return x;
			else {
				double m = a._m * y._m;
				if (e > x._e + 53){
					x._e = e;
					x._m = - m;
					return x;
				}
				else {
					long ee = x._e - e;
					if (ee >=0){
						x._m -= ldexp( m, -ee);
					}
					else {
						x._m = ldexp( x._m, ee) -m ;
						x._e = e;
					}
					int t;
					x._m = frexp(x._m, &t);
					x._e+=t;
					return x;
				}
			}
		}

		static LargeDouble& mul(LargeDouble &x, const LargeDouble &y, const LargeDouble &z)
		{
			x._e = y._e + z._e;
			x._m = y._m + z._m;
			int t;
			x._m = frexp(x._m, &t);
			x._e+=t;
			return x;
		}

		static LargeDouble& negin(LargeDouble &x)
		{
			x._m = -x._m;
			return x;
		}
	};


	std::ostream& operator<< (std::ostream& os, LargeDouble& x)
	{
		integer tmp;
		x.convert(tmp);
		return os<<tmp;
	}


}// end of namespace LinBox

#endif //__LINBOX_ternary_lattice_H

