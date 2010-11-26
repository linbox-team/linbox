/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/** 
 * examples/qchar.C
 *
 * Copyright (C) 2007, 2010 A. Urbanska
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */
#include <linbox/field/PID-integer.h>
#include <linbox/field/modular-double.h>
#include <linbox/field/gmp-rational.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/sparse.h>
#include "linbox/solutions/charpoly.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/ring/givaro-polynomial.h"
#include "linbox/element/givaro-polynomial.h"

using namespace LinBox;
using namespace std;

typedef PID_integer Integers;
typedef Integers::Element Integer;
typedef Modular<double > Field;
typedef Field::Element Element;
typedef Vector<Integers>::Dense DVector;
typedef GMPRationalField Rationals;
typedef Rationals::Element Quotient;
typedef DenseMatrix<Integers > Blackbox;
typedef SparseMatrix<Rationals > RBlackbox;

//#define CONT_FR

typedef UserTimer myTimer;

myTimer tRationalModular;

template <class Field, class Polynomial>
std::ostream& printPolynomial (std::ostream& out, const Field &F, const Polynomial &v) 
{
	for (int i = v.size () - 1; i >= 0; i--) {
		F.write (out, v[i]);
		if (i > 0)
			out << " X^" << i << " + ";
		out << "\n";
	}
	return out;
}

template <class Blackbox, class MyMethod>
struct PrecRationalModularCharpoly {       
  const Blackbox &A;
  const MyMethod &M;
  const Integer &prec;
	
PrecRationalModularCharpoly(const Integer& detPrec, const Blackbox& b, const MyMethod& n) 
		: prec(detPrec), A(b), M(n) {}
	
  template<typename Polynomial, typename Field>
  Polynomial& operator()(Polynomial& P, const Field& F) const {
  
    typedef typename Blackbox::template rebind<Field>::other FBlackbox;
    FBlackbox * Ap;
    
tRationalModular.clear();
tRationalModular.start();
    MatrixHom::map(Ap, A, F);
tRationalModular.stop();
    //minpoly(P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
    charpoly( P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
    //minpolySymmetric(P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
    integer p;
    F.characteristic(p);
    cout<<"Valence "<<p<<" = P[P.size()-1]" << P[P.size(0)-1] << " R-M time " ;
    tRationalModular.print(cout);
    cout << "\n" ;
     
    typename Field::Element fprec;
    F.init(fprec, prec);
    //printPolynomial (std::cout, F, P); 
    delete Ap;

    for (int i=0; i < P.size(); ++i)
      F.mulin(P[i],fprec);
    //std::cout<<"prec*Det(A) mod "<<p<<" = P[0]" << P[0] << "\n";
    return P;
  }            
};

template <class Blackbox, class MyMethod>
struct PrecRationalModularMinpoly {
  const Blackbox &A;
  const MyMethod &M;
  const DVector &prec;


PrecRationalModularMinpoly(const DVector& detPrec, const Blackbox& b, const MyMethod& n)
                : prec(detPrec), A(b), M(n) {}

  template<typename Polynomial, typename Field>
  Polynomial& operator()(Polynomial& P, const Field& F) const {

    typedef typename Blackbox::template rebind<Field>::other FBlackbox;
    FBlackbox * Ap;

tRationalModular.clear();
tRationalModular.start();
    MatrixHom::map(Ap, A, F);
tRationalModular.stop();
    //minpoly(P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
    //charpoly( P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
    minpolySymmetric(P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
    integer p;
    F.characteristic(p);
    //cout<<"Det(A) mod "<<p<<" = P[0]" << P[0] << " R-M time " ;
    cout<<"Valence "<<p<<" = P[P.size()-1]" << P[P.size()-1] << " R-M time " ;
    tRationalModular.print(cout);
    cout << "\n" ;

    typename Field::Element fprec;
    //F.init(fprec, prec);
    //printPolynomial (std::cout, F, P);
    delete Ap;

	for (int i=0; i < P.size()-1; ++i) {
		F.init(fprec, prec[i]);
		F.mulin(P[i],fprec);
	}
    //std::cout<<"prec*Det(A) mod "<<p<<" = P[0]" << P[0] << "\n";
	//printPolynomial(cout,F,P);
    return P;
  }
};


void continuedFractionIn(double x,Integer& num, Integer& den, double epsi,const size_t s, Integer den_bound);
template <class RMatrix>
void generate_precRatMat(string& filename, RMatrix& M, DVector& den, Integer& denPrec);
void i_vti_v(RBlackbox& Res, DenseMatrix<Rationals >& M, DVector& den, Integer& detPrec);

int main (int argc, char** argv)
{

if ((argc < 3) || (argc > 4) ) {
  cout << "Usage: qchar n matrix" << endl; 
	return -1;
}

size_t n = atoi(argv[1]);
cout << "Matrix size " << n<< endl;
string filename;
filename=argv[2];

Integers Z;
Rationals Q;

Integer detPrec = 1;
Integer detNum;
Integer detDen;

Blackbox DM(Z,n,n);
DenseMatrix<Rationals > In(Q,n,n);
RBlackbox M(Q,n,n);
DVector V(n,1); 

//	generate_precRatMat(filename, In, V, detPrec);
//i_vti_v(M,In, V, detPrec);
  generate_precRatMat(filename, M, V, detPrec);
  cout << "detPrec is " << detPrec  << "\n";
  typedef GivPolynomialRing<Integers,Dense> IntPolRing;
  IntPolRing::Element c_A;

		Integer max = 1,min=0;
		for (int i=0; i < n; ++i) {
		  for (int j=0; j < n; ++j) {
		    Integer a;
		    Q.convert(a,M.getEntry(i,j));
		  	//      cerr<<"it="<<(*it)<<endl;
		    if (max < a)
				max = a; 
		    if (min > a)
				min = a;
		  }
		}
		if (max<-min) 
			max=-min;
		else max = max+1;
		double hadamarcp = (double)n/2.0*(log(double(n))+2*log(double(max))+0.21163275)/log(2.0);

		
		cout << "had" << hadamarcp << "\n";
		cout << "had2" << (Integer)hadamarcp*detPrec << "\n";

  RandomPrimeIterator genprime( 26-(int)ceil(log((double)M.rowdim())*0.7213475205)); 
  ChineseRemainder< EarlyMultipCRA<Field  > > cra(3UL);
  typedef Method::Hybrid MyMethod; 
  MyMethod Met;
	//PrecRationalModularCharpoly <RBlackbox  ,MyMethod> iteration (detPrec, M, Met);
  PrecRationalModularMinpoly< RBlackbox  ,MyMethod> iteration(V, M, Met);
  cra (c_A, iteration, genprime);
  printPolynomial(cout, Z, c_A);
}

/* for a rational number num/den construct a continued fraction approximation:
    - with den bounded by den_bound
AND - number of steps bounded by s
NEW fraction replaces num/den
*/

void continuedFractionIn(Integer& num, Integer& den, double epsi,const size_t s, Integer den_bound) 
{
  //den_bound = 100;
  Integer a=num;
  Integer b=den;
  Integer t;
  //double y;
  Integer f[s];
  f[0] = a/b;// y = f[0];
  int i=1;
  while (1) {//abs(y-x)>=epsi) {
    //cout << "a :" << a << " b: " << b << "\n"; 
    t = a%b;
    a = b;
    b = t;
    f[i] = a/b;
    //y = (double)f[i];
    num = 1;
    den = f[i];
    for (int j=i-1; j > 0; --j) {
      Integer tmp = num;
      num = den;
      den = f[j]*den+tmp;
      //y = (double)f[j]+1/y;
    }
    //y = (double)f[0]+1/y;
    num = den*f[0]+num;
    if (den >= den_bound) break;
    ++i;
    if (i >= s) {
      break;
    }
  }
}

void i_vti_v(RBlackbox& Res, DenseMatrix<Rationals >& M, DVector& den, Integer& detPrec)
{
	detPrec=1;
	Rationals Q;
	int n = M.coldim();
	//DVector denV(n,1);
	for (int i=0; i < den.size();++i) den[i]=1;
	for (int i=0; i < n; ++i) {
		for (int j=0; j <=i ; ++j) {
			Integer q_num =0;
			Integer q_den =1;
			for (int k=0; k < n; ++k) {
				Integer deno_ik, nume_ik, deno_jk, nume_jk, deno, nume;
				Q.get_den (deno_ik, M.getEntry(k,i));
				Q.get_num (nume_ik, M.getEntry(k,i));
				Q.get_den (deno_jk, M.getEntry(k,j));
				Q.get_num (nume_jk, M.getEntry(k,j));
				
				if (i==k) {
					nume_ik = deno_ik-nume_ik;	
				} else {
					nume_ik = -nume_ik;
				}

				if (j==k) {
					nume_jk = deno_jk-nume_jk;
				} else {
					nume_jk = -nume_jk;
				}
				//cout << nume_ik << nume_jk;

			        deno = deno_ik*deno_jk;
				nume = nume_ik*nume_jk;
				//cout << q_num << "/" << q_den << " " ;
				//cout << nume << "/" << deno << " " ;
				if (deno==q_den) q_num+=nume;
				else {
					if (nume != 0) {
						Integer g = gcd(q_den, deno);
						q_num = q_num*deno/g+nume*q_den/g;
						q_den = q_den*deno/g;
					}
				}
				//cout << q_num << "/" << q_den << " " ;
			}
			if (q_num != 0) {
				Quotient q(q_num,q_den);
				if (i!=j) {
					den[i]=lcm(den[i], q_den);
					den[j]=lcm(den[j], q_den);
					Res.setEntry(i,j,q);
                                	Res.setEntry(j,i,q);
					cout << i << " " << j << " " << q_num << "/" << q_den << "\n";
					cout << j << " " << i << " " << q_num << "/" << q_den << "\n";
				} else {
					den[i]=lcm(den[i], q_den);
					Res.setEntry(i,j,q);
					cout << i << " " << j << " " << q_num << "/" << q_den << "\n";
				}
			}
		}
	}
	Integer d = 1;
	for (int i=0; i < den.size(); ++i) {
		detPrec *=den[i];
        	d = lcm(d, den[i]);
	}
        den[den.size()-1]=d;
        for (int i=den.size()-2; i >=0; --i) {
                den[i]=d*den[i+1];
        }
	
}

template <class RMatrix>
void generate_precRatMat(string& filename, RMatrix& M, DVector& den, Integer& denPrec)
{
int prec=10;
  denPrec = 1;
  ifstream is(filename.c_str(), std::ios::in);
  int m,n;
  char c;
  is >> m >> n;
  do is >> c; while (isspace (c));
  if ((m > M.rowdim()) || (n > M.coldim())) {
    cout << m << " " << n << " " ;
    cout << "Wrong matrix size\n";
    return;
  }

  cout << "Reading matrix\n";

  den.resize(m,1);
  while (1) {//temp fix 
	  //while (is >> m) {//temp fix
	  //if (m==0) break;//temp fix
    is >> m;//temp fix
    is >> n;
    //cout << m << n << "\n";
    //    long double x;
    //    is >> x;
    char xstr[500];
    char numstr[500];
    //    cout << x << " ";
    //    sprintf(xstr, "%100f", x);
    //    cout << xstr << " " ;
    //    is >> c; int i=0;
    //while (c!= '\n') {
    //  xstr[i]=c;
    //  ++i;
    //  is >> c;
    //  if (i>=199) {
    //	xstr[i]=0;
    //	break;
    //      }
    //    }
    //is >> c;
    is.getline(xstr,500);

    Integer ten=1;
    if (strchr(xstr, '/') != NULL) {
	    //cout << "xstr " << xstr << "\n";
	    strcpy(numstr, strtok(xstr, "/"));
	    char tenstr[500];
	    strcpy(tenstr, strtok(NULL, "/"));
	    if (prec < strlen(tenstr)) {
		    int cut = strlen(tenstr)-prec;
		    tenstr[prec]=0;
		    //c = numstr[strlen(numstr)-2];
		    numstr[strlen(numstr)-cut]=0;//temp round down
	    }
	    Integer tmp(tenstr);
	    ten = tmp;
	    //cout << numstr << "/"  << ten << "\n";
    }
    else 
    {

    //    cout << "xstr" << xstr << " " ;
    int i=0;
    int j=0;
    c=xstr[i]; ++i;
    while (isspace (c)) {
	if (i >=strlen(xstr)) break;
      c=xstr[i];
      ++i;
    }
	if (c=='-') {
		numstr[j]=c;++j;
        	if (i < strlen(xstr)) {
                	c = xstr[i];
                	++i;
		}
	}
	if (c=='0') {
	  //cout << "zero";
		if (i < strlen(xstr)) {
      		c = xstr[i];
      		++i;
	} 
    } //else cout << " not 0" << c; 
    while (strchr("0123456789",c) != NULL) {
      //cout << "number";
      numstr[j]=c; ++j;
	if (i >=strlen(xstr)) break;
      c=xstr[i];
      ++i; 
    }
    if (c=='.') {
    	if (i < strlen(xstr)) {
      		c = xstr[i];
      	//cout << "c" << c ;
      		++i;
      		while (strchr("0123456789",c) != NULL) {
			numstr[j]=c;++j;
			if (i >=strlen(xstr)) break;
			c=xstr[i];
			++i; 
		}	
	} else {
                cout << "wrong number format at .";
		break;
      	}     
    }
    numstr[j]=0;
    //cout << "num" << numstr << " " ;
    
    size_t dplaces=0;
    int exp=0;
    j=0;
    c = xstr[j]; ++j;
    while (c != '.') {
      if (j >= strlen(xstr)) break;
      c = xstr[j];
	++j;
    }
    if (j < strlen(xstr)) {
	    c = xstr[j]; ++j;
	    while (c != 'e') {
		    ++dplaces;
		    ten *= 10; 
		    if (j >= strlen(xstr)) break;
		    c = xstr[j];
			++j;
	    }
    }
	if (j < strlen(xstr)) {	
    		sscanf(xstr+j, "%d", &exp);
    		j=j-2;c = xstr[j];
    		while (c=='0') {
      			ten/=10;
      			--j;
     		 	c = xstr[j];
		}
    	}

    dplaces -=exp;
    if (exp > 0) {
      for (int i=0; i < exp; ++i) ten /=10;
    } else if (exp < 0) { 
      for (int i=0; i < exp; ++i) ten *=10;
    }
    //double nume = x * (double)ten;
    }
    Integer num (numstr); 
    Integer g;
 
#ifdef CONT_FR
    double epsi = 1/(double)ten*10;
    size_t s=10;
    continuedFractionIn(num, ten, epsi, s, ten);
#endif			
    //cout << m << " " << n << " " << num << "/" << ten << "\n"; 
    g = gcd(num,ten);
    //g =1;
    Quotient q(num/g,ten/g);
    if ((m==0)&&(n==0)&&(num==0)) break;//temp fix
    //M.setEntry(m-1,n-1,q);//temp fix
    M.setEntry(m,n,q);
    den[m-1] = lcm(den[m-1],ten/g);
  }
	Integer d = 1;
  for (int i=0; i < den.size(); ++i) {
    denPrec *=den[i];
	d = lcm(d, den[i]);
  }
	den[den.size()-1]=d;
	cout << d << "\n";
	for (int i=den.size()-2; i >=0; --i) {
		den[i]=d*den[i+1];
		cout << den[i] << "\n";
	}
}

