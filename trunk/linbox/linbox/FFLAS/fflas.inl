/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/FFLAS/fflas.inl
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
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
 *
 */

#ifndef __FFLAS_INL
#define __FFLAS_INL


#ifdef BLAS_AVAILABLE
extern "C" {
#include "cblas.h"
}
#endif

#include <math.h>

#include "linbox/FFLAS/FFFMMBLAS.h"
#include "linbox/integer.h"
#include "linbox/util/field-axpy.h"

namespace LinBox {


#ifdef BLAS_AVAILABLE
	
	template <class Field>
	inline typename Field::Element*  Field_dgemm	(const Field& F,
							 int m, int n, int k,
							 int alpha,
							 typename Field::Element * B,
							 int ldb,
							 typename Field::Element * C,
							 int ldc,
							 int beta,
							 typename Field::Element * A,
							 int lda,
							 int nbe=0) 
		
	{ return FFFMMBLAS() (F,m,n,k,alpha,B,ldb,C,ldc,beta,A,lda,nbe);}
	

	// A=B.T^(-1)
	// T is an upper triangular matrix with unit on diagonal
	template <class Field>
	inline void Field_trsm_up_right_unit (const Field& F,
					      int m, int n,
					      typename Field::Element * B, int ldb,
					      typename Field::Element * T, int ldt,
					      typename Field::Element * A, int lda) {
		
		LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(n)*pow(integer(2),n-1)*pow(p,n);
		LinBox::integer norm=pow(p,n);
		
		if (norm > max_double )
			{
				int n2= n >>1;
				int n1= n-n2;
				
				Field_trsm_up_right_unit (F,m,n1,B,ldb,T,ldt,A,lda);
				Field_dgemm (F,m,n2,n1,-1,A,lda,T+n1,ldt,1,B+n1,ldb);
				Field_trsm_up_right_unit (F,m,n2,B+n1,ldb,T+n1*ldt+n1,ldt,A+n1,lda);		       
			}
		else {		
			double Td[n*n];
			double Bd[m*n];			
			
			MatGFq2MatDouble_Triangular (F,n,n,Td,n,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);			
			cblas_dtrsm(CblasRowMajor,CblasRight,CblasUpper, CblasNoTrans,CblasUnit,m,n,1.0,Td,n,Bd,n);
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);			
		}
	}
		
	// A=B.T^(-1)
	// T is an upper triangular matrix. 
	// be careful , at the end of this function T is modified (it should not be modified) !!!
	template <class Field>
	inline void Field_trsm_up_right (const Field& F,
					 int m, int n,
					 typename Field::Element * B, int ldb,
					 typename Field::Element * T, int ldt,
					 typename Field::Element * A, int lda) {
		
		LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(n)*pow(integer(2),n-1)*pow(p,n);
		LinBox::integer norm=pow(p,n);

		if (norm > max_double )
			{
				int n2= n >>1;
				int n1= n-n2;
						
				Field_trsm_up_right (F,m,n1,B,ldb,T,ldt,A,lda);
				Field_dgemm (F,m,n2,n1,-1,A,lda,T+n1,ldt,1,B+n1,ldb);
				Field_trsm_up_right (F,m,n2,B+n1,ldb,T+n1*ldt+n1,ldt,A+n1,lda);		       
			}
		else {		
		
			double Td[n*n];
			double Bd[m*n];

			typedef typename Field::Element Element;	
			Element tmp[n];	 
			// normalisation of T
			for (int i=0;i<n;i++) {
				F.inv(tmp[i],*(T+i*ldt+i));
				for (int j=i;j<n;j++)
					F.mulin(*(T+j+ldt*i),tmp[i]);
			}
				
			MatGFq2MatDouble_Triangular (F,n,n,Td,n,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);
			cblas_dtrsm(CblasRowMajor,CblasRight,CblasUpper, CblasNoTrans,CblasUnit,m,n,1.0,Td,n,Bd,n);
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);
				
			// denormalisation of the result A
			for (int i=0;i<m;i++) 
				for (int j=0;j<n;j++)
					F.mulin(*(A+i*lda+j),tmp[j]);		
		}

	}
	
	// A=T^(-1).B
	// T is a lower triangular matrix with unit on diagonal
	template <class Field>
	void Field_trsm_low_left_unit (const Field& F,
				       int m, int n,
				       typename Field::Element * B, int ldb,
				       typename Field::Element * T, int ldt,
				       typename Field::Element * A, int lda) {
		
		LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(m)*pow(integer(2),m-1)*pow(p,m);
		LinBox::integer norm=pow(p,m);

		if (norm > max_double )
			{
				int m2= m >>1;
				int m1= m-m2;
						
				Field_trsm_low_left_unit (F,m1,n,B,ldb,T,ldt,A,lda);
				Field_dgemm (F,m2,n,m1,-1,T+m1*ldt,ldt,A,lda,1,B+m1*ldb,ldb);
				Field_trsm_low_left_unit (F,m2,n,B+m1*ldb,ldb,T+m1*ldt+m1,ldt,A+m1*lda,lda);		       
			}
		else {		
		
			double Td[m*m];
			double Bd[m*n];		
						
			MatGFq2MatDouble_Triangular_Low (F,m,m,Td,m,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);		
			cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower, CblasNoTrans,CblasUnit,m,n,1.0,Td,m,Bd,n);			
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);
		}
	}
	
	// A=T^(-1).B
	// T is a lower triangular matrix.
	// Be careful , T is modified at the end of the function
	template <class Field>
	void Field_trsm_low_left (const Field& F,
				  int m, int n,
				  typename Field::Element * B, int ldb,
				  typename Field::Element * T, int ldt,
				  typename Field::Element * A, int lda) {
		
		LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(m)*pow(integer(2),m-1)*pow(p,m);
		LinBox::integer norm=pow(p,m);

		if (norm > max_double )
			{
				int m2= m >>1;
				int m1= m-m2;
						
				Field_trsm_low_left (F,m1,n,B,ldb,T,ldt,A,lda);
				Field_dgemm (F,m2,n,m1,-1,T+m1*ldt,ldt,A,lda,1,B+m1*ldb,ldb);
				Field_trsm_low_left (F,m2,n,B+m1*ldb,ldb,T+m1*ldt+m1,ldt,A+m1*lda,lda);		       
			}
		else {		
		
			double Td[m*m];
			double Bd[m*n];		
						
			typedef typename Field::Element Element;	
			Element tmp[m];	 
			// normalisation of T
			for (int i=0;i<m;i++) {
				F.inv(tmp[i],*(T+i*ldt+i));
				for (int j=i;j<m;j++)
					F.mulin(*(T+j+ldt*i),tmp[j]);
			}

			MatGFq2MatDouble_Triangular_Low (F,m,m,Td,m,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);		
			cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower, CblasNoTrans,CblasUnit,m,n,1.0,Td,m,Bd,n);			
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);

			// denormalisation of the result A
			for (int i=0;i<m;i++) 
				for (int j=0;j<n;j++)
					F.mulin(*(A+i*lda+j),tmp[i]);
		}
	}

	// A=B.T^(-1)
	// T is a lower triangular matrix with unit on diagonal
	template <class Field>
	inline void Field_trsm_low_right_unit (const Field& F,
					       int m, int n,
					       typename Field::Element * B, int ldb,
					       typename Field::Element * T, int ldt,
					       typename Field::Element * A, int lda) {
		
		LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(n)*pow(integer(2),n-1)*pow(p,n);
		LinBox::integer norm=pow(p,n);
		
		if (norm > max_double )
			{
				int n2= n >>1;
				int n1= n-n2;
				
				Field_trsm_low_right_unit (F,m,n2,B+n1,ldb,T+n1*ldt+n1,ldt,A+n1,lda);
				Field_dgemm (F,m,n1,n2,-1,A+n1,lda,T+n1*ldt,ldt,1,B,ldb);
				Field_trsm_low_right_unit (F,m,n1,B,ldb,T,ldt,A,lda);		       
			}
		else {		
			double Td[n*n];
			double Bd[m*n];			
			
			MatGFq2MatDouble_Triangular_Low (F,n,n,Td,n,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);			
			cblas_dtrsm(CblasRowMajor,CblasRight,CblasLower, CblasNoTrans,CblasUnit,m,n,1.0,Td,n,Bd,n);
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);			
		}
	}

	// A=B.T^(-1)
	// T is a lower triangular matrix.
	template <class Field>
	inline void Field_trsm_low_right (const Field& F,
					  int m, int n,
					  typename Field::Element * B, int ldb,
					  typename Field::Element * T, int ldt,
					  typename Field::Element * A, int lda) {
		
		LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(n)*pow(integer(2),n-1)*pow(p,n);
		LinBox::integer norm=pow(p,n);
		
		if (norm > max_double )
			{
				int n2= n >>1;
				int n1= n-n2;
				
				Field_trsm_low_right (F,m,n2,B+n1,ldb,T+n1*ldt+n2,ldt,A+n1,lda);
				Field_dgemm (F,m,n1,n2,-1,A+n1,lda,T+n1*ldt,ldt,1,B,ldb);
				Field_trsm_low_right (F,m,n1,B,ldb,T,ldt,A,lda);		       
			}
		else {		
			double Td[n*n];
			double Bd[m*n];	
		
			typedef typename Field::Element Element;	
			Element tmp[n];	 
			// normalisation of T
			for (int i=0;i<n;i++) {
				F.inv(tmp[i],*(T+i*ldt+i));
				for (int j=i;j>=0;j--)
					F.mulin(*(T+j+ldt*i),tmp[i]);
			}
			
			MatGFq2MatDouble_Triangular_Low (F,n,n,Td,n,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);			
			cblas_dtrsm(CblasRowMajor,CblasRight,CblasLower, CblasNoTrans,CblasUnit,m,n,1.0,Td,n,Bd,n);
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);	

			// denormalisation of the result A
			for (int i=0;i<m;i++) 
				for (int j=0;j<n;j++)
					F.mulin(*(A+i*lda+j),tmp[j]);		
		}
	}
	

	// A=T^(-1).B
	// T is an upper triangular matrix with unit on diagonal
	template <class Field>
	void Field_trsm_up_left_unit (const Field& F,
				      int m, int n,
				      typename Field::Element * B, int ldb,
				      typename Field::Element * T, int ldt,
				      typename Field::Element * A, int lda) {
		
		LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(m)*pow(integer(2),m-1)*pow(p,m);
		LinBox::integer norm=pow(p,m);

		if (norm > max_double )
			{
				int m2= m >>1;
				int m1= m-m2;
				
				Field_trsm_up_left_unit (F,m2,n,B+m1*ldb,ldb,T+m1*ldt+m1,ldt,A+m1*lda,lda);
				Field_dgemm (F,m1,n,m2,-1,T+m1,ldt,A+m1*lda,lda,1,B,ldb);
				Field_trsm_up_left_unit (F,m1,n,B,ldb,T,ldt,A,lda);		       
			}
		else {		
			
			double Td[m*m];
			double Bd[m*n];		
						
			MatGFq2MatDouble_Triangular (F,m,m,Td,m,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);		
			cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper, CblasNoTrans,CblasUnit,m,n,1.0,Td,m,Bd,n);			
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);
		}
	}
	

	// A=T^(-1).B
	// T is an upper triangular matrix.
	// Be careful , T is modified at the end of the function
	template <class Field>
	void Field_trsm_up_left (const Field& F,
				  int m, int n,
				  typename Field::Element * B, int ldb,
				  typename Field::Element * T, int ldt,
				  typename Field::Element * A, int lda) {
			LinBox::integer p;
		LinBox::integer max_double("9007199254740991");	F.characteristic(p);
		//LinBox::integer norm= integer(n)*pow(integer(n-1),n>>1)*pow(p,n);		
		//LinBox::integer norm= integer(m)*pow(integer(2),m-1)*pow(p,m);
		LinBox::integer norm=pow(p,m);

		if (norm > max_double )
			{
				int m2= m >>1;
				int m1= m-m2;
				
				Field_trsm_up_left (F,m2,n,B+m1*ldb,ldb,T+m1*ldt+m1,ldt,A+m1*lda,lda);
				Field_dgemm (F,m1,n,m2,-1,T+m1,ldt,A+m1*lda,lda,1,B,ldb);
				Field_trsm_up_left (F,m1,n,B,ldb,T,ldt,A,lda);		       
			}
		else {		
			
			double Td[m*m];
			double Bd[m*n];	

			typedef typename Field::Element Element;	
			Element tmp[m];	 
			// normalisation of T
			for (int i=0;i<m;i++) {
				F.inv(tmp[i],*(T+i*ldt+i));
				for (int j=i;j>=0;j--)
					F.mulin(*(T+i+ldt*j),tmp[i]);
			}	
						
			MatGFq2MatDouble_Triangular (F,m,m,Td,m,T,ldt);
			MatGFq2MatDouble (F,m,n,Bd,n,B,ldb);		
			cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper, CblasNoTrans,CblasUnit,m,n,1.0,Td,m,Bd,n);			
			MatDouble2MatGFq (F,m,n,A,lda,Bd,n);

			// denormalisation of the result A
			for (int i=0;i<m;i++) 
				for (int j=0;j<n;j++)
					F.mulin(*(A+i*lda+j),tmp[i]);
		}
	
	}

	// Main function calling each specific function (depending on tr,un,si)
	// for field_trsm function

	template <class Field>
	inline void Field_trsm (const Field& F,
				int m, int n,
				typename Field::Element * B, int ldb,
				typename Field::Element * T, int ldt,
				typename Field::Element * A, int lda,
				Triangular  tr,
				Unitary     un,
				Side        si) {

		int choice=0;
		if (tr == UPPER)  choice+=1;
		if (un == UNIT)   choice+=2;
		else              choice+=4;
		if (si == LEFT)   choice+=8;
		else              choice+=16;

		switch (choice) {
		case 10:
			Field_trsm_low_left_unit (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		case 11:
			Field_trsm_up_left_unit (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		case 12:
			Field_trsm_low_left (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		case 13:
			Field_trsm_up_left (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		case 18:
			Field_trsm_low_right_unit (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		case 19:
			Field_trsm_up_right_unit (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		case 20:
			Field_trsm_low_right (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		case 21:
			Field_trsm_up_right (F,m,n,B,ldb,T,ldt,A,lda);
			break;
		};
	}


#endif
       		                              

} // end of namespace LinBox

#endif
