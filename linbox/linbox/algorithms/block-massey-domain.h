/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/block-massey-domain.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
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


/* 
 * Be careful with the masseydomain function. It's a function template by a matrix type which needs specific function.
 * The inerface of a such matrix is given by the class DenseMatrix provided in the file "dense.h"
 *
 */



#ifndef __MASSEY_BLOCK_DOMAIN_H
#define __MASSEY_BLOCK_DOMAIN_H

#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include "linbox/blackbox/dense.h"
#include "linbox/matrix/matrix-domain.h"


namespace LinBox 
{


#define DEFAULT_EARLY_TERM_THRESHOLD 20
#define DEFAULT_ADDITIONAL_ITERATION 2



	template<class Field, class Sequence>
	class MasseyBlockDomain {
	private:
		Sequence            *_container;
		Field                _F;
		VectorDomain<Field>  _VD;
		unsigned long         EARLY_TERM_THRESHOLD;
		MatrixDomain<Field>  _MD;

	public:
		typedef typename Field::Element Element;

		//-- Constructors
		MasseyBlockDomain (unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
			: _container           (), 
			  _F                   (), 
			  _VD                  (_F),
			  EARLY_TERM_THRESHOLD (ett_default)
		{}

		MasseyBlockDomain (const MasseyBlockDomain<Field, Sequence> &M, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
			: _container           (M._container), 
			  _F                   (M._F), 
			  _VD                  (M._F),
			  EARLY_TERM_THRESHOLD (ett_default)
		{}

		MasseyBlockDomain (Sequence *D, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
			: _container           (D), 
			  _F                   (D->getField ()),
		  _VD                  (D->getField ()),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}
  
	MasseyBlockDomain (Sequence *MD, const Field &F, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (MD), 
		  _F                   (F),
		  _VD                  (F),
		  EARLY_TERM_THRESHOLD (ett_default) 
		{}

        //-- Principal method
	template<class Matrix>
	void operator () (vector<Matrix> &C, bool full_poly = false) {
		masseyblock (C, full_poly);
	};

        //-- Domains access
	const Field &getField    () const { return _F; }
	Sequence    *getSequence () const { return _container; }

    private:
	// -----------------------------------------------
	// Polynomial emulation
	// Only container aspects of polynomials
	// AND degree and valuation are needed !
	// -----------------------------------------------


	// -------------------------------------------------------------------
	// Berlekamp/Massey algorithm by block with Massey's Sequence generation
	// -------------------------------------------------------------------

	template<class Matrix>
	  vector<long>& non_zero_row(vector<long>& Res, Matrix& M,const vector<long>& V,
				      long j) const
		{	    
			long nbr_elem=0;
			long degree = 10000000 ;
			Element tmp;
			_F.init(tmp,0UL);

			for (long i=0;i<M.rowdim();i++)
				{		
					M.getEntry(tmp,i,j);
					if (!_F.isZero(tmp))
						{
							Res.push_back(i);
							nbr_elem++;			    
					    
							if  (V[i] < degree)
								{			      
									degree = V[i];
									swap(Res[nbr_elem-1],Res[0]);
								}	     		
						}
				}	    	    

	   
			return Res;
		}
	    
		template<class Matrix>
		ostream& PrintMapleMatrix(Matrix& M ,ostream& os) {
			
			long current_row = 1;
			long current_col = 1;
			
			typename Matrix::RowIterator p;
			os<<"Matrix( [";
			for (p = M.rowBegin (); p != M.rowEnd (); ++p) {
				typename Matrix::ConstRow::const_iterator pe;
				os<<"[ ";
				for (pe = p->begin (); pe != p->end (); ++pe) {
					if (current_col != M.coldim()) {
						_F.write (os, *pe)<<",";
						current_col++;
					}
					else {
						_F.write (os, *pe)<<"]";	
						current_col=1;
					}
				}
				if (current_row != M.rowdim()) {
					os<<",";
					current_row++;
				}
				else 
					os<<"] )";
			}
			return os;
		}
		

		template<class Matrix>
		void PrintMapleMatrixPolynomial(vector<Matrix>& P,ostream& out)
		{
			for (long i=0;i<P.size()-1;i++)
				PrintMapleMatrix(P[i],out)<<" , ";
			PrintMapleMatrix(P[P.size()-1],out);
		}
		

		

        template<class Matrix>
	long masseyblock (vector<Matrix> &C, bool full_poly = false) { 


		typedef typename Matrix::RowIterator RowIterator;
		typedef typename Matrix::ConstRowIterator ConstRowIterator;

		const long END = _container->size () + (full_poly ? DEFAULT_ADDITIONAL_ITERATION:0) ;
		const long n = END >> 1;
		
		long p,q;
		p = _container->row();
		q = _container->col();
		
		commentator.start ("Massey", "LinBox::MasseyDomain::massey", END);

		// ====================================================
		// Sequence and iterator initialization
		// ====================================================
		
		// Initialization of the iterator
		typename Sequence::const_iterator _iter (_container->begin ());

		// Reservation of memory for the entire sequence
		vector<Matrix> S (END + 1,Matrix(_F,p,q));
		
		// Definition of useful data
		Element one;
		_F.init(one,1UL);
		const Matrix ZeroB(_F,q,p);
	        Matrix ZeroC(_F,p,p);
		Matrix One(_F,q,p);
		for(long i=0;i<min(q,p);i++)
			One.setEntry(i,i,one);
		

		// Reservation of memory for  matrix polynomial 		
		C = vector<Matrix> (n+1,ZeroC); C.resize(1,ZeroC);

		// Initialization of the current minimal matrix polynomial to the constant polynomial Identity.
		for(long i=0;i<p;i++)
			C[0].setEntry(i,i,one);

		// Definition and memory reservation of the last non-modified minimal matrix polynomial
		vector<Matrix> B (n+1,ZeroB); B.resize (1,ZeroB); 
		B[0] = One;
		
		// Definition of the discrepancy matrices
		Matrix d(_F,p,q);
		Matrix b(_F,q,q);

		// initialization of the last non zero discrepancy
		for (long i=0;i<q;i++)
			b.setEntry(i,i,one);

		// Useful Data to make the determinant of the minimal matrix polynomial as a monic polynomial
		Element monic;
		_F.init(monic,1UL);
		
		// tips to not compute at each iteration the multiplication by x of the pivot row.
		long x[q];
		for (long i=0;i<q;i++)
			x[i]=1;
		
		// Initialization of the row degree of B and C .
		vector<long> B_deg(q,0),C_deg(p,0);	

		
		for (long N = 0; N < END && x[0] < EARLY_TERM_THRESHOLD; ++N, ++_iter) {
			if (!(N % 1000)) 
				commentator.progress (N);
		  
			// Get the next coefficient in the sequence		 
			S[N] = *_iter;
		  

			// Discrepancy Computation
			_MD.mul(d,C[0],S[N]);
			Matrix tmp(_F,p,q);
			for (long i = 1;i< C.size();i++)
				{
					_MD.mul(tmp,C[i],S[N - i]);
					_MD.addin(d,tmp);
				}
		  
			// the discrepancy is zero so multiply B by x.
			if (_MD.isZero(d)) {
				for (long i=0;i<q;i++) x[i]++;
			}
			// the discrepancy is not zero
 			else {
				for (long j=0;j<q;j++) {
					
					//Computation of the non zero rows index of C
					vector<long> _row_index;
					non_zero_row(_row_index,d,C_deg,j);
				
					if ( _row_index.size()) {
						if ( C_deg[_row_index[0]] > B_deg[j]+x[j])
							{ //case where the pivot is in B .
								/*
								 * Elimination of each row of C by the pivot and multiply the pivot row by x
								 */
								for (long ii=0;ii<_row_index.size();ii++) {
									// Computation of - d[i,j] / b[j,j]
									Element tmp,elt;
									_F.init(tmp);
									_F.init(elt);
									_F.neg(tmp,d.getEntry(elt,_row_index[ii],j));	
									_F.divin(tmp,b.getEntry(elt,j,j));
							  
									// Computation of C[i] = C[i] - (d[i,j] / b[j,j]) B[j] 
									for (long k=x[j];k<B_deg[j]+1+x[j];k++) {
										RowIterator lineC=C[k].rowBegin ()+_row_index[ii];
										RowIterator lineB=B[k-x[j]].rowBegin ()+j;
										_VD.axpyin(*lineC,tmp, *lineB);
									}
							  
									// Computation of d[i] = d[i] - (d[i,j] / b[j,j]) b[j]
									RowIterator lined = d.rowBegin ()+_row_index[ii];
									RowIterator lineb = b.rowBegin ()+j;
									_VD.axpyin(*lined,tmp,*lineb);			
								}
								x[j]++;

							}
						else { // case where the pivot is in C .
							/*
							 * Elimination of each row in C different from the pivot.
							 */
							for (long ii=1;ii<_row_index.size();ii++) {
								
								
								// computation of - d[i,j]/ d[i0,j]  .
								Element tmp,elt;
								_F.init(tmp);
								_F.init(elt);
								_F.neg(tmp,d.getEntry(elt,_row_index[ii],j));
								_F.divin(tmp,d.getEntry(elt,_row_index[0],j));
						  
								// computation of C[ii] = C[ii] - (d[ii,j]/ d[i0,j]) C[i0]
								for (long k=0;k<C_deg[_row_index[0]]+1;k++) {
									RowIterator lineC1= C[k].rowBegin()+_row_index[ii];
									RowIterator lineC2= C[k].rowBegin()+_row_index[0];					       
									_VD.axpyin(*lineC1,tmp,*lineC2 ); 
								}					      
								// computation of d[ii] = d[ii] - (d[ii,j]/ d[i0,j]) d[i0]
								RowIterator lined1 = d.rowBegin ()+_row_index[ii];
								RowIterator lined2 = d.rowBegin ()+_row_index[0];
								_VD.axpyin(*lined1,tmp,*lined2);
							}		
							
							/* 
							 * Elimination of the j'th row of B with the pivot 
							 */
							// computation of - b[j,j]/d[i0,j]
							Element tmp2,elt2;
							_F.init(tmp2);
							_F.init(elt2);
							_F.neg(tmp2,b.getEntry(elt2,j,j));
							_F.divin(tmp2,d.getEntry(elt2,_row_index[0],j));
			      
							_F.mulin(monic,tmp2);
			       			       
							
			       
							// resizing B and C due to the elimination above
							if (B.size() < B_deg[j]+x[j]+1) {
								long Bsize = B.size();
								B.resize(B_deg[j]+x[j]+1,ZeroB);
								for (long k=Bsize;k<B_deg[j]+x[j]+1;k++)
									B[k] = ZeroB;
							}
			       
							if (C.size() < B_deg[j]+x[j]+1) {				    
								long Csize = C.size();
								C.resize(B_deg[j]+x[j]+1,ZeroC);
								for (long k=Csize;k<B_deg[j]+x[j]+1;k++)
									C[k] = ZeroC;
							}
			       
							// computation of B[j] = B[j] - (b[j,j]/d[i0,j]) C[i0] 
							// and swaping B[j] and C[i0]
							for (long k= B_deg[j]+x[j];k>x[j]-1;k--) {
								RowIterator lineB1= B[k].rowBegin()+j;
								RowIterator lineB2= B[k-x[j]].rowBegin()+j;
								RowIterator lineC= C[k].rowBegin()+_row_index[0];
								_VD.axpy(*lineB1,tmp2,*lineC,*lineB2);
								swap_ranges ((*lineB1).begin(),(*lineB1).end(),(*lineC).begin());
							}
			       
							for (long k=x[j]-1;k>=0;k--) {
								RowIterator lineB= B[k].rowBegin()+j;
								RowIterator lineC= C[k].rowBegin()+_row_index[0];
								_VD.mul(*lineB,*lineC,tmp2);
								//swap(lineB,lineC);
								swap_ranges ((*lineB).begin(),(*lineB).end(),(*lineC).begin());
							}
			       
			       
							// computation of b[j] = b[j]  - (b[j,j]/d[i0,j]) d[i0]
							RowIterator lineb= b.rowBegin()+j;
							RowIterator lined= d.rowBegin()+_row_index[0];
							_VD.axpyin(*lineb,tmp2,*lined);
			       
							// exchange b[j] and d[i0]			      
							//swap (lineb,lined);
							swap_ranges ((*lineb).begin(),(*lineb).end(),(*lined).begin());
			       
							// exchange C_deg[i0] and B_deg[j]
							long toto;
							toto = C_deg[_row_index[0]];
							C_deg[_row_index[0]] = B_deg[j]+x[j];
							B_deg[j] = toto;
							// increasing the pivot row by x.
							x[j]=1;
			      			      			       

						}
					}
					else { x[j]++;} 
				}
			}
		}				    				    				  				     
				
		// Get the reversal polynomial of C according to the degree of each rows and make his determinant monic.
		
		for (long i=0;i<p;i++) {
			for (long j=0;j<(C_deg[i]>>1)+1;j++) {
				RowIterator lineC1= C[j].rowBegin()+i;
				RowIterator lineC2= C[C_deg[i]-j].rowBegin()+i;
				//swap(lineC1,lineC2);
				swap_ranges ((*lineC1).begin(),(*lineC1).end(),(*lineC2).begin());
			}
			
		}
		
		_F.invin(monic); 
		
		for (long k=0;k< C_deg[0]+1;k++) {
			RowIterator lineC1= C[k].rowBegin();
			RowIterator lineC2= C[k].rowBegin();
			_VD.mul(*lineC1,*lineC2,monic);
		}
		
				
		commentator.stop ("Done", "Done", "LinBox::MasseyDomain::massey");
				

		PrintMapleMatrixPolynomial(C ,cout);
			

		return C.size();
	}
		


		
	

public:
	// ---------------------------------------------
	// Massey Block
	//

		// this function compute the minimal matrix polynomial as a vector of matrices 
		template<class Matrix>
		void minpoly_block (std::vector<Matrix> &phi, bool full_poly = true) {
			masseyblock<Matrix>(phi,1);
		}


		// this function is just a test of the masseyblock function
		long test_massey(long& dp)
		{
			vector<DenseMatrix<Field> > phi;
			return dp = masseyblock<DenseMatrix<Field> >(phi,1);
		}
		
	};
	
} // end of namespace LinBox

#endif // __MASSEY_DOMAIN_H
