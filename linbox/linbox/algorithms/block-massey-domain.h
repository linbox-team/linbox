
/* linbox/algorithms/block-massey-domain.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@lirmm.fr
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */



#ifndef __LINBOX_massey_block_domain_H
#define __LINBOX_massey_block_domain_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"
#include "linbox/field/unparametric.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/matrix/factorized-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/algorithms/sigma-basis.h"


#include "linbox/util/timer.h"

// #define  __CHECK_RESULT
// #define __DEBUG_MAPLE
// #define __CHECK_LOOP
// #define __PRINT_MINPOLY
// #define __CHECK_DISCREPANCY
// #define __CHECK_TRANSFORMATION
// #define __CHECK_SIGMA_RESULT
// #define __PRINT_SEQUENCE
// #define __PRINT_SIGMABASE

//#define _BM_TIMING
#define DEFAULT_BLOCK_EARLY_TERM_THRESHOLD 10

namespace LinBox
{
	template<class Field, class Coefficient>
		void write_maple(const Field& F, const std::vector<Coefficient> & P)
		{
			std::cout<<"Matrix([";
			for (size_t i=0;i< P[0].rowdim();++i){
				std::cout<<"[";
				for (size_t j=0;j< P[0].coldim();++j){
					F.write(std::cout,P[0].getEntry(i,j));
					for (size_t k=1;k<P.size();++k){
						std::cout<<"+ x^"<<k<<"*";
						F.write(std::cout,P[k].getEntry(i,j));
					}
					if (j != P[0].coldim()-1)
						std::cout<<",";
				}
				if (i != P[0].rowdim()-1)
					std::cout<<"],";
				else
					std::cout<<"]";
			}
			std::cout<<"]);\n";
		}





	/** Compute the linear generator of a sequence of matrices.
	 *
	 * This class encapsulates the functionality required for computing
	 * the block minimal polynomial of a matrix.
	 * @bib
	 * Giorgi, Jeannerod, Villard algorithm from ISSAC'03
	 */
	template<class _Field, class _Sequence>
	class BlockMasseyDomain {

	public:
		typedef _Field                           Field;
		typedef typename Field::Element        Element;
		typedef _Sequence                     Sequence;
		typedef BlasMatrix<Field>          Coefficient;
                typedef BlasSubmatrix<Field>         CoeffView;


	private:
		Sequence                          *_container;
		Field                                  _field;
		BlasMatrixDomain<Field>                  _BMD;
		MatrixDomain<Field>                       _MD;
		unsigned long            EARLY_TERM_THRESHOLD;


	public:

#ifdef _BM_TIMING
		mutable Timer   ttGetMinPoly;			mutable Timer tGetMinPoly;
		mutable Timer	ttNewDiscrepancy;		mutable Timer tNewDiscrepancy;
		mutable Timer	ttShiftSigma;			mutable Timer tShiftSigma;
		mutable Timer   ttApplyPerm;			mutable Timer tApplyPerm;
		mutable Timer   ttUpdateSigma;			mutable Timer tUpdateSigma;
		mutable Timer   ttInverseL;			mutable Timer tInverseL;
		mutable Timer   ttGetPermutation;		mutable Timer tGetPermutation;
		mutable Timer   ttLQUP;				mutable Timer tLQUP;
		mutable Timer   ttDiscrepancy;			mutable Timer tDiscrepancy;
		mutable Timer   ttGetCoeff;			mutable Timer tGetCoeff;
		mutable Timer   ttCheckSequence;		mutable Timer tCheckSequence;
		mutable Timer   ttSetup;			mutable Timer tSetup;
		mutable Timer   ttMBasis;			mutable Timer tMBasis;
		mutable Timer   ttUpdateSerie;			mutable Timer tUpdateSerie;
		mutable Timer   ttBasisMultiplication;	        mutable Timer tBasisMultiplication;
		mutable Timer   ttCopyingData;			mutable Timer tCopyingData;
		mutable Timer   Total;

		void clearTimer()
		{
			ttGetMinPoly.clear();
			ttNewDiscrepancy.clear();
			ttShiftSigma.clear();
			ttApplyPerm.clear();
			ttUpdateSigma.clear();
			ttInverseL.clear();
			ttGetPermutation.clear();
			ttLQUP.clear();
			ttDiscrepancy.clear();
			ttGetCoeff.clear();
			ttCheckSequence.clear();
			ttSetup.clear();
			ttMBasis.clear();
			ttUpdateSerie.clear();
			ttBasisMultiplication.clear();
			ttCopyingData.clear(),
			Total.clear();
		}

		void print(Timer& T, const  char* timer, const char* title)
		{
			if (&T != &Total)
				Total+=T;
			if (T.count() > 0) {
				std::cout<<title<<": "<<timer;
				for (int i=(int)strlen(timer); i<28; i++)
					std::cout << ' ';
				std::cout<<T<<std::endl;
			}
		}

		void printTimer()
		{
			print(ttSetup, "Setup", "direct");
			print(ttCheckSequence, "Rank of Seq[0]", "direct");
			print(ttGetCoeff, "Compute sequence", "direct");
			print(ttDiscrepancy, "Compute Discrepancy", "direct");
			print(ttLQUP, "LQUP","direct");
			print(ttGetPermutation, "Compute Permutation", "direct");
			print(ttApplyPerm, "Apply Permutation", "direct");
			print(ttInverseL, "Inverse of L", "direct");
			print(ttUpdateSigma, "Update Sigma", "direct");
			print(ttShiftSigma, "Shift Sigma by x", "direct");
			print(ttNewDiscrepancy, "Keep half Discrepancy", "direct");
			print(ttMBasis, "MBasis computation", "recursive");
			print(ttUpdateSerie, "Updating Power Serie", "recursive");
			print(ttBasisMultiplication, "Basis Multiplication", "recursive");
			print(ttCopyingData, "Copying Data", "recursive");
			print(Total, "Total", "");
			std::cout<<std::endl<<std::endl;
		}
#endif


		BlockMasseyDomain (const BlockMasseyDomain<Field, Sequence> &Mat, unsigned long ett_default = DEFAULT_BLOCK_EARLY_TERM_THRESHOLD) :
			_container(Mat._container), _field(Mat._field), _BMD(Mat._field),
			_MD(Mat._field),  EARLY_TERM_THRESHOLD (ett_default)
		{
#ifdef _BM_TIMING
			clearTimer();
#endif
 
		}

		BlockMasseyDomain (Sequence *D, unsigned long ett_default = DEFAULT_BLOCK_EARLY_TERM_THRESHOLD) :
			_container(D), _field(D->getField ()), _BMD(D->getField ()), _MD(D->getField ()), EARLY_TERM_THRESHOLD (ett_default)
		{
#ifdef _BM_TIMING
			clearTimer();
#endif
		}


		// field of the domain
		const Field &getField    () const
		{ return _field; }

		// sequence of the domain
		Sequence *getSequence () const
		{ return _container; }

		// left minimal generating polynomial of the sequence
		void left_minpoly  (std::vector<Coefficient> &P)
		{
			masseyblock_left(P);
		}

		void left_minpoly_rec  (std::vector<Coefficient> &P)
		{
			masseyblock_left_rec(P);
		}


		// left minimal generating polynomial  of the sequence, keep track on degree
		void left_minpoly (std::vector<Coefficient> &phi, std::vector<size_t> &degree)
		{
			degree = masseyblock_left(phi);
		}

		void left_minpoly_rec  (std::vector<Coefficient> &P, std::vector<size_t> &degree)
		{
			degree = masseyblock_left_rec(P);
		}


		// right minimal generating polynomial of the sequence
		void right_minpoly (std::vector<Coefficient> &P) { masseyblock_right(P);}


	private:

	
 
		std::vector<size_t> masseyblock_left (std::vector<Coefficient> &P)
		{
                        std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

			const size_t length = _container->size ();
			const size_t m = _container->rowdim();
			const size_t n = _container->coldim();

			// ====================================================
			// Sequence and iterator initialization
			// ====================================================

			// Initialization of the sequence iterator
			typename Sequence::const_iterator _iter (_container->begin ());

			// Reservation of memory for the entire sequence
                        const Coefficient Zeromn(_field,m,n);
			std::vector<Coefficient> S (length,Zeromn);


			Coefficient Unit(_field,m+n,m);
			const Coefficient Zero(_field,m+n,m);
                        
			for (size_t i=0;i<m;i++)
				Unit.setEntry(i,i,_field.one);
			size_t min_mn=(m <n)? m :n;

			// initialization of discrepancy
			Coefficient Discrepancy(_field,m+n,n);
			for (size_t i=0;i<n;i++)
				Discrepancy.setEntry(i+m,i,_field.one);

			// initialization of sigma base
			std::vector<Coefficient> SigmaBase(1, Unit);

			// initialization of order of sigma base's rows
			std::vector<long> order(m+n,1);
			for (size_t i=0;i<m;++i)
				order[i]=0;

			// initialisation of degree of sigma base's rows
			std::vector<long> degree(m+n,0);
			for (size_t i=0;i<m;++i)
				degree[i]=0;

			// The first sequence element should be of full rank
			// this is due to the strategy which say that we can compute
			// only the first column of the approximation of [ S(x) Id]^T
			// since the other colums have always lower degree.
			if (_BMD.rank(*_iter)< min_mn)
				throw PreconditionFailed (__func__, __LINE__, "Bad random Blocks, abort\n");

			unsigned long early_stop=0;
			long NN;
			for (NN = 0; (NN < (long)length) && (early_stop < EARLY_TERM_THRESHOLD) ; ++NN, ++_iter) {

				// Get the next coefficient in the sequence
				S[NN]=*_iter;

				/*
				 * Compute the new discrepancy (just updating the first m rows)
				 */
				// view of m first rows of SigmaBasis[0]
				CoeffView Sigma(SigmaBase[0],0,0,m,m);

				// view of m first rows of Discrepancy
				CoeffView Discr(Discrepancy,0,0,m,n);

				_BMD.mul(Discr,Sigma,S[NN]);
                                
				for (size_t i=1;i<SigmaBase.size();i++){
					CoeffView  Sigmaview(SigmaBase[i],0,0,m,m);
					_BMD.axpyin(Discr,Sigmaview,S[NN-i]);
				}

				typename CoeffView::Iterator _iter_Discr = Discr.Begin();
				while (_iter_Discr != Discr.End() && (_field.isZero(*_iter_Discr)))
					++_iter_Discr;
                                if (_iter_Discr!=Discr.End())
                                        early_stop=0;
                                else 
                                        early_stop++;
                               
				// maybe there is something to do here
				// increase the last n rows of orders
				// multiply by X the last n rows of SigmaBase
				//if (_iter_Discr != Discr.End())
				
                                /*
                                Coefficient ZeroD(_field,m+n,n);
                                if (_MD.areEqual(Discrepancy,ZeroD))
                                        early_stop=0;
                                else 
                                        early_stop++;
                                */

				// Computation of the permutation BPerm1 such that BPerm1.order is in increasing order.
				// order=Perm.order
				//! @todo factorize this in \c BlasPermutation.
				std::vector<size_t> Perm1(m+n);
				for (size_t i=0;i<m+n;++i)
					Perm1[i]=i;
				if (NN>=2) {
					for (size_t i=0;i<m+n;++i) {
						size_t idx_min=i;
						for (size_t j=i+1;j<m+n;++j)
							if (order[j]< order[idx_min])
								idx_min=j;
						std::swap(order[i],order[idx_min]);
						Perm1[i]=idx_min;
					}
				}
				BlasPermutation<size_t> BPerm1(Perm1);

				// Discrepancy= BPerm1.Discrepancy
				_BMD.mulin_right(BPerm1,Discrepancy);


#ifdef __CHECK_DISCREPANCY

				report<<"Discrepancy"<<NN<<":=Matrix(";
				Discrepancy.write(report)<<");"<<std::endl;
#endif

				// Computation of the LQUP decomposition of the discrepancy
				Coefficient CopyDiscr(Discrepancy);
				BlasPermutation<size_t> Pp (CopyDiscr.coldim());
				BlasPermutation<size_t> Qt (CopyDiscr.rowdim());
				LQUPMatrix<Field> LQUP(CopyDiscr,Pp,Qt);

				// Get the matrix L of LQUP decomposition
				TriangularBlasMatrix<Field> L(_field,m+n,m+n, LinBoxTag::Lower, LinBoxTag::Unit );
				LQUP.getL(L);

				// Get the tranposed  permutation of Q from LQUP
				// BlasPermutation<size_t> Qt=LQUP.getQ();


				// Computation of permutations BPerm2 such that the last n rows of BPerm2.Qt.Discrepancy are non zero.
				std::vector<size_t> Perm2(m+n);
				for (size_t i=0;i<n;++i)
					Perm2[i]=m+i;
				for (size_t i=n;i<m+n;++i)
					Perm2[i]=i;
				BlasPermutation<size_t> BPerm2(Perm2);

				// compute the inverse of L
				TriangularBlasMatrix<Field> invL (_field,m+n,m+n, LinBoxTag::Lower,LinBoxTag::Unit);
				FFPACK::trinv_left((typename Field::Father_t)_field,m+n,L.getPointer(),L.getStride(),invL.getWritePointer(),invL.getStride());

#ifdef 	__CHECK_TRANSFORMATION
				report<<"invL"<<NN<<":=Matrix(";
				invL.write(report)<<");"<<std::endl;

#endif
				// SigmaBase =  BPerm2.Qt. L^(-1) . BPerm1 . SigmaBase
				for (size_t i=0;i<SigmaBase.size();i++) {
					_BMD.mulin_right(BPerm1,SigmaBase[i]);
					_BMD.mulin_right(invL,SigmaBase[i]);
					_BMD.mulin_right(Qt,SigmaBase[i]);
					_BMD.mulin_right(BPerm2,SigmaBase[i]);
				}

				// Apply BPerm2 and Qt to the vector of order and increase by 1 the last n rows
				UnparametricField<long> UF(0);
				// What?  
				BlasMatrixDomain<UnparametricField<long> > BMDUF(UF);
				BMDUF.mulin_right(Qt,order);
				BMDUF.mulin_right(BPerm2,order);
				BMDUF.mulin_right(BPerm1,degree);
				BMDUF.mulin_right(Qt,degree);
				BMDUF.mulin_right(BPerm2,degree);
				for (size_t i=m;i<m+n;++i){
					order[i]++;
					degree[i]++;
				}

				// Multiplying the last n row of SigmaBase by x.
				long max_degree=degree[m];
				for (size_t i=m+1;i<m+n;++i) {
					if (degree[i]>max_degree)
						max_degree=degree[i];
				}
				size_t size=SigmaBase.size();
				if (SigmaBase.size()<= (size_t)max_degree)
				{
					SigmaBase.resize(size+1,Zero);
					//report << size << std::endl;
					size++;
				}
				//report << "size going in" << size << std::endl;
				for (int i= (int)size-2;i>=0;i--)
					for (size_t j=0;j<n;j++)
						for (size_t k=0;k<n;++k){						
							_field.assign(SigmaBase[i+1].refEntry(m+j,k), SigmaBase[i].getEntry(m+j,k));
                                                        
						}
                                
				for (size_t j=0;j<n;j++)
					for (size_t k=0;k<n;++k)
						_field.assign(SigmaBase[0].refEntry(m+j,k),_field.zero);
                                
#ifdef __DEBUG_MAPLE
				report<<"\n\nSigmaBase"<<NN<<":= ";
				write_maple(_field,SigmaBase);

				report<<"order"<<NN<<":=<";
				for (size_t i=0;i<m+n;++i){
					report<<order[i];
					if (i!=m+n-1) report<<",";
				}
				report<<">;"<<std::endl;
				report<<"degree"<<NN<<":=<";
				for (size_t i=0;i<m+n;++i){
					report<<degree[i];
					if (i!=m+n-1) report<<",";
				}
				report<<">;"<<std::endl;

#endif

#ifdef __CHECK_LOOP
				report<<"\nCheck validity of current SigmaBase\n";
				report<<"SigmaBase size: "<<SigmaBase.size()<<std::endl;
				report<<"Sequence size:  "<<NN+1<<std::endl;
				size_t min_t = (SigmaBase.size() > NN+1)? NN+1: SigmaBase.size();
				for (size_t i=min_t - 1 ; i<NN+1; ++i){
					Coefficient Disc(_field,m+n,n);
					for (size_t j=0;j<min_t ;++j)
						_BMD.axpyin(Disc,SigmaBase[j],S[i-j]);
					Disc.write(report)<<std::endl;
				}
#endif

				// Discrepancy= BPerm2.U.Pp from LQUP
				Coefficient U(_field,m+n,n);
				TriangularBlasMatrix<Field> trU(U,LinBoxTag::Upper,LinBoxTag::NonUnit);
				LQUP.getU(trU);
				//Discrepancy=U;
				// BlasPermutation<size_t> Pp= LQUP.getP();
				_BMD.mul(Discrepancy,trU, Pp);
				_BMD.mulin_right(BPerm2,Discrepancy);
			}

                        if ( early_stop == EARLY_TERM_THRESHOLD)
				report<<"Early termination is used: stop at "<<NN<<" from "<<length<<" iterations\n\n";

#ifdef __PRINT_SEQUENCE
			report<<"\n\nSequence:= ";
			write_maple(_field,S);
#endif



#ifdef __CHECK_SIGMA_RESULT
			report<<"Check SigmaBase application\n";
			for (size_t i=SigmaBase.size()-1 ;i< length ;++i){
				Coefficient res(_field,m+n,n);
				for (size_t k=0;k<SigmaBase.size();++k)
					_BMD.axpyin(res,SigmaBase[k],S[i-k]);
				res.write(report)<<std::endl;
			}

#endif

			// Get the reverse matrix polynomial of the first m rows of SigmaBase according to degree.
			degree=order;
			long max= *std::max_element(degree.begin(),degree.end());

			Coefficient tmp(_field,m,m);
                        P = std::vector<Coefficient> (max+1,tmp);
			for (size_t i=0;i<m;i++)
				for (long j=0;j<=degree[i];j++)
					for (size_t k=0;k<m;k++)
						_field.assign(P[degree[i]-j].refEntry(i,k), SigmaBase[j].getEntry(i,k));
#ifdef __CHECK_RESULT
			report<<"Check minimal polynomial application\n";
			bool valid=true;
			for (size_t i=0;i< NN - P.size();++i){
				Coefficient res(_field,m,n);
				_BMD.mul(res,P[0],S[i]);
				for (size_t k=1,j=i+1;k<P.size();++k,++j)
					_BMD.axpyin(res,P[k],S[j]);
				for (size_t j=0;j<m*n;++j)
					if (!_field.isZero(*(res.getPointer()+j)))
						valid= false;
			}
			if (valid)
				report<<"minpoly is correct\n";
			else
				report<<"minpoly is wrong\n";
#endif

#ifdef __PRINT_MINPOLY
			report<<"MinPoly:=";
			write_maple(_field,P);
			Coefficient Mat(*_container->getBB());
			report<<"A:=Matrix(";
			Mat.write(report);
#endif

			std::vector<size_t> deg(m);
			for (size_t i=0;i<m;++i)
				deg[i]=(size_t)degree[i];
			
			return deg;
		}


		std::vector<size_t> masseyblock_left_rec (std::vector<Coefficient> &P)
		{
#ifdef __CHECK_RESULT
			std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
#endif
			// Get information of the Sequence (U.A^i.V)
			size_t length = _container->size();
			size_t m, n;
			m = _container->rowdim();
			n = _container->coldim();

			// Set some useful constant
			Element one;
			_field.init(one,1UL);
			const Coefficient Zero(_field,2*m,2*m);
                        const Coefficient Zeromn(_field,2*m,n);

			// Make the Power Serie from  Sequence (U.A^i.V) and Identity
			//_container->recompute(); // make sure sequence is already computed
			std::vector<Coefficient> PowerSerie(length,Zeromn);
			typename Sequence::const_iterator _iter (_container->begin ());
			for (size_t i=0;i< length; ++i, ++_iter){
				for (size_t j=0;j<m;++j)
					for (size_t k=0;k<n;++k)
						PowerSerie[i].setEntry(j,k, (*_iter).getEntry(j,k));
			}
			for (size_t j=0;j<n;++j)
				PowerSerie[0].setEntry(m+j, j, one);
#ifdef __PRINT_SEQUENCE
			report<<"PowerSerie:=";
			write_maple(_field,PowerSerie);
#endif


			// Set the defect to [0 ... 0 1 ... 1]^T
			std::vector<size_t> defect(2*m,0);
			for (size_t i=m;i< 2*m;++i)
				defect[i]=1;

			// Prepare SigmaBase
			std::vector<Coefficient> SigmaBase(length,Zero);

			// Compute Sigma Base up to the order length - 1
			SigmaBasis<Field> SB(_field, PowerSerie);
			SB.left_basis(SigmaBase, length-1, defect);

			// take the m rows which have lowest defect
			// compute permutation such that first m rows have lowest defect
			std::vector<size_t> Perm(2*m);
			for (size_t i=0;i<2*m;++i)
				Perm[i]=i;
			for (size_t i=0;i<2*m;++i) {
				size_t idx_min=i;
				for (size_t j=i+1;j<2*m;++j)
					if (defect[j]< defect[idx_min])
						idx_min=j;
				std::swap(defect[i],defect[idx_min]);
				Perm[i]=idx_min;
			}
			BlasPermutation<size_t> BPerm(Perm);

			// Apply BPerm to the Sigma Base
			for (size_t i=0;i<SigmaBase.size();++i)
				_BMD.mulin_right(BPerm,SigmaBase[i]);

#ifdef __PRINT_SIGMABASE
                        report<<"order is "<<length-1<<endl;
			report<<"SigmaBase:=";
			write_maple(_field,SigmaBase);
#endif
			// Compute the reverse polynomial of SigmaBase according to defect of each row
			size_t max=defect[0];
			for (size_t i=0;i<m;++i)
				if (defect[i] > max)
					max=defect[i];

                        const Coefficient tmp(_field,m,m);
			P = std::vector<Coefficient> (max+1,tmp);
			for (size_t i=0;i<m;i++)
				for (size_t j=0;j<=defect[i];j++)
					for (size_t k=0;k<m;k++)
						_field.assign(P[defect[i]-j].refEntry(i,k), SigmaBase[j].getEntry(i,k));

#ifdef __CHECK_RESULT
			report<<"Check minimal polynomial application\n";
			bool valid=true;
			for (size_t i=0;i< length - P.size();++i){
				Coefficient res(_field,m,n);
				Coefficient Power(PowerSerie[i],0,0,m,n);
				_BMD.mul(res,P[0],Power);
				for (size_t k=1,j=i+1;k<P.size();++k,++j){
					Coefficient Powerview(PowerSerie[j],0,0,m,n);
					_BMD.axpyin(res,P[k],Powerview);
				}
				for (size_t j=0;j<m*n;++j)
					if (!_field.isZero(*(res.getPointer()+j)))
						valid= false;
                        }
			if (valid)
				report<<"minpoly is correct\n";
			else
				report<<"minpoly is wrong\n";
#endif

#ifdef __PRINT_MINPOLY
			report<<"MinPoly:=";
			write_maple(_field,P);
#endif
			std::vector<size_t> degree(m);
			for (size_t i=0;i<m;++i)
				degree[i] = defect[i];
			return degree;
		}

	}; //end of class BlockMasseyDomain

} // end of namespace LinBox

#endif // __LINBOX_massey_block_domain_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

