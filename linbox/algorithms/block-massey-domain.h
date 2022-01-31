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
#include <givaro/zring.h>
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/factorized-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"
#include "linbox/matrix/polynomial-matrix.h"

#include "linbox/util/timer.h"

//  #define  __CHECK_RESULT
//  #define __DEBUG_MAPLE
//  #define __CHECK_LOOP
//  #define __CHECK_DISCREPANCY
// #define __CHECK_TRANSFORMATION
//  #define __CHECK_SIGMA_RESULT
//#define __PRINT_SEQUENCE
//#define __PRINT_SIGMABASE
//#define __PRINT_MINPOLY

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
        typedef BlasSubmatrix<Coefficient>   CoeffView;


	private:
		Sequence                          *_container;
		const Field                           *_field;
		BlasMatrixDomain<Field>                  _BMD;
		MatrixDomain<Field>                       _MD;
		size_t            EARLY_TERM_THRESHOLD;


	public:

#ifdef _BM_TIMING
		mutable Timer   ttGetMinPoly;			mutable Timer tGetMinPoly;
		mutable Timer	ttNewDiscrepancy;		mutable Timer tNewDiscrepancy;
		mutable Timer	ttShiftSigma;			mutable Timer tShiftSigma;
		mutable Timer   ttApplyPerm;			mutable Timer tApplyPerm;
		mutable Timer   ttUpdateSigma;			mutable Timer tUpdateSigma;
		mutable Timer   ttInverseL;			mutable Timer tInverseL;
		mutable Timer   ttGetPermutation;		mutable Timer tGetPermutation;
		mutable Timer   ttPLUQ;				mutable Timer tPLUQ;
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
			ttPLUQ.clear();
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
			print(ttPLUQ, "PLUQ","direct");
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


		BlockMasseyDomain (const BlockMasseyDomain<Field, Sequence> &Mat, size_t ett_default = DEFAULT_BLOCK_EARLY_TERM_THRESHOLD) :
			_container(Mat._container), _field(Mat._field), _BMD(Mat.field()),
			_MD(Mat.field()),  EARLY_TERM_THRESHOLD (ett_default)
		{
#ifdef _BM_TIMING
			clearTimer();
#endif

		}

		BlockMasseyDomain (Sequence *D, size_t ett_default = DEFAULT_BLOCK_EARLY_TERM_THRESHOLD) :
			_container(D), _field(&(D->field ())), _BMD(D->field ()), _MD(D->field ()), EARLY_TERM_THRESHOLD (ett_default)
		{
#ifdef _BM_TIMING
			clearTimer();
#endif
		}


		// field of the domain
		const Field &field    () const
		{ return *_field; }

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
            const Coefficient Zeromn(field(),m,n);
			std::vector<Coefficient> S (length,Zeromn);


			Coefficient Unit(field(),m+n,m);
			const Coefficient Zero(field(),m+n,m);

			for (size_t i=0;i<m;i++)
				Unit.setEntry(i,i,field().one);
			size_t min_mn=(m <n)? m :n;

			// initialization of discrepancy
			Coefficient Discrepancy(field(),m+n,n);
			for (size_t i=0;i<n;i++)
				Discrepancy.setEntry(i+m,i,field().one);

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

			size_t early_stop=0;
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
				while (_iter_Discr != Discr.End() && (field().isZero(*_iter_Discr)))
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
                  Coefficient ZeroD(field(),m+n,n);
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

				// Computation of the PLUQ decomposition of the discrepancy
				Coefficient CopyDiscr(Discrepancy);
				BlasPermutation<size_t> Pp (CopyDiscr.coldim());
				BlasPermutation<size_t> Qt (CopyDiscr.rowdim());
				PLUQMatrix<Field> PLUQ(CopyDiscr,Pp,Qt);

				// Get the matrix L of PLUQ decomposition
				TriangularBlasMatrix<Field> L(field(),m+n,m+n, Tag::Shape::Lower, Tag::Diag::Unit );
				PLUQ.getL(L);

				// Get the tranposed  permutation of Q from PLUQ
				// BlasPermutation<size_t> Qt=PLUQ.getQ();


				// Computation of permutations BPerm2 such that the last n rows of BPerm2.Qt.Discrepancy are non zero.
				std::vector<size_t> Perm2(m+n);
				for (size_t i=0;i<n;++i)
					Perm2[i]=m+i;
				for (size_t i=n;i<m+n;++i)
					Perm2[i]=i;
				BlasPermutation<size_t> BPerm2(Perm2);

				// compute the inverse of L
				TriangularBlasMatrix<Field> invL (field(),m+n,m+n, Tag::Shape::Lower,Tag::Diag::Unit);
				FFPACK::trinv_left(field(),m+n,L.getPointer(),L.getStride(),invL.getPointer(),invL.getStride());

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
				Givaro::ZRing<long> UF(0);
				// What?
				BlasMatrixDomain<Givaro::ZRing<long> > BMDUF(UF);
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
							field().assign(SigmaBase[i+1].refEntry(m+j,k), SigmaBase[i].getEntry(m+j,k));

						}

				for (size_t j=0;j<n;j++)
					for (size_t k=0;k<n;++k)
						field().assign(SigmaBase[0].refEntry(m+j,k),field().zero);

#ifdef __DEBUG_MAPLE
				report<<"\n\nSigmaBase"<<NN<<":= ";
				write_maple(field(),SigmaBase);

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
					Coefficient Disc(field(),m+n,n);
					for (size_t j=0;j<min_t ;++j)
						_BMD.axpyin(Disc,SigmaBase[j],S[i-j]);
					Disc.write(report)<<std::endl;
				}
#endif

				// Discrepancy= BPerm2.U.Pp from PLUQ
				Coefficient U(field(),m+n,n);
				TriangularBlasMatrix<Field> trU(U,Tag::Shape::Upper,Tag::Diag::NonUnit);
				PLUQ.getU(trU);
				//Discrepancy=U;
				// BlasPermutation<size_t> Pp= PLUQ.getP();
				_BMD.mul(Discrepancy,trU, Pp);
				_BMD.mulin_right(BPerm2,Discrepancy);
			}

            if ( early_stop == EARLY_TERM_THRESHOLD)
				report<<"Early termination is used: stop at "<<NN<<" from "<<length<<" iterations\n\n";

#ifdef __PRINT_SEQUENCE
			report<<"\n\nSequence:= ";
			write_maple(field(),S);
#endif



#ifdef __CHECK_SIGMA_RESULT
			report<<"Check SigmaBase application\n";
			for (size_t i=SigmaBase.size()-1 ;i< length ;++i){
				Coefficient res(field(),m+n,n);
				for (size_t k=0;k<SigmaBase.size();++k)
					_BMD.axpyin(res,SigmaBase[k],S[i-k]);
				res.write(report)<<std::endl;
			}

#endif

			// Get the reverse matrix polynomial of the first m rows of SigmaBase according to degree.
			degree=order;
			long max= *std::max_element(degree.begin(),degree.end());

			Coefficient tmp(field(),m,m);
            P = std::vector<Coefficient> (max+1,tmp);
			for (size_t i=0;i<m;i++)
				for (long j=0;j<=degree[i];j++)
					for (size_t k=0;k<m;k++)
						field().assign(P[degree[i]-j].refEntry(i,k), SigmaBase[j].getEntry(i,k));
#ifdef __CHECK_RESULT
			report<<"Check minimal polynomial application\n";
			bool valid=true;
			for (size_t i=0;i< NN - P.size();++i){
				Coefficient res(field(),m,n);
				_BMD.mul(res,P[0],S[i]);
				for (size_t k=1,j=i+1;k<P.size();++k,++j)
					_BMD.axpyin(res,P[k],S[j]);
				for (size_t j=0;j<m*n;++j)
					if (!field().isZero(*(res.getPointer()+j)))
						valid= false;
			}
			if (valid)
				report<<"minpoly is correct\n";
			else
				report<<"minpoly is wrong\n";
#endif

#ifdef __PRINT_MINPOLY
			report<<"MinPoly:=";
			write_maple(field(),P);
			Coefficient Mat(*_container->getBB());
			report<<"A:=Matrix(";
			Mat.write(report);
#endif

			std::vector<size_t> deg(m);
			for (size_t i=0;i<m;++i)
				deg[i]=(size_t)degree[i];

			return deg;
		}


		std::vector<size_t> masseyblock_left_rec (std::vector<Coefficient> &lingen)
		{
#if (defined  __CHECK_RESULT) or (defined __PRINT_SEQUENCE)
			//std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
            std::ostream& report = std::cout;
#endif
			// Get information of the Sequence (U.A^i.V)
			size_t length = _container->size();
			size_t m, n, mn;
			m = _container->rowdim();
			n = _container->coldim();

			// Set some useful constant
			const Coefficient Zero(field(),2*m,2*m);
            const Coefficient Zeromn(field(),2*m,n);

            mn=m+n;
                        
                        
			// Make the Power Serie from  Sequence (U.A^i.V) and Identity
			//_container->recompute(); // make sure sequence is already computed			
            typedef PolynomialMatrix<Field, PMType::polfirst> PMatrix;
            PMatrix PowerSerie(field(),mn,n,length);
                        
			typename Sequence::const_iterator _iter (_container->begin ());
			for (size_t i=0;i< length; ++i, ++_iter)
				for (size_t j=0;j<m;++j)
					for (size_t k=0;k<n;++k)
						field().assign(PowerSerie.ref(j,k,i), (*_iter).getEntry(j,k));                        
			for (size_t j=0;j<n;++j)
				field().assign(PowerSerie.ref(m+j,j,0),field().one);
#ifdef __PRINT_SEQUENCE
			report<<"PowerSerie:=";
			PowerSerie.write(report);
#endif

            // set the shift to [ 0 .. 0 1 .. 1]
            std::vector<size_t> shift(mn,0);
            std::fill(shift.begin()+m,shift.end(),1);

			// Prepare SigmaBase
			PMatrix SigmaBase(field(),mn,mn,length);

			// Compute OrderBasis up to the order length 
            OrderBasis<Field> SB(field());
            SB.PM_Basis(SigmaBase, PowerSerie, length, shift);


			// take the m rows which have lowest defect
			// compute permutation such that first m rows have lowest defect
			std::vector<size_t> Perm(mn);
			for (size_t i=0;i<mn;++i)
				Perm[i]=i;
			for (size_t i=0;i<mn;++i) {
				size_t idx_min=i;
				for (size_t j=i+1;j<m+n;++j)
					if (shift[j]< shift[idx_min])
						idx_min=j;
				std::swap(shift[i],shift[idx_min]);
				Perm[i]=idx_min;
			}

            // convert to polynomial of matrices
            PolynomialMatrix<Field, PMType::matfirst> Sigma (field(), mn,mn,SigmaBase.size());                        
            Sigma.copy(SigmaBase);

            BlasPermutation<size_t> BPerm(Perm);

			// Apply BPerm to the Sigma Base
			for (size_t i=0;i<Sigma.size();++i){
                auto Sigmai=Sigma[i];
                _BMD.mulin_right(BPerm,Sigmai);
            }

                        
#ifdef __PRINT_SIGMABASE
            report<<"order is "<<length-1<<std::endl;
			   report<<"SigmaBase:=";
            Sigma.write(report);
            report<<"shift:=[";
            std::ostream_iterator<int> out_it (report,", ");
            std::copy ( shift.begin(), shift.end(), out_it );
            report<<"];\n";
                        
#endif

            // Compute the reverse polynomial of Sigma according to row shift 
            size_t max= *std::max_element(shift.begin(),shift.begin()+m);
            //PMatrix lingen(field(),m,m,max+1);
            Coefficient Zeromm(field(),m,m);
            lingen.resize(max+1,Zeromm);
            for (size_t i=0;i<m;i++)
                for (size_t j=0;j<=shift[i];j++)
                    for (size_t k=0;k<m;k++)
                        field().assign(lingen[shift[i]-j].refEntry(i,k), Sigma.ref(i,k,j));
                     
#ifdef __CHECK_RESULT
			report<<"Check minimal polynomial application\n";
			bool valid=true;
			for (size_t i=0;i< length - lingen.size();++i){
				Coefficient res(field(),m,n);
				Coefficient Power(PowerSerie[i],0,0,m,n);
				_BMD.mul(res,lingen[0],Power);
				for (size_t k=1,j=i+1;k<lingen.size();++k,++j){
					Coefficient Powerview(PowerSerie[j],0,0,m,n);
					_BMD.axpyin(res,lingen[k],Powerview);
				}
				for (size_t j=0;j<m*n;++j)
					if (!field().isZero(*(res.getPointer()+j)))
						valid= false;
            }
			if (valid)
				report<<"minpoly is correct\n";
			else
				report<<"minpoly is wrong\n";
#endif

#ifdef __PRINT_MINPOLY
			report<<"MinPoly:=";
            write_maple(field(),lingen);
#endif
			   std::vector<size_t> degree(m);
			   for (size_t i=0;i<m;++i)
			   	degree[i] = shift[i];
			   return degree;
   	}

	}; //end of class BlockMasseyDomain

} // end of namespace LinBox

#endif // __LINBOX_massey_block_domain_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
