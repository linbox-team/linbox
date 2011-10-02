/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * examples/omp_block_rank.C
 *
 * Copyright (C) 2010 J-G Dumas
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

/*! @file examples/omp_block_rank.C
 * @example  examples/omp_block_rank.C
 * @ingroup examples
 * @brief Block Wiedemann Rank with OpenMP
 */

#include <iostream>
#include <fstream>
#include <omp.h>
#include <givaro/givtimer.h>
#include <givaro/givpoly1crt.h>

#include "linbox-config.h"

// **********************************************************
// Variable globale pour fixer le générateurs des FFT primes
struct FFTSeeder {
	unsigned long seed;
	FFTSeeder(unsigned long s=0) : seed(s) {}
	void setseed(unsigned long s=0) { seed=s; }
	unsigned long operator()() const { return this->seed; }
};
FFTSeeder  FFTgenerator;
#define FFT_PRIME_SEED FFTgenerator()
// **********************************************************


#include "linbox/field/givaro-gfq.h"
#define LINBOX_EXTENSION_DEGREE_MAX 20
#include "linbox/field/givaro-extension.h"
#include "linbox/field/modular-double.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/trace.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/algorithms/sigma-basis.h"
#include "linbox/algorithms/block-massey-domain.h"

template<class Field>
void extractLeftSigma(const Field &F,
		      std::vector<LinBox::BlasMatrix<typename Field::Element> >              &S,
		      std::vector<LinBox::BlasMatrix<typename Field::Element> >      &SigmaBase,
		      std::vector<size_t>                       &defect,
		      size_t                                      block)
{

	typedef typename Field::Element Element;
	LinBox::BlasMatrixDomain<Field> _BMD(F);
	// take the block rows which have lowest defect
	// compute permutation such that first block rows have lowest defect
	std::vector<size_t> Perm(2*block);
	for (size_t i=0;i<2*block;++i)
		Perm[i]=i;
	for (size_t i=0;i<2*block;++i) {
		size_t idx_min=i;
		for (size_t j=i+1;j<2*block;++j)
			if (defect[j]< defect[idx_min])
				idx_min=j;
		std::swap(defect[i],defect[idx_min]);
		Perm[i]=idx_min;
	}
	LinBox::BlasPermutation<size_t>  BPerm(Perm);

	// Apply BPerm to the Sigma Base
	for (size_t i=0;i<SigmaBase.size();++i)
		_BMD.mulin_right(BPerm,SigmaBase[i]);

	size_t max=defect[0];
	for (size_t i=0;i<block;++i)
		if (defect[i] > max)
			max=defect[i];

	// prepare S to receive the sigma base
	const LinBox::BlasMatrix<Element> Zero(block,block);
	S.resize(max+1, Zero);

	// extract the sigma base
	//for (size_t k=0;k<S.size();++k){
	//	for(size_t i=0;i<block;++i)
	//		for (size_t j=0;j<block;++j)
	//			S[k].setEntry(i,j, SigmaBase[k].getEntry(i,j));
	//}

	// extract the reverse sigma base
	for(size_t i=0;i<block;++i)
		for (size_t j=0;j<=defect[i];++j){
			for (size_t k=0;k<block;++k)
				S[defect[i]-j].setEntry(i,k, SigmaBase[j].getEntry(i,k));

		}


}

template<class Field>
void write_sigma(const Field &_F, const char* name, const std::vector<LinBox::BlasMatrix<typename Field::Element> > & P)
{
	size_t m,n;
	m = P[0].rowdim();
	n = P[0].coldim();
	std::cerr<<name<<":=[";
	for (size_t k=0;k<P.size()-1;++k){
		std::cerr<<"Matrix([";
		for (size_t i=0;i<m-1;++i){
			std::cerr<<"[";
			for (size_t j=0;j<n-1;++j)
				_F.write(std::cerr,P[k].getEntry(i,j))<<",";
			_F.write(std::cerr, P[k].getEntry(i,n-1))<<"] , ";
		}
		std::cerr<<"[";
		for (size_t j=0;j<n-1;++j)
			_F.write(std::cerr,P[k].getEntry(m-1,j))<<",";
		_F.write(std::cerr, P[k].getEntry(m-1,n-1))<<"]]) , ";
	}

	std::cerr<<"Matrix([";
	for (size_t i=0;i<m-1;++i){
		std::cerr<<"[";
		for (size_t j=0;j<n-1;++j)
			_F.write(std::cerr,P[P.size()-1].getEntry(i,j))<<",";
		_F.write(std::cerr, P[P.size()-1].getEntry(i,n-1))<<"] , ";
	}
	std::cerr<<"[";
	for (size_t j=0;j<n-1;++j)
		_F.write(std::cerr,P[P.size()-1].getEntry(m-1,j))<<",";
	_F.write(std::cerr, P[P.size()-1].getEntry(m-1,n-1))<<"]])]; \n";
}


template<class Container, class Field>
void scalarmulin(Container& C, const Field& F, const typename Field::Element& x)
{
	for(typename Container::Iterator it=C.Begin();it!=C.End();++it)
		F.mulin(*it, x);
}

template<class Container1, class Field, class Container2>
void contaddin(Container1& C, const Field& F, const Container2& V)
{
	typename Container1::Iterator cit=C.Begin();
	typename Container2::ConstIterator vit=V.Begin();
	for( ; cit!=C.End(); ++cit, ++vit)
		F.addin(*cit, *vit);
}

template<class Field, class Array, class Matrix>
void EvalPolyMat(Array& EvalDets, const Field& F, const LinBox::BlasMatrixDomain<Field>& D, const std::vector<Matrix>& matminpol, const Array& Points)
{
	const long nump = Points.size();
	std::cerr << "num procs: " << omp_get_num_procs() << std::endl;
	std::cerr << "max threads: " << omp_get_max_threads() << std::endl;
	std::cerr << "eval points: " << nump << std::endl;
#pragma omp parallel for schedule(static)
	for(int i=0; i<nump; ++i) {
		const long degree = matminpol.size()-1;
		Matrix mat=matminpol[degree];
		for(int j=degree-1;j>=0;--j) {
			scalarmulin(mat, F, Points[i]);
			contaddin(mat, F, matminpol[j]);
		}
		EvalDets[i] = D.det(mat);
	}

	//     for(int i=0; i<nump; ++i) {
	//         std::cerr << Points[i] << "\t: " << EvalDets[i] << std::endl;
	//     }
}

#include "linbox/algorithms/whisart_trace.h"


#if 0 // now in LinBox
template<class Field, class BB>
void WhisartTrace(
		  typename Field::Element& trace,
		  const Field& F,
		  const LinBox::Diagonal<Field>& ExtD,
		  const BB& A,
		  const LinBox::Diagonal<Field>& InD) {
	// Trace of ExtD B InD B^T ExtD
	// is sum ExtD_i^2 B_{i,j} InD_j
	F.init(trace, 0);
	for(typename BB::ConstIndexedIterator it = A.IndexedBegin();
	    it != A.IndexedEnd(); ++it) {
		typename Field::Element tmp,e,i; F.init(tmp);F.init(e);F.init(i);
		F.mul(tmp,it.value(),it.value());
		ExtD.getEntry(e, it.rowIndex(),it.rowIndex());
		InD.getEntry(i, it.colIndex(),it.colIndex());
		F.mulin(tmp,e);
		F.mulin(tmp,e);
		F.mulin(tmp,i);
		F.addin(trace, tmp);
	}
}

template<class Field, class BB>
void WhisartTraceTranspose(
			   typename Field::Element& trace,
			   const Field& F,
			   const LinBox::Diagonal<Field>& ExtD,
			   const BB& A,
			   const LinBox::Diagonal<Field>& InD) {
	// Trace of ExtD B^T  InD B ExtD
	// is sum ExtD_j^2 B_{i,j} InD_i
	F.init(trace, 0);
	for(typename BB::ConstIndexedIterator it = A.IndexedBegin();
	    it != A.IndexedEnd(); ++it) {
		typename Field::Element tmp,e,i; F.init(tmp);F.init(e);F.init(i);
		F.mul(tmp,it.value(),it.value());
		ExtD.getEntry(e, it.colIndex(),it.colIndex());
		InD.getEntry(i, it.rowIndex(),it.rowIndex());
		F.mulin(tmp,e);
		F.mulin(tmp,e);
		F.mulin(tmp,i);
		F.addin(trace, tmp);
	}
}
#endif

using namespace Givaro;


template<class Field>
int OMP_BLOCK_RANK_main (const Field& F, int argc, char **argv)
{
	LinBox::commentator.setMaxDetailLevel (-1);
	LinBox::commentator.setMaxDepth (-1);
	LinBox::commentator.setReportStream (std::cerr);

	OMPTimer chrono1,chrono2,chrono3,chrono4; chrono1.clear(); chrono2.clear(); chrono3.clear(); chrono4.clear();

	Integer c; F.cardinality(c);
	unsigned long seed = (argc>4?atoi(argv[4]):0);
	FFTgenerator.setseed(seed);

	typename Field::RandIter generator (F,c,seed);
	LinBox::VectorDomain<Field> VD(F);

	std::ifstream input (argv[1]);
	LinBox::MatrixStream<Field> ms( F, input );
	typedef LinBox::SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq > Blackbox;
	typedef LinBox::BlasBlackbox<Field> Block_t;

	Blackbox B (ms);

	std::cerr << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
	long M = B.rowdim() ;
	long N = B.coldim();
	long R = (M<N?M:N);
	long S = (M>N?M:N);

	int nb = (argc>3 ? atoi(argv[3]) : omp_get_max_threads() );
	std::cerr << "block size: " << nb << std::endl;

	chrono1.start();
	std::vector< std::vector< typename Field::Element > > LV(nb);
	for (typename std::vector< std::vector< typename Field::Element > >::iterator it = LV.begin(); it != LV.end(); ++it) {
		it->resize(R);
		for(typename std::vector<typename Field::Element>::iterator vit=it->begin(); vit != it->end(); ++vit)
			generator.random( *vit );
		// 		VD.write(output, *it) << std::endl;
	}

	Block_t LM(F, nb, R);
	for (size_t i=0; i < LM.rowdim(); ++i) {
		for (size_t j=0; j < LM.coldim(); ++j) {
			LM.setEntry(i,j, LV[i][j]);
		}
	}

	chrono1.stop();

	std::cerr << "Generated " << nb << ' ' << R << "-vectors" << std::endl;
	std::cerr << chrono1 << std::endl;

	chrono2.start();

	typedef LinBox::BlasMatrix<typename Field::Element>        Matrix;
	typedef std::vector<Matrix>   Polynomial;


	const Matrix SigmaZero(2*nb,2*nb);
	const Matrix SerieZero(2*nb,nb);

	long d = 4+(R<<1)/nb;

	// define the serie and the sigmabase
	Polynomial Serie(d, SerieZero);
	Polynomial Sigma(d, SigmaZero);


	std::vector<typename Field::Element> d1(S);
	for (int i = 0; i < S; i++)
		do generator.random (d1[i]); while (F.isZero (d1[i]));
	LinBox::Diagonal<Field> D1 (F, d1);
	std::vector<typename Field::Element> d2(R);
	for (int i = 0; i < R; i++)
		do generator.random (d2[i]); while (F.isZero (d2[i]));
	LinBox::Diagonal<Field> D2 (F, d2);

	std::cerr << "num procs: " << omp_get_num_procs() << std::endl;
	std::cerr << "max threads: " << omp_get_max_threads() << std::endl;

	if (M>N) {
#pragma omp parallel for firstprivate(B) schedule(static)
		for(int j=0; j<nb; ++j) {
			std::vector< typename Field::Element > v(S),u(S),w(R);
			std::vector< typename Field::Element > colonne(nb);

			LM.apply(colonne, LV[j]);
			for(int i=0;i<nb;++i)
				Serie[0].setEntry(i,j,colonne[i]);

			for(int k=1;k<d;++k) {
				// BlackBox Apply
				D2.apply(w,LV[j]);
				B.apply(v,w);
				D1.apply(u,v);
				B.applyTranspose(w,u);
				D2.apply(LV[j],w);

				// Dot products
				LM.apply(colonne, LV[j]);
				for(int i=0;i<nb;++i)
					Serie[k].setEntry(i,j,colonne[i]);
			}
			std::cerr << "Thread[" << omp_get_thread_num() << "]: Done BTB-Serie[k][" << j << ']' << std::endl;
		}
	}
	else {
#pragma omp parallel for firstprivate(B) schedule(static)
		for(int j=0; j<nb; ++j) {
			std::vector< typename Field::Element > v(S),u(S),w(R);
			std::vector< typename Field::Element > colonne(nb);

			LM.apply(colonne, LV[j]);
			for(int i=0;i<nb;++i)
				Serie[0].setEntry(i,j,colonne[i]);

			for(int k=1;k<d;++k) {
				// BlackBox Apply
				D2.apply(w,LV[j]);
				B.applyTranspose(v,w);
				D1.apply(u,v);
				B.apply(w,u);
				D2.apply(LV[j],w);

				// Dot products
				LM.apply(colonne, LV[j]);
				for(int i=0;i<nb;++i)
					Serie[k].setEntry(i,j,colonne[i]);

			}
			std::cerr << "Thread[" << omp_get_thread_num() << "]: Done BBT-Serie[k][" << j << ']' << std::endl;

		}
	}

	chrono2.stop();
	std::cerr << "Computed a degree " << d << ' ' << nb << 'x' << nb << "-series" << std::endl;
	std::cerr << chrono2 << std::endl;
	// write_sigma(F, "serie", Serie);

	chrono3.start();
	// append Identity to the serie
	for (int i=0;i<nb;++i)
		F.init(Serie[0].refEntry(nb+i,i), 1);

	// define defect
	std::vector<size_t> defect(2*nb,0);
	for (int i=nb;i<2*nb;++i){
		defect[i]=1;
	}


	LinBox::SigmaBasis<Field> SB(F, Serie);
	std::cerr<<"blockminpoly computation... ";
	SB.PM_Basis(Sigma, Serie, d-1, defect);
	std::cerr<<"done\n";
	// write_sigma(F, "serie", Serie);
	// write_sigma(F, "sigma", Sigma);

	std::cerr<<"extracting bminpoly... ";
	std::vector<Matrix> LS2;
	extractLeftSigma(F, LS2, Sigma, defect, nb);
	std::cerr<<"done\n";
	std::cerr<<"Rank of the highest degree coefficient...";
	unsigned long rdeg;
	LinBox::BlasMatrixDomain<Field> D(F);
	rdeg = D.rank(LS2[LS2.size()-1]);
	typename Field::Element d0,de;
	d0 = D.det(LS2[0]);
	de = D.det(LS2[LS2.size()-2]);

	if( ! ( F.isZero(d0) || F.isZero(de) ) )
	{
		int rank =  ((LS2.size()-2)*(LS2[0].rowdim())+rdeg);

		chrono3.stop();
		std::cerr<<"is " << rdeg << ", done.\n";
		std::cerr << "Estimated rank: " << rank << std::endl;
	}
	else {
		chrono3.stop();
		std::cerr<<  "\n*** WARNING *** Insufficient information, interpolation required. You might also try again with a larger field.\n";
		F.write(std::cerr<<  "det(bm[0]): ", d0) << std::endl;
		std::cerr<<  "rk(bm[0]): " << D.rank(LS2[0])<< std::endl;
		F.write(std::cerr<<  "det(bm[" << (LS2.size()-2) << "]): ", de) << std::endl;
		std::cerr<<  "rk(bm[" << (LS2.size()-2) << "]): "<< D.rank(LS2[LS2.size()-2]) << std::endl;
		std::cerr<<  "rk(bm[" << (LS2.size()-1) << "]): "<< rdeg << std::endl;
		long def=0;
		if (F.isZero(d0)) ++def;
		if (F.isZero(de)) ++def;
		int rank =  ((LS2.size()-2-def)*(LS2[0].rowdim())+rdeg);
		std::cerr<< "*** VERY ROUGH *** rank approximation  " << rank << std::endl;
	}
	std::cerr << "recursive PMBasis CPU time (s)  : " << chrono3.usertime() << std::endl<<std::endl;
	std::cerr << chrono3 << std::endl;

	chrono4.start();


	std::cerr<<"Interpolation of matrix minpoly determinant ...";

	// write_sigma(F, "bminpoly", LS2);

	typedef Poly1CRT< LinBox::GivaroField<Field> >  PolyCRT;
	typedef typename PolyCRT::array_T VScal;

	VScal EvalPoints( LS2.size()*nb );
	for(typename VScal::iterator itp = EvalPoints.begin(); itp != EvalPoints.end(); ++itp) {
		do {
			do generator.random (*itp); while (F.isZero (*itp));
		} while ( (std::find(EvalPoints.begin(), itp, *itp) != itp)) ;
	}

	VScal EvalDets( EvalPoints.size() );
	EvalPolyMat(EvalDets, F, D, LS2, EvalPoints);


	PolyCRT Interpolator(F, EvalPoints, "Y");
	typename PolyCRT::Element Determinant;

	Interpolator.RnsToRing(Determinant, EvalDets);
	chrono4.stop();

	std::cerr << "done\n";


	//     Interpolator.write(std::cerr << "Determinant of MinPoly: ", Determinant) << std::endl;

	Degree deg; Interpolator.getpolydom().degree(deg, Determinant);
	Degree val; Interpolator.getpolydom().val(val, Determinant);


	F.write(std::cerr, Determinant[0]) << " + ";
	if (val > 0) F.write(std::cerr<< "... + ", Determinant[val.value()]) << "Y^" << val << " + ";
	std::cerr << "... + ";
	if (Determinant.size() >= 2) {
		F.write(std::cerr, Determinant[Determinant.size()-2]) << "Y^" << (deg-1) << " + ";
	}
	F.write(std::cerr, Determinant[Determinant.size()-1]) << "Y^" << deg << std::endl;


	typename Field::Element t, p2; F.init(p2, 0UL);
	if (Determinant.size() >= 2) {
		F.neg(p2, Determinant[ Determinant.size()-2]);
		F.divin(p2, Determinant[ Determinant.size()-1]);
	}

	if (M>N) {
		typedef LinBox::Transpose<Blackbox> AT_BB;
		typedef LinBox::Compose< LinBox::Diagonal<Field>, AT_BB > DAT_BB ;
		typedef LinBox::Compose< DAT_BB, LinBox::Diagonal<Field> > DATD_BB;
		typedef LinBox::Compose< DATD_BB, Blackbox> DATDA_BB;
		typedef LinBox::Compose< DATDA_BB, LinBox::Diagonal<Field> > DATDAD_BB;

		AT_BB 		AT(&B);
		DAT_BB 		B1(&D2, &AT);	// B1 = D2 B^T
		DATD_BB 	B2(&B1, &D1);	// B2 = B1 D1 = D2 B^T D1
		DATDA_BB 	B3(&B2, &B);	// B3 = B2 B  = D2 B^T D1 B
		DATDAD_BB	B4(&B3, &D2);	// B4 = B3 D2 = D2 B^T D1 B D2
		std::cerr << "B4: " << B4.rowdim() << "x" << B4.coldim() << std::endl;

		//         trace(t, B4);
		LinBox::WhisartTraceTranspose(t, F, D2, B, D1);

		F.write(std::cerr << "Trace D2 B^T D1 B D2: ", t) << std::endl;

	}
	else {
		typedef LinBox::Compose< LinBox::Diagonal<Field>, Blackbox > DA_BB;
		typedef LinBox::Compose< DA_BB, LinBox::Diagonal<Field> > DAD_BB;
		typedef LinBox::Transpose<Blackbox> AT_BB;
		typedef LinBox::Compose< DAD_BB, AT_BB> DADAT_BB;
		typedef LinBox::Compose< DADAT_BB, LinBox::Diagonal<Field> > DADATD_BB;

		DA_BB 		B1(&D2, &B);	// B1 = D2 B
		DAD_BB 		B2(&B1, &D1);	// B2 = B1 D1 = D2 B D1
		AT_BB 		AT(&B);
		DADAT_BB 	B3(&B2, &AT);	// B3 = B2 B^T= D2 B D1 B^T
		DADATD_BB	B4(&B3, &D2);	// B4 = B3 D2 = D2 B D1 B^T D2
		std::cerr << "B4: " << B4.rowdim() << "x" << B4.coldim() << std::endl;

		//         trace(t, B4);
		LinBox::WhisartTrace(t, F, D2, B, D1);

		F.write(std::cerr << "Trace D2 B D1 B^T D2: ", t) << std::endl;
	}

	if (! F.areEqual( t, p2 )) {
		std::cerr << "*** FAILURE (" << (deg-val) << ") ***" << std::endl;
		F.write(std::cerr << "Trace: ", t) << std::endl;
		F.write(std::cerr << "Minpo: ", p2) << std::endl;
	}
	else {
		//         std::cerr << "Degree - valuation: " << (deg-val) << std::endl;
		std::cerr << "MONTE CARLO RANK: " << (deg-val) << std::endl;
	}
	std::cerr << chrono4 << std::endl;



	std::cerr << "Rank |\t Time-genvect\t Time-seq\t Time-SB\t Time-Interp |\t Total-Time" << std::endl;
	std::cout << (deg-val) << std::scientific << std::setprecision(3)
	<< " |\t" << chrono1.usertime()
	<< '\t' << chrono2.usertime()
	<< '\t' << chrono3.usertime()
	<< '\t' << chrono4.usertime()
	<< " |\t" << (chrono1.usertime()+chrono2.usertime()+chrono3.usertime()+chrono4.usertime())
	<< std::endl;
	return 0;
}

int main (int argc, char **argv)
{
	int c = (argc>2 ? atoi(argv[2]) : 65521);
	unsigned long extend = (unsigned long)FF_EXPONENT_MAX(c,(int)LINBOX_EXTENSION_DEGREE_MAX);
	if (extend > 1) {
		std::cerr << "*** WARNING *** would be best using an extension field of degree " << extend << std::endl;
	}
#if 0
	if (extend > 1) {
		typedef LinBox::GivaroGfq Field;
		Field EF( (unsigned long)c, extend);
		EF.write(std::cerr << "Using an extension field ") << std::endl;
		return OMP_BLOCK_RANK_main(EF,argc,argv);
	}
#endif
	// else {
	typedef LinBox::Modular<double> Field;
	Field F(c);
	return OMP_BLOCK_RANK_main(F,argc,argv);
	//     }

}

