
/* tests/test-bmseq.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by George Yuhasz <yuhasz@gmail.com>
 *
 * --------------------------------------------------------
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <set>
#include <list>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/algorithms/bm-seq.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/matrix/sparse.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;
using namespace std;


int main (int argc, char **argv)
{
	bool pass = true;

	static integer n = 5;
	static integer q = 101;
	//static integer q = 5;
	static int iterations = 2; // was 100
	//static int N = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",   TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator().start("bmseq test suite", "BlasMatrix");

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);



	Field::Element one, zero;
	F.init(one,1);
	F.init(zero,0);
	MatrixDomain<Field> MD(F);
#if 0
	BlasMatrix<Field> D(F,2,2);
	BlasMatrix<Field> zero24(F,2,4);
	for(size_t i=0; i<2; i++)
		D.setEntry(i,i,one);
	D.setEntry(1,0,one);
	BM_Seq<Field> seq(2,D);
	BlasMatrix<Field> S2(F,2,2);
	BM_Seq<Field>::BM_iterator bmit(seq, 0), bmit2(seq.BM_begin());
	bmit.setDelta(4);
	BM_Seq<Field>::BM_iterator::TerminationState check = bmit.state();
	while(!check.IsGeneratorFound() && !check.IsSequenceExceeded()){
					bmit++;
					check = bmit.state();
	}
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	bmit++;
	check = bmit.state();
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	MD.add(S2,D,D);
	seq.push_back(S2);
	bmit++;
	check = bmit.state();
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	MD.addin(S2,D);
	seq.push_back(S2);
	bmit++;
	check = bmit.state();
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	if(check.IsGeneratorFound())
					report << "Generator Found" << endl;
	report << "mu = " << bmit.get_mu() << endl;
	report << "sigma = " << bmit.get_sigma() << endl;
	report << "beta = " << bmit.get_beta() << endl;
	BM_Seq<Field>::BM_iterator::TerminationState check2 = bmit2.state();
	while(!check2.IsGeneratorFound() && !check2.IsSequenceExceeded()){
					++bmit2;
					check2 = bmit2.state();
	}
	if(bmit==bmit2)
					report << "Iterators are equal" << endl;
	if(bmit2==seq.BM_end())
					report << "bmit2 is equal to end" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit2->begin(); it != bmit2->end(); it++)
					(*it).write(report);
	BM_Seq<Field>::BM_iterator bmit3 = seq.BM_begin();
	bmit3 = bmit;
	if(bmit==bmit3)
					report << "Iterators are equal" << endl;
	vector<BlasMatrix<Field> >gen(bmit.GetGenerator());
	int d = bmit.get_mu();
	for(int j = 0; j <= d; j++)

					gen[j].write(report);
#endif
#if 1
	size_t r,d,c;
	d=10;
	r=2;
	c=2;
	SparseMatrix<Field> A(F,d,d);
	BlasMatrix<Field> U(F,r,d);
	BlasMatrix<Field> W(F,d,c-1);
	BlasMatrix<Field> V(F,d,c);
	MatrixDomain<Field>::Submatrix V2(V,0,1,d,c-1);
	std::vector<Field::Element> b(d);
	Field::RandIter rand(F);
	for(size_t i=0; i<d;i++)
			A.setEntry(i,d-i-1,one);
	/*
	for(size_t i=0; i<d;i++)
		for(size_t j=0; j<d; j++)
			rand.random(A.refEntry(i,j));
	*/
	for(size_t i=0; i<r;i++)
		for(size_t j=0; j<d; j++)
			rand.random(U.refEntry(i,j));
	for(size_t i=0; i<d;i++)
		for(size_t j=0; j<c-1; j++){
			rand.random(W.refEntry(i,j));
	}
	MD.mul(V2,A,W);
	for(size_t i=0; i<d; i++){
		rand.random(b[i]);
		V.setEntry(i,0,b[i]);
}
	report << endl << "W" << endl;
	W.write(report);
	report << endl << "U" << endl;
	U.write(report);
	report << endl << "V" << endl;
	V.write(report);
	BlackboxBlockContainer<Field, SparseMatrix<Field> > blockseq(&A,F,U,V);
	BlackboxBlockContainer<Field, SparseMatrix<Field> >::const_iterator contiter(blockseq.begin());
	BM_Seq<Field> seq(F,r,c);
	seq.push_back(*contiter);
	BM_Seq<Field>::BM_iterator bmit(seq.BM_begin());
	bmit.setDelta(d);
	BM_Seq<Field>::BM_iterator::TerminationState check = bmit.state();
	while(!check.IsGeneratorFound() ){
		bmit++;
		check = bmit.state();
		if(check.IsSequenceExceeded()){
			++contiter;
			seq.push_back(*contiter);
		}
	}
	if(check.IsGeneratorFound())
		report << "Generator Found" << endl;
	report << "mu = " << bmit.get_mu() << endl;
	report << "sigma = " << bmit.get_sigma() << endl;
	report << "beta = " << bmit.get_beta() << endl;
	vector<BlasMatrix<Field> >gen(bmit.GetGenerator());
	int mu = bmit.get_mu();
	report << "The generator is " << endl;
	for(size_t i=0; i<mu+1; i++){
		report << "Degree " << i << endl;
		gen[i].write(report);
		report << endl;
}
	size_t idx = 0;
	if(F.isZero(gen[0].getEntry(0,0))){
		size_t i = 1;
		while(i<c && F.isZero(gen[0].getEntry(0,i)))
			i++;
		if(i==c)
			throw LinboxError(" block minpoly: matrix seems to be singular - abort");
		else
			idx=i;
	}
	report << "The nonzero index is " << idx << endl;
	
	BlasMatrix<Field> AV(V);
	BlasMatrix<Field> xm(F,d,1);
	for(size_t i = 1; i < mu+1; i++){
		MatrixDomain<Field>::Submatrix gencol(gen[i],0,idx,d,1);
		BlasMatrix<Field> AVgencol(F,d,1);
		MD.mul(AVgencol,AV,gencol);
		report << endl << "AVgencol at degree " << i << endl;
		AVgencol.write(report);
		MD.addin(xm, AVgencol);
		report << endl << "xm" << endl;
		xm.write(report);
		MD.leftMulin(A,AV);	
		report << endl << "AV" << endl;
		AV.write(report);
		report << endl;
	}

	for(size_t i = 1; i < c; i++){
		MatrixDomain<Field>::Submatrix Wcol(W,0,i-1,d,1);
		BlasMatrix<Field> Wcolgen0(F,d,1);
		MD.mul(Wcolgen0, Wcol, gen[0].getEntry(i,idx));
		MD.addin(xm,Wcolgen0);
		
	}
	MD.negin(xm);
	Field::Element gen0inv;
	MD.mulin(xm, F.inv(gen0inv, gen[0].getEntry(0,idx)));
	
	report << "The computed solution is" << endl;
	xm.write(report);		
	BlasMatrix<Field> Axm(F,r,1);
	BlasMatrix<Field> UAxm(F,r,1);
	MD.mul(Axm, A,xm);
	MD.mul(UAxm, U,Axm);
	report << "UA times solution is" << endl;
	UAxm.write(report);
	BlasMatrix<Field> Uym(F,r,1);
	MatrixDomain<Field>::Submatrix ym(V,0,0,d,1);
	MD.mul(Uym, U, ym);
	report << "U times rhs is" << endl;
	Uym.write(report);
	


#endif
	commentator().stop("bm-seq test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

