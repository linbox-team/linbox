/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * examples/grid_reduce.C
 *
 * Copyright (C) 2008, 2010 A. Urbanska
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

/*! @file examples/grid_reduce.C
 * @example examples/grid_reduce.C
 * @ingroup examples
 * @brief undocumented
 */

#include <iostream>
#include <sstream>
#include<fstream>

#include "linbox/field/PID-integer.h"
#include "linbox/field/modular-balanced-int.h"
#include "linbox/field/modular-int.h"
#include "linbox/util/timer.h"
#include "linbox-config.h"
#include "linbox/matrix/grid.h"

using namespace std;
using namespace LinBox;

int main(int argc, char* argv[])
{

	srand(time(NULL));

	if (argc < 5) {

		cout << "usage: " << argv[0] << " alg m n S source format \n"  << endl;

		cout << "alg = `adaptive', `ilio', `local', or `2local', \n"
		<< "m is modulus (ignored by 2local, adaptive), "
		<< "n is matrix order, \n"
		<< "S is reduced rows size\n"
		<< "source is `random', `random-rough', `fib', `tref', or a filename \n"
		<< "format is `dense' or `sparse' (if matrix from a file)\n"
		<< "compile with -DBIG if you want big integers used.\n";

		return 0;
	}

	string algo = argv[1];

	int m = atoi(argv[2]);

	int n = atoi(argv[3]);

	int S = atoi(argv[4]);
	string src = argv[5];
	string file = src;

	string format = (argc >= 6 ? argv[5] : "");

	UserTimer T, TT;

	string out(src);
	out=out+"z";
	if (algo == "reduceT") out=out+"T";

	typedef PID_integer Ints;
	//typedef ModularBalanced<int> Ints;
	//typedef Modular<int> Ints;
	//Ints Z(m);
	Ints Z;

	integer p;
	cout << "Computation modulo " << Z.characteristic(p) << "\n"<<flush;

	ifstream in (src.c_str(), ios::in);
	ofstream os (out.c_str(), ios::out);

	std::vector<int> mC;
	std::vector<int> mR;
#if 0
	   ifstream mRow ("mR2", ios::in);
	   char c='['; mRow >> c; cout << c;
	   while (c!='[') {
	   mRow >> c;
	   cout << c;
	   }
	   int mark;
	   while( mRow >> mark) {
	   mC.push_back(mark);
	   }
	   mRow.close();
	   mRow.open("mC2", ios::in) ;
	   c='['; mRow >> c; cout << c;
	   while (c!='[') {
	   mRow >> c;
	   cout << c;
	   }
	   for (int i=0; i< mC.size(); ++i) {
	   mRow >> mark;
	   if(mark==1) mC[i] = mark;
	   }
#endif
	cout << src << "\n"<< flush ;
	TT.clear();
	TT.start();
	Grid<Ints, Ints::Element> A(Z,in, mR, mC);
	TT.stop();
	cout << "Reading the matrix in ";
	TT.print(cout);
	cout << " seconds\n" << flush;

	int rank =0;
	TT.clear();
	TT.start();
	A.reduce(rank, S, mR, mC, os);
	TT.stop();
	cout << "Reducing the matrix in ";
	TT.print(cout);
	cout << " seconds\n" << flush;

	cout << "The rank is " << rank << "\n" << flush;
	A.write(os);
	//cout << "mR: \n[";
	//for (int i=0; i < mR.size() ; ++i) {
	//	cout << mR[i] << " ";
	//}
	os << "mC: \n[";
	for (int i=0; i < mC.size() ; ++i) {
		os << mC[i] << " ";
	}
	os << "]\n";
	os << "mR: \n[";
	for (int i=0; i < mR.size() ; ++i) {
		os<< mR[i] << " ";
	}
	os <<endl <<flush;


	int matrix=1;
	TT.clear();
	TT.start();
	while  (matrix+1 < argc-5) {//argc - "grid-reduce reduce m n S", while takes 2 matrices
		TT.stop();
		cout << "Time other:";
		TT.print(cout);
		cout << "seconds\n" << flush;
		TT.clear();
		TT.start();
		//mC = mR;	//mR plays the role of mC
		mC.resize(0);
		char c='['; in >> c; cout << c;
		while (c!='[') {
			in >> c;
			cout << c;
		}
		TT.stop();
		cout << "Time vector exchange in";
		TT.print(cout);
		cout << "seconds\n" << flush;
		TT.clear();
		TT.start();
		cout << "reading additional columns\n" << flush;
		for (int i=0; i < mR.size(); ++i) {
			int mark;
			in >> mark;
			if (mR[i]!=1) mR[i] = mark;
			//if (mark>0) cout << mark << " ";
		}
		cout << "columns read in ";
		TT.stop();
		TT.print(cout);
		cout << "seconds\n"  << flush;

		cout << "Time file preparation in";
		TT.clear();
		TT.start();
		in.close();
		src = argv[5+matrix];
		cout << src << "\n" << flush;
		string out3(src);
		out3= out3+"z";
		if (algo == "reduceT") out3=out3+"T";
		//++matrix;
		in.open(src.c_str(), std::ios::in);
		os.close();
		os.open(out3.c_str(), std::ios::out);
		TT.stop();
		TT.print(cout);
		cout << "seconds\n"  << flush;
		TT.clear();
		TT.start();
		A.read(Z,in, mC, mR);
		//Grid<Ints, Ints::Element> A2(Z,in, mR, mC);
		TT.stop();
		cout << "Reading the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		rank =0;
		TT.clear();
		TT.start();
		if (matrix==1) A.reduce(rank, 1000000*S, mC, mR, os); //22(1,9)
		if (matrix==9) A.reduce(rank, 1000000*S, mC, mR, os); //22(1,9)
		if (matrix==7) A.reduce(rank, 4.7*S, mC, mR, os); //20(3,7)change1.6
		if (matrix==5) A.reduce(rank, 1000000*S, mC, mR, os); //18(5,5)change3
		//if (matrix==1) A.reduce(rank, S,mC,mR,os);
		else A.reduce(rank, S, mC, mR, os);
		TT.stop();
		cout << "Reducing the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		cout << "The rank is " << rank << "\n" << flush;
		TT.clear();
		TT.start();
		A.write(os);
		TT.stop();
		cout << "Writing the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		TT.clear();
		TT.start();
		os << "mC:\n[";
		for (int i=0; i < mR.size() ; ++i) {
			os << mR[i] << " ";
		}
		os << "\n";
		os << "mR: \n[";
		for (int i=0; i < mC.size() ; ++i) {
			os<< mC[i] << " ";
		}
		os <<endl <<flush;
		TT.stop();
		cout << "Writing mC vector in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		//mC = mR;
		mR.resize(0);
		c='['; in >> c; cout << c;
		while (c!='[') {
			in >> c;
			cout << c;
		}
		TT.stop();
		cout << "Time vector exchange in";
		TT.print(cout);
		cout << "seconds\n" << flush;
		TT.clear();
		TT.start();
		cout << "reading additional columns\n" << flush;
		for (int i=0; i < mC.size(); ++i) {
			int mark;
			in >> mark;
			if (mC[i]!=1) mC[i] = mark;
			//if (mark>0) cout << mark << " ";
		}
		cout << "columns read in ";
		TT.stop();
		TT.print(cout);
		cout << "seconds\n"  << flush;

		cout << "Time file preparation in";
		TT.clear();
		TT.start();
		in.close();
		src = argv[6+matrix];
		cout << src << "\n" << flush;
		string out4(src);
		out4= out4+"z";
		if (algo == "reduceT") out4=out4+"T";
		matrix+=2;
		in.open(src.c_str(), std::ios::in);
		os.close();
		os.open(out4.c_str(), std::ios::out);
		TT.stop();
		TT.print(cout);
		cout << "seconds\n"  << flush;
		TT.clear();
		TT.start();

		A.read(Z,in, mR, mC);
		//Grid<Ints, Ints::Element> A2(Z,in, mR, mC);
		TT.stop();
		cout << "Reading the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		int rank =0;
		TT.clear();
		TT.start();
		if (matrix==3) A.reduce(rank, 6.8*S, mR, mC, os); //21(3,9)
		if (matrix==9) A.reduce(rank, 6.8*S, mR, mC, os); //21(3,9)
		if (matrix==7) A.reduce(rank, 4.6*S, mR, mC, os); //19(5,7)nochange1.6
		A.reduce(rank, S, mR, mC, os);
		TT.stop();
		cout << "Reducing the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		cout << "The rank is " << rank << "\n" << flush;
		TT.clear();
		TT.start();
		A.write(os);
		TT.stop();
		cout << "Writing the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		TT.clear();
		TT.start();
		os << "mC:\n[";
		for (int i=0; i < mC.size() ; ++i) {
			os << mC[i] << " ";
		}
		os << "\n";
		os << "mR: \n[";
		for (int i=0; i < mR.size() ; ++i) {
			os<< mR[i] << " ";
		}
		os <<endl <<flush;
		TT.stop();
		cout << "Writing mC vector in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		TT.clear();
		TT.start();
	}

	if (matrix < argc-5) {
		TT.stop();
		cout << "Time other:";
		TT.print(cout);
		cout << "seconds\n" << flush;
		TT.clear();
		TT.start();
		//mC = mR;      //mR plays the role of mC
		mC.resize(0);
		char c='['; in >> c; cout << c;
		while (c!='[') {
			in >> c;
			cout << c;
		}
		TT.stop();
		cout << "Time vector exchange in";
		TT.print(cout);
		cout << "seconds\n" << flush;
		TT.clear();
		TT.start();
		cout << "reading additional columns\n" << flush;
		for (int i=0; i < mR.size(); ++i) {
			int mark;
			in >> mark;
			if (mR[i]!=1) mR[i] = mark;
			//if (mark>0) cout << mark << " ";
		}
		cout << "columns read in ";
		TT.stop();
		TT.print(cout);
		cout << "seconds\n"  << flush;

		cout << "Time file preparation in";
		TT.clear();
		TT.start();
		in.close();
		src = argv[5+matrix];
		cout << src << "\n" << flush;
		string out2(src);
		out2= out2+"z";
		if (algo == "reduceT") out2=out2+"T";
		//++matrix;
		in.open(src.c_str(), std::ios::in);
		os.close();
		os.open(out2.c_str(), std::ios::out);
		TT.stop();
		TT.print(cout);
		cout << "seconds\n"  << flush;
		TT.clear();
		TT.start();
		A.read(Z,in, mC, mR);
		//Grid<Ints, Ints::Element> A2(Z,in, mR, mC);
		TT.stop();
		cout << "Reading the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		int rank =0;
		TT.clear();
		TT.start();
		A.reduce(rank, S, mC, mR, os);
		TT.stop();
		cout << "Reducing the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		cout << "The rank is " << rank << "\n" << flush;
		TT.clear();
		TT.start();
		A.write(os);
		TT.stop();
		cout << "Writing the matrix in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

		TT.clear();
		TT.start();
		os << "mC:\n[";
		for (int i=0; i < mR.size() ; ++i) {
			os << mR[i] << " ";
		}
		os << "\n";
		os << "mC: \n[";
		for (int i=0; i < mC.size() ; ++i) {
			os<< mC[i] << " ";
		}
		os <<endl <<flush;
		TT.stop();
		cout << "Writing mC vector in ";
		TT.print(cout);
		cout << " seconds\n" << flush;

	}
	TT.stop();

	return 0 ;
}

