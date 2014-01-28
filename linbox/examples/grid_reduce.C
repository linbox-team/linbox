
/*
 * examples/grid_reduce.C
 *
 * Copyright (C) 2008, 2010 A. Urbanska
  ========LICENCE========
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
 */

/*! @file examples/grid_reduce.C
 * @example examples/grid_reduce.C
 * @ingroup examples
 * @brief undocumented
 */

#include <iostream>
#include <sstream>
#include<fstream>

#include <linbox/linbox-config.h>
#include <linbox/field/PID-integer.h>
#include <linbox/field/modular.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/grid.h>

using namespace LinBox;

int main(int argc, char* argv[])
{

	srand(time(NULL));

	if (argc < 5) {

		std::cout << "usage: " << argv[0] << " alg m n S source format " << std::endl  << std::endl;

		std::cout << "alg = `adaptive', `ilio', `local', or `2local', " << std::endl
		<< "m is modulus (ignored by 2local, adaptive), "
		<< "n is matrix order, " << std::endl
		<< "S is reduced rows size" << std::endl
		<< "source is `random', `random-rough', `fib', `tref', or a filename " << std::endl
		<< "format is `dense' or `sparse' (if matrix from a file)" << std::endl
		<< "compile with -DBIG if you want big integers used." << std::endl;

		return 0;
	}

	std::string algo = argv[1];

	// int m = atoi(argv[2]);

	// int n = atoi(argv[3]);

	int S = atoi(argv[4]);
	std::string src = argv[5];
	// std::string file = src;

	// std::string format = (argc >= 6 ? argv[5] : "");

	UserTimer T, TT;

	std::string out(src);
	out=out+"z";
	if (algo == "reduceT") out=out+"T";

	typedef PID_integer Ints;
	//typedef ModularBalanced<int> Ints;
	//typedef Modular<int> Ints;
	//Ints Z(m);
	Ints Z;

	integer p;
	std::cout << "Computation modulo " << Z.characteristic(p) << "" << std::endl<<std::flush;

	std::ifstream in (src.c_str(), std::ios::in);
	std::ofstream os (out.c_str(), std::ios::out);

	std::vector<int> mC;
	std::vector<int> mR;
#if 0
	   ifstream mRow ("mR2", ios::in);
	   char c='['; mRow >> c; std::cout << c;
	   while (c!='[') {
	   mRow >> c;
	   std::cout << c;
	   }
	   int mark;
	   while( mRow >> mark) {
	   mC.push_back(mark);
	   }
	   mRow.close();
	   mRow.open("mC2", ios::in) ;
	   c='['; mRow >> c; std::cout << c;
	   while (c!='[') {
	   mRow >> c;
	   std::cout << c;
	   }
	   for (int i=0; i< mC.size(); ++i) {
	   mRow >> mark;
	   if(mark==1) mC[i] = mark;
	   }
#endif
	std::cout << src << "" << std::endl<< std::flush ;
	TT.clear();
	TT.start();
	Grid<Ints, Ints::Element> A(Z,in, mR, mC);
	TT.stop();
	std::cout << "Reading the matrix in ";
	TT.print(std::cout);
	std::cout << " seconds" << std::endl << std::flush;

	int rank =0;
	TT.clear();
	TT.start();
	A.reduce(rank, S, mR, mC, os);
	TT.stop();
	std::cout << "Reducing the matrix in ";
	TT.print(std::cout);
	std::cout << " seconds" << std::endl << std::flush;

	std::cout << "The rank is " << rank << "" << std::endl << std::flush;
	A.write(os);
	//std::cout << "mR: \n[";
	//for (int i=0; i < mR.size() ; ++i) {
	//	std::cout << mR[i] << " ";
	//}
	os << "mC: \n[";
	for (int i=0; i < mC.size() ; ++i) {
		os << mC[i] << " ";
	}
	os << "]" << std::endl;
	os << "mR: \n[";
	for (int i=0; i < mR.size() ; ++i) {
		os<< mR[i] << " ";
	}
	os <<std::endl <<std::flush;


	int matrix=1;
	TT.clear();
	TT.start();
	while  (matrix+1 < argc-5) {//argc - "grid-reduce reduce m n S", while takes 2 matrices
		TT.stop();
		std::cout << "Time other:";
		TT.print(std::cout);
		std::cout << "seconds" << std::endl << std::flush;
		TT.clear();
		TT.start();
		//mC = mR;	//mR plays the role of mC
		mC.resize(0);
		char c='['; in >> c; std::cout << c;
		while (c!='[') {
			in >> c;
			std::cout << c;
		}
		TT.stop();
		std::cout << "Time vector exchange in";
		TT.print(std::cout);
		std::cout << "seconds" << std::endl << std::flush;
		TT.clear();
		TT.start();
		std::cout << "reading additional columns" << std::endl << std::flush;
		for (int i=0; i < mR.size(); ++i) {
			int mark;
			in >> mark;
			if (mR[i]!=1) mR[i] = mark;
			//if (mark>0) std::cout << mark << " ";
		}
		std::cout << "columns read in ";
		TT.stop();
		TT.print(std::cout);
		std::cout << "seconds" << std::endl  << std::flush;

		std::cout << "Time file preparation in";
		TT.clear();
		TT.start();
		in.close();
		src = argv[5+matrix];
		std::cout << src << "" << std::endl << std::flush;
		std::string out3(src);
		out3= out3+"z";
		if (algo == "reduceT") out3=out3+"T";
		//++matrix;
		in.open(src.c_str(), std::ios::in);
		os.close();
		os.open(out3.c_str(), std::ios::out);
		TT.stop();
		TT.print(std::cout);
		std::cout << "seconds" << std::endl  << std::flush;
		TT.clear();
		TT.start();
		A.read(Z,in, mC, mR);
		//Grid<Ints, Ints::Element> A2(Z,in, mR, mC);
		TT.stop();
		std::cout << "Reading the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

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
		std::cout << "Reducing the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		std::cout << "The rank is " << rank << "" << std::endl << std::flush;
		TT.clear();
		TT.start();
		A.write(os);
		TT.stop();
		std::cout << "Writing the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		TT.clear();
		TT.start();
		os << "mC:\n[";
		for (int i=0; i < mR.size() ; ++i) {
			os << mR[i] << " ";
		}
		os << "" << std::endl;
		os << "mR: \n[";
		for (int i=0; i < mC.size() ; ++i) {
			os<< mC[i] << " ";
		}
		os <<std::endl <<std::flush;
		TT.stop();
		std::cout << "Writing mC vector in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		//mC = mR;
		mR.resize(0);
		c='['; in >> c; std::cout << c;
		while (c!='[') {
			in >> c;
			std::cout << c;
		}
		TT.stop();
		std::cout << "Time vector exchange in";
		TT.print(std::cout);
		std::cout << "seconds" << std::endl << std::flush;
		TT.clear();
		TT.start();
		std::cout << "reading additional columns" << std::endl << std::flush;
		for (int i=0; i < mC.size(); ++i) {
			int mark;
			in >> mark;
			if (mC[i]!=1) mC[i] = mark;
			//if (mark>0) std::cout << mark << " ";
		}
		std::cout << "columns read in ";
		TT.stop();
		TT.print(std::cout);
		std::cout << "seconds" << std::endl  << std::flush;

		std::cout << "Time file preparation in";
		TT.clear();
		TT.start();
		in.close();
		src = argv[6+matrix];
		std::cout << src << "" << std::endl << std::flush;
		std::string out4(src);
		out4= out4+"z";
		if (algo == "reduceT") out4=out4+"T";
		matrix+=2;
		in.open(src.c_str(), std::ios::in);
		os.close();
		os.open(out4.c_str(), std::ios::out);
		TT.stop();
		TT.print(std::cout);
		std::cout << "seconds" << std::endl  << std::flush;
		TT.clear();
		TT.start();

		A.read(Z,in, mR, mC);
		//Grid<Ints, Ints::Element> A2(Z,in, mR, mC);
		TT.stop();
		std::cout << "Reading the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		int rank =0;
		TT.clear();
		TT.start();
		if (matrix==3) A.reduce(rank, 6.8*S, mR, mC, os); //21(3,9)
		if (matrix==9) A.reduce(rank, 6.8*S, mR, mC, os); //21(3,9)
		if (matrix==7) A.reduce(rank, 4.6*S, mR, mC, os); //19(5,7)nochange1.6
		A.reduce(rank, S, mR, mC, os);
		TT.stop();
		std::cout << "Reducing the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		std::cout << "The rank is " << rank << "" << std::endl << std::flush;
		TT.clear();
		TT.start();
		A.write(os);
		TT.stop();
		std::cout << "Writing the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		TT.clear();
		TT.start();
		os << "mC:\n[";
		for (int i=0; i < mC.size() ; ++i) {
			os << mC[i] << " ";
		}
		os << "" << std::endl;
		os << "mR: \n[";
		for (int i=0; i < mR.size() ; ++i) {
			os<< mR[i] << " ";
		}
		os <<std::endl <<std::flush;
		TT.stop();
		std::cout << "Writing mC vector in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		TT.clear();
		TT.start();
	}

	if (matrix < argc-5) {
		TT.stop();
		std::cout << "Time other:";
		TT.print(std::cout);
		std::cout << "seconds" << std::endl << std::flush;
		TT.clear();
		TT.start();
		//mC = mR;      //mR plays the role of mC
		mC.resize(0);
		char c='['; in >> c; std::cout << c;
		while (c!='[') {
			in >> c;
			std::cout << c;
		}
		TT.stop();
		std::cout << "Time vector exchange in";
		TT.print(std::cout);
		std::cout << "seconds" << std::endl << std::flush;
		TT.clear();
		TT.start();
		std::cout << "reading additional columns" << std::endl << std::flush;
		for (int i=0; i < mR.size(); ++i) {
			int mark;
			in >> mark;
			if (mR[i]!=1) mR[i] = mark;
			//if (mark>0) std::cout << mark << " ";
		}
		std::cout << "columns read in ";
		TT.stop();
		TT.print(std::cout);
		std::cout << "seconds" << std::endl  << std::flush;

		std::cout << "Time file preparation in";
		TT.clear();
		TT.start();
		in.close();
		src = argv[5+matrix];
		std::cout << src << "" << std::endl << std::flush;
		std::string out2(src);
		out2= out2+"z";
		if (algo == "reduceT") out2=out2+"T";
		//++matrix;
		in.open(src.c_str(), std::ios::in);
		os.close();
		os.open(out2.c_str(), std::ios::out);
		TT.stop();
		TT.print(std::cout);
		std::cout << "seconds" << std::endl  << std::flush;
		TT.clear();
		TT.start();
		A.read(Z,in, mC, mR);
		//Grid<Ints, Ints::Element> A2(Z,in, mR, mC);
		TT.stop();
		std::cout << "Reading the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		int rank =0;
		TT.clear();
		TT.start();
		A.reduce(rank, S, mC, mR, os);
		TT.stop();
		std::cout << "Reducing the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		std::cout << "The rank is " << rank << "" << std::endl << std::flush;
		TT.clear();
		TT.start();
		A.write(os);
		TT.stop();
		std::cout << "Writing the matrix in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

		TT.clear();
		TT.start();
		os << "mC:\n[";
		for (int i=0; i < mR.size() ; ++i) {
			os << mR[i] << " ";
		}
		os << "" << std::endl;
		os << "mC: \n[";
		for (int i=0; i < mC.size() ; ++i) {
			os<< mC[i] << " ";
		}
		os <<std::endl <<std::flush;
		TT.stop();
		std::cout << "Writing mC vector in ";
		TT.print(std::cout);
		std::cout << " seconds" << std::endl << std::flush;

	}
	TT.stop();

	return 0 ;
}


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
