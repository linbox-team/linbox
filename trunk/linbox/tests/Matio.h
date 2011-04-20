/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) LinBox
 *
 *
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

#ifndef __LINBOX_matio_H
#define __LINBOX_matio_H

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include "fflas-ffpack/fflas/fflas.h"



// Reading and writing matrices over double

#if 0
// Reading a matrice from a (eventually zipped) file
double * read_dbl(char * mat_file,int* tni,int* tnj)
{
	char *UT, *File_Name;
	int is_gzipped = 0;
	size_t s = strlen(mat_file);
	double* X;
	if ((mat_file[--s] == 'z') &&
	    (mat_file[--s] == 'g') &&
	    (mat_file[--s] == '.')) {
		is_gzipped = 1;
		File_Name = "/tmp/bbXXXXXX_";
		mkstemp(File_Name);
		UT = new char[s+34+strlen(File_Name)];
		sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
		system(UT);
		sprintf(UT,"\\rm %s", File_Name);
	} else
		File_Name = mat_file;

	FILE* FileDes = fopen(File_Name, "r");
	if (FileDes != NULL) {
		char * tmp = new char[200];// usigned long tni, tnj;
		fscanf(FileDes,"%d %d %s\n",tni, tnj, &tmp) ;
		int n=*tni;
		int p=*tnj;
		X = new double[n*p];
		for (int i=0;i<n*p;++i)
			X[i] = (double) 0;
		long i,j; long val;
		fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
		while(i && j) {
			X[p*(i-1)+j-1] = (double) val;
			fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
		}
	}

	fclose(FileDes);
	if (is_gzipped) system(UT);
	return X;
}

// Displays a matrix
std::ostream& write_dbl(std::ostream& c,
			double* E,
			int n, int m, int id)
{

	for (int i = 0; i<n;++i){
		for (int j=0; j<m;++j)
			c << *(E+j+id*i) << " ";
		c << std::endl;
	}
	return c << std::endl;
}
#endif
// Reading and writing matrices over field

// Reading a matrice from a (eventually zipped) file
template<class Field>
typename Field::Element * read_field(const Field& F,char * mat_file,int* tni,int* tnj)
{
	char *UT, *File_Name;
	int is_gzipped = 0;
	size_t s = strlen(mat_file);
	typename Field::Element zero;
	F.init(zero,0);
	typename Field::Element * X;
	if ((mat_file[--s] == 'z') &&
	    (mat_file[--s] == 'g') &&
	    (mat_file[--s] == '.')) {
		is_gzipped = 1;
		File_Name = "/tmp/bbXXXXXX_";
		mkstemp(File_Name);
		UT = new char[s+34+strlen(File_Name)];
		sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
		system(UT);
		sprintf(UT,"\\rm %s", File_Name);
	} else
		File_Name = mat_file;
	FILE* FileDes = fopen(File_Name, "r");
	if (FileDes != NULL) {
		char * tmp = new char[200];// usigned long tni, tnj;
		fscanf(FileDes,"%d %d %s\n",tni, tnj, &tmp) ;
		int n=*tni;
		int p=*tnj;
		X = new typename Field::Element[n*p];
		for (int i=0;i<n*p;++i)
			X[i] = zero;
		long i,j; long val;
		fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
		while(i && j) {
			F.init(X[p*(i-1)+j-1],val);
			fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
		}
	}

	fclose(FileDes);
	if (is_gzipped) system(UT);
	return X;
}

template<class Field>
void read_field4(const Field& F,char * mat_file,int* tni,int* tnj,
		 typename Field::Element *& NW,typename Field::Element *& NE,
		 typename Field::Element *& SW,typename Field::Element *& SE)
{
	char *UT, *File_Name;
	int is_gzipped = 0;
	size_t s = strlen(mat_file);
	typename Field::Element zero;
	F.init(zero,0);
	typename Field::Element * X;
	if ((mat_file[--s] == 'z') &&
	    (mat_file[--s] == 'g') &&
	    (mat_file[--s] == '.')) {
		is_gzipped = 1;
		File_Name = "/tmp/bbXXXXXX_";
		mkstemp(File_Name);
		UT = new char[s+34+strlen(File_Name)];
		sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
		system(UT);
		sprintf(UT,"\\rm %s", File_Name);
	} else
		File_Name = mat_file;
	FILE* FileDes = fopen(File_Name, "r");
	if (FileDes != NULL) {
		char * tmp = new char[200];// usigned long tni, tnj;
		fscanf(FileDes,"%d %d %s\n",tni, tnj, &tmp) ;
		int n=*tni;
		int p=*tnj;
		int no2= n>>1;
		int po2 = p>>1;
		NW = new typename Field::Element[no2*po2];
		NE = new typename Field::Element[no2*(p-po2)];
		SW = new typename Field::Element[(n-no2)*po2];
		SE = new typename Field::Element[(n-no2)*(p-po2)];

		for (int i=0;i<no2*po2;++i)
			NW[i] = zero;
		for (int i=0;i<no2*(p-po2);++i)
			NE[i] = zero;
		for (int i=0;i<(n-no2)*po2;++i)
			SW[i] = zero;
		for (int i=0;i<(n-no2)*(p-po2);++i)
			SE[i] = zero;
		long i,j; long val;
		fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
		while(i && j) {
			if (i<=no2){
				if (j<=po2){
					F.init(NW[po2*(i-1)+j-1],val);
					fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
				}
				else{
					F.init(NE[po2*(i-1)+j-1-po2],val);
					fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
				}
			}
			else{
				if (j<=po2){
					F.init(SW[(p-po2)*(i-1-no2)+j-1],val);
					fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
				}
				else{
					F.init(SE[(p-po2)*(i-1-no2)+j-1-po2],val);
					fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
				}
			}
		}
		//    *A1 = NW;
		//*A2 = NE;
		//*A3 = SW;
		//*A4 = SE;

	}

	fclose(FileDes);
	if (is_gzipped) system(UT);
}

// Displays a matrix
template<class Field>
std::ostream& write_field(const Field& F,std::ostream& c,
			  const typename Field::Element* E,
			  int n, int m, int id, bool mapleFormat = false)
{

	double tmp;
	if (mapleFormat) c << "Matrix(" << n <<',' << m << ", [" ;
	for (int i = 0; i<n;++i){
		if (mapleFormat) c << '[';
		for (int j=0; j<m;++j){
			F.convert(tmp,*(E+j+id*i));
			c << tmp;
			if (mapleFormat && j<m-1) c << ',';
			c << ' ';
		}
		if (mapleFormat) c << ']';
		if (mapleFormat && i<n-1) c << ',';
		if (!mapleFormat) c << std::endl;
	}
	if (mapleFormat) c << "])";
	return c ;
}

// Displays a triangular matrix
//! @todo let the user choose to convert to a non destructive format (not double but long or Integer...)
template<class Field>
std::ostream& write_field(const Field& F,std::ostream& c,
			  const LinBox::FFLAS::FFLAS_UPLO uplo, const LinBox::FFLAS::FFLAS_DIAG unit,
			  const typename Field::Element* E,
			  int n, int m, int id, bool mapleFormat = false)
{

	double tmp;
	if (mapleFormat) c << "Matrix(" << n <<',' << m << ",[";
	for (int i = 0; i<n;++i){
		if (mapleFormat) c << '[';
		// under diag
		for (int j=0; j<i ;++j){
			if (uplo == LinBox::FFLAS::FflasLower)
				F.convert(tmp,*(E+j+id*i));
			else tmp = 0 ;
			c << tmp;
			if (mapleFormat && j<m-1) c << ',';
			c << ' ';
		}
		// on diag
		if (unit == LinBox::FFLAS::FflasNonUnit)
			F.convert(tmp,*(E+i+id*i));
		else
			tmp = 1.;
		c << tmp;
		if (mapleFormat && i<m-1) c << ',';
		c << ' ';
		// over diag
		for (int j=i+1; j<m;++j){
			if (uplo == LinBox::FFLAS::FflasUpper)
				F.convert(tmp,*(E+j+id*i));
			else
				tmp = 0 ;
			c << tmp;
			if (mapleFormat && j<m-1) c << ',';
			c << ' ';
		}
		if (mapleFormat) c << ']';
		if (mapleFormat && i<n-1) c << ',';
		if (!mapleFormat) c << std::endl;
	}
	if (mapleFormat) c << "])";
	return c ;
}

#endif //__LINBOX_matio_H
