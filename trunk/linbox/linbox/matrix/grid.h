/* Copyright (C) LinBox
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_matrix_grid_H
#define __LINBOX_matrix_grid_H

#include <iostream>

#include <vector>
#include <queue>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"


namespace LinBox
{

	template <class Element>
	class ijElement {

	public:
		size_t i;
		size_t j;
		Element x;

		ijElement(size_t ii, size_t jj, Element xx) {i=ii;j=jj;x=xx;}
		ijElement(ijElement& Y) {i = Y.i; j = Y.j; x = Y.x;}
		~ijElement() {};
		//Element& setElement(Element& y) {return x=y;}
	};

	template <class Element>
	class GridElement {

	public:
		GridElement* prev;
		GridElement* next;
		GridElement* down;
		GridElement* up;
		ijElement<Element> X;

		GridElement(ijElement<Element>& Y) :
			X(Y)
		{
			prev = NULL;
			next = NULL;
			down = NULL;
			up = NULL;
		}

		GridElement(size_t i, size_t j, Element x) :
			X(i,j,x)
		{
			prev = NULL;
			next = NULL;
			down = NULL;
			up = NULL;
		}

		~GridElement()
		{
			//prev = NULL; next = NULL; down = NULL; up = NUll;
		}

		Element& getX() {return X.x;}
		size_t getI() {return X.i;}
		size_t getJ() {return X.j;}

		Element& setElement(Element& y) {return X.x=y;}
		ijElement<Element>& setijElement(size_t i, size_t j, Element& x) {X.i=i; X.j = j; X.x = x; return X;}
	};

	template <class Field, class Element>
	class Grid {

	public:
		Field _field;
		size_t _n;
		size_t _m;
		std::vector<int> rowOcc;
		//std::vector<Element> rowGcd;
		std::vector<int> colOcc;
		//std::vector<Element> colGcd;

		std::queue<size_t> Q; //change to priority queue?
		std::vector<GridElement<Element>*> A; 	//std::vector of row heads
		std::vector<GridElement<Element>*> AT;	//std::vector of column heads

		/*
		 * Creates Grid from file using read procedure
		 * file is sorted in sms format
		 */
		Grid (Field F, std::istream& in, std::vector<int>& mR, std::vector<int>& mC) :
			_field(F)
		{
			read(F, in, mR, mC);
			in >> _m >> _n;
		}


		/*
		 * reads the matrix from the file, deletes mC[i]=1 columns
		 * sets row/col Occ and Gcd
		 * file is sorted in sms format
		 */
		void read(Field F, std::istream& in, std::vector<int>& mR, std::vector<int>& mC)
		{
			in >> _m >> _n;
			/* !!! */
			//int t= _m;
			//_m = _n;
			//_n = t;
			rowOcc.resize(_m,0);
			colOcc.resize(_n,0);
			//rowGcd.resize(_m,0);
			//colGcd.resize(_n,0);
			A.resize(_m, NULL);
			AT.resize(_n, NULL);

			mR.resize(_m);
			mC.resize(_n);
			char c;
			do in >> c; while (isspace (c));

			std::vector<GridElement<Element>*> ends(_n, NULL);

			if (c == 'M') {
				// size_t i, i_prev, j, j_prev;
				size_t i = 0, j = 0 ;
				// i_prev = 0;
				// j_prev = 0;

				Element x, x_prev; x_prev = 0; x=0;
				while (in >> i) {
					in >> j;
					/* !!! */
					//int t = i; i=j; j=t;
					if ( (i > _m) || (j > _n) ) {
						std::cout << "InvalidMatrixInput \n"<< std::flush;
						return;
					}
					F.read(in, x);
					if ((j > 0) && (F.isZero(x))) continue;
					if ((j > 0) && (mC[j-1]==1)){
						continue;       //mC are to be ignored;
					}
					if ((i!=0) && (j !=0 )) {
						ijElement<Element> X(i-1,j-1,x);
						ends[j-1] = addElement(ends[j-1], X);
					}
					if ((i==0) && (j ==0 )) break;
				}
				for (int k=0; k < rowOcc.size(); ++k) {
					if ((rowOcc[k]==1) && (/*colGcd[A[k]->getJ()]*/1==abs(A[k]->getX()))) {
						//std::cout << "1 row" << i_prev << "\n";
						//std::cout << "Adding " << k+1 << " to Q\n" << std::flush;
						Q.push(k);
					}
				}
			}
			mC.clear();
			mC.resize(_n,0);

		}

		/*
		 * deletes the non-null element Aij
		 */
		GridElement<Element>* deleteElement(GridElement<Element>* Aij)
		{
			if (Aij->prev != NULL) Aij->prev->next = Aij->next;
			else A[Aij->getI()] = Aij->next;
			if (Aij->next != NULL) Aij->next->prev = Aij->prev;
			if (Aij->down != NULL) Aij->down->up = Aij->up;
			else AT[Aij->getJ()] = Aij->up;
			if (Aij->up != NULL)   Aij->up->down = Aij->down;
			--rowOcc[Aij->getI()];
			if (rowOcc[Aij->getI()]==1) {
				//std::cout << "Adding " << Aij->getI()+1 << " to Q\n" << std::flush;
				if (1/*colGcd[A[Aij->getI()]->getJ()]*/ == abs (A[Aij->getI()]->getX()))
					Q.push(Aij->getI());
			}
			--colOcc[Aij->getJ()];
			GridElement<Element>* tmp = Aij->up;
			delete Aij;
			return tmp;
		}

		void deleteColumn (size_t j)
		{
			colOcc[j] = 0;
			GridElement<Element>* tmp = AT[j];
			while (tmp != NULL) {
				--rowOcc[tmp->getI()];
				if (rowOcc[tmp->getI()] ==1) {
					//std::cout << "Adding " << tmp->getI()+1 << " to Q\n" << std::flush;
					if ( 1/*colGcd[A[tmp->getI()]->getJ()]*/ == abs (A[tmp->getI()]->getX()))
						Q.push(tmp->getI());
				}
				if (tmp->prev != NULL) tmp->prev->next = tmp->next;
				else A[tmp->getI()] = tmp->next;
				if (tmp->next != NULL) tmp->next->prev = tmp->prev;
				GridElement<Element>* tmp2 = tmp->up;
				delete tmp;
				tmp = tmp2;
			}
		}

		void deleteRow (size_t i)
		{
			rowOcc[i] = 0;
			GridElement<Element>* tmp = A[i];
			while (tmp != NULL) {
				--colOcc[tmp->getJ()];
				if (tmp->down != NULL) tmp->down->up = tmp->up;
				else AT[tmp->getJ()] = tmp->up;
				if (tmp->up != NULL) tmp->up->down = tmp->down;
				GridElement<Element>* tmp2 = tmp->next;
				delete tmp;
				tmp = tmp2;
			}
		}

		/*
		 * adds a new Grid Element after (at ''up'') the given GridElement
		 * if lower==NULL Y should be the lowest element of the column
		 * returns a pointer to the new element
		 */
		GridElement<Element>* addElement(GridElement<Element>* lower, ijElement<Element>& Y)
		{
			++rowOcc[Y.i];
			++colOcc[Y.j];
			//rowGcd[Y.i]=gcd(rowGcd[Y.i], Y.x);
			//colGcd[Y.j]=gcd(colGcd[Y.j], Y.x);
			GridElement<Element>* X = new GridElement<Element>(Y);
			if (lower != NULL) {
				GridElement<Element>* tmp = lower->up;
				lower->up = X;
				X->down = lower;
				X->up = tmp;
				if (tmp != NULL) tmp->down = X;
				tmp = A[Y.i];
				A[Y.i] = X;
				X->next = tmp;
				if (tmp != NULL) tmp->prev = X;
			}
			else {
				GridElement<Element>* tmp = AT[Y.j];
				AT[Y.j] = X;
				X->up = tmp;
				if (tmp != NULL) tmp->down = X;
				tmp = A[Y.i];
				A[Y.i] = X;
				X->next = tmp;
				if (tmp != NULL) tmp->prev = X;
			}
			return X;
		}

		/*
		 * recursive procedure to reduce the grid
		 * returns precomputed rank
		 */

		int reduce(int& rank, int S, std::vector<int>& mR, std::vector<int>& mC, std::ostream& os)
		{
			std::cout << "rank at begin reduce " << rank << "\n" << std::flush;
			while (!Q.empty()) {
				size_t i=Q.front();
				Q.pop();
				if (rowOcc[i]==1) {//if row not deleted
					size_t j = A[i]->getJ();
					Element x; _field.init(x, A[i]->getX());
					if (mC[j]==1) {//if col deleted (not needed as we delete at once)
						mR[i]=2;
					}
					else if (abs(x)==1/*colGcd[j]*/) {
						//std::cout << "Row/column "<< i+1 << "," << j+1<< "," << x << "to reduce\n" <<std::flush;
						if (abs(x) > 1) std::cout << "adds " << x << "to the diagonal\n"<< std::flush;
						if (mC[j] !=1 ) {//if col not deleted then delete, update rank
							++rank;
							mC[j]=1;
							mR[i]=1;
							deleteColumn(j);//do not mark 0 rows
							//deleteRow(i);//not needed
						}
					} //else std::cout << "NOT Row/column "<< i+1 << "," << j+1<< "," << x << "to reduce\n" <<std::flush;
				}
			}

			std::cout << "Rank at end reduce/begin elimination" << rank <<"\n" << std::flush;

			size_t ini=0;

			size_t row2;
			size_t j1,j2;
			Element x,x1,x2;
			x1=1;x2=1;x=0;

			while (1) {
			bool pivotFound;
				pivotFound =false;
				if (ini>=mR.size()) ini=0;
				for (size_t i=ini; i <  mR.size(); ++i) {
					if (rowOcc[i]==2) {
						//std::cout << "2 row found " << i <<"\n" << std::flush;
						j1 = A[i]->getJ();
						_field.init(x1,A[i]->getX());
						j2 = A[i]->next->getJ();
						_field.init(x2,A[i]->next->getX());
						if ((abs(x1)==1/*rowGcd[i]*/) && (abs(x1)== 1/*colGcd[j1]*/)) {
							if ((abs(x2)== 1/*rowGcd[i]*/) && (abs(x2)==1/*colGcd[j2]*/) && (colOcc[j1]> colOcc[j2])) {
								size_t jj = j1;
								Element xx = x1;
								j1 =j2;x1=x2;
								j2 = jj;x2=xx;
							}
							pivotFound = true;
							row2=i;
							ini = i+1; break;
						}
						else
							if ((abs(x2)==1/*rowGcd[i]*/) && (abs(x2)== 1/*colGcd[j2]*/)) {
								size_t jj = j1;
								Element xx = x1;
								j1 =j2;x1=x2;
								j2 = jj;x2=xx;
								pivotFound = true;
								row2=i;
								ini = i+1;break;
							}

					}
				}

				//pivotFound =false;
				if (pivotFound) {
					if (abs(x1)>1) std::cout << "adds " << x1 << "to the diagonal\n" << std::flush;
					mC[j1]=1;
					++rank;
					mR[row2]=1;
					_field.init(x, -x2/x1);
					//std::cout << "found " << j1 << "," << x2 << "," << j2 << "," << x2 << "\n"<< std::flush;
					////std::cout << "reducing column " << j2+1 << " by  (" << j1+1 << "," << x << "}\n" << std::flush;
					//std::cout << "row " << row2 << "\n"<<std::flush;

					GridElement<Element>* p1=AT[j1];
					GridElement<Element>* p2=NULL;
					GridElement<Element>* p2next=AT[j2];

					while (p1 != NULL) {
						while (p2next != NULL) {
							if (p2next->getI() >= p1->getI()) break;
							p2 = p2next;
							p2next = p2next->up;
						}
						Element y; _field.init(y, p1->getX());
						Element z;
						if ((p2next!= NULL) && (p2next->getI()==p1->getI()) ) {
							//std::cout << "updating " << p2next->getI() << " row" << std::flush;
							_field.init (z, x*y)	;
							_field.addin(z, p2next->getX());
							if (z==0) {
								//std::cout << ".....deleting \n" << std::flush;
								p2next = deleteElement(p2next);
							}
							else {
								//std::cout << ".....new value\n" << std::flush;
								p2next->setElement(z);
								p2 = p2next;
								p2next = p2next->up;
							}
						}
						else {
							//std::cout << "adding " << p1->getI() << " row\n" << std::flush;
							_field.init(z, x*y);
							ijElement<Element> X(p1->getI(), j2, z);
							p2=addElement(p2, X);
						}

						p1 = deleteElement(p1);
					}

					if (!Q.empty()) {
						size_t i = Q.front();
						while (rowOcc[i] != 1) {
							Q.pop();
							if (Q.empty()) break;
							i = Q.front();
						}
						if (!Q.empty()) {
							//std::cout << "Adding " << i+1 << " to Q \n" << std::flush;
							reduce(rank, S, mR, mC,os);
							break;
						}
					}
				}
				else {
					std::cout << "Elimination of " << S-1<< " rows at rank " <<rank << "\n" << std::flush;

					size_t i =0;
					// bool tworow=false;
					int init_rank = rank;
					while  (1) {
						while (i < _m) {
							if (rowOcc[i]==0) ++i;
							else if (rowOcc[i]<S)  break;
							else ++i;
						}
						if (i>=_m) break;
						//std::cout << "Elimination of row " << i+1 << std::endl << std::flush;
						//std::cout << "RowGcd=" << rowGcd[i] << "\n" << std::flush;
						GridElement<Element>* r_p = A[i];
						int min_col = _m+1;
						Element xx1;
						int j=-1;
						std::vector<std::pair< size_t, GridElement<Element>*> > j_pts ;
						std::vector<std::pair< Element, GridElement<Element>*> > jnext_pts ;
						while (r_p != NULL) {
							if ((colOcc[r_p->getJ()] < _m+1) && (colOcc[r_p->getJ()] < min_col)) {
								if ((abs(r_p->getX()) == 1/*rowGcd[i]*/) && (abs(r_p->getX())==1/*colGcd[r_p->getJ()]*/)) {
									if (j != -1) {
										j_pts.push_back(std::pair<size_t, GridElement<Element>*>  (j, NULL) );
										jnext_pts.push_back(std::pair<Element, GridElement<Element>*>  (xx1, AT[j]));
									}
									min_col = colOcc[r_p->getJ()];
									_field.init(xx1, r_p->getX());
									j = r_p->getJ();
								}
								else {
									j_pts.push_back( std::pair<size_t, GridElement<Element>*>  (r_p->getJ(), NULL));
									jnext_pts.push_back( std::pair<Element, GridElement<Element>*>  (r_p->getX(), AT[r_p->getJ()]));
								}
							}
							else {
								j_pts.push_back(std::pair<size_t, GridElement<Element>*>  (r_p->getJ(), NULL));
								jnext_pts.push_back(std::pair<Element, GridElement<Element>*>  (r_p->getX(), AT[r_p->getJ()]));
							}

							r_p = r_p->next;
						}

						j_pts.resize(rowOcc[i]-1);
						jnext_pts.resize(rowOcc[i]-1);
						if (j < 0) {
							++i; continue;
						}

						if ((abs(xx1) == 1/*rowGcd[i]*/) && (abs(xx1)==1/*colGcd[j]*/)) {
							if (abs(xx1) > 1) std::cout << "adds " << xx1 << "to the diagonal\n"<<std::flush;
							//std::cout << "Eliminating row "<< i+1 << " by column "<< j+1 << "\n";
							//std::cout << "Element x" << x << "\n" << std::flush;

							mR[i] = 1;
							mC[j] = 1;
							++rank;

							GridElement<Element>* p1=AT[j];

							while (p1 != NULL) {
								typename std::vector<std::pair< size_t, GridElement<Element>*> >::iterator p2 = j_pts.begin();
								typename std::vector<std::pair< Element, GridElement<Element>*> >::iterator p2next = jnext_pts.begin();
								for (; p2 != j_pts.end(); ++p2, ++p2next) {
									//std::cout << "eliminating column " << p2->first << "\n" << std::flush;
									Element x; _field.init(x, -(p2next->first)/xx1);

									while (p2next->second != NULL) {
										if (p2next->second->getI() >= p1->getI()) break;
										p2->second = p2next->second;
										p2next->second = p2next->second->up;
									}
									Element y; _field.init(y, p1->getX());
									Element z;
									if ((p2next->second!= NULL) && (p2next->second->getI()==p1->getI()) ) {
										//std::cout << "updating " << p2next->second->getI() << " row" << std::flush;
										_field.init (z, x*y)        ;
										_field.addin(z, p2next->second->getX());
										if (z==0) {
											//std::cout << ".....deleting \n" << std::flush;
											p2next->second = deleteElement(p2next->second);
										}
										else {
											//std::cout << ".....new value\n" << std::flush;
											p2next->second->setElement(z);
											p2->second = p2next->second;
											p2next->second = p2next->second->up;
										}
									}
									else {
										//std::cout << "adding " << p1->getI() << " row\n" << std::flush;
										_field.init(z, x*y);
										ijElement<Element> X(p1->getI(), p2->first, z);
										p2->second=addElement(p2->second, X);
									}
								}
								p1 = deleteElement(p1);
							}
						}
						++i;
						if (!Q.empty()) {
							size_t k = Q.front();
							while (rowOcc[k] != 1) {
								Q.pop();
								if (Q.empty()) break;
								k = Q.front();
							}
							if (!Q.empty()) {
								//std::cout << "breaking at " << i << " row\n" << std::flush;
								break;
							}
						}
					}

					if ((!Q.empty()) || (rank - init_rank > 100)) {
						reduce(rank, S, mR, mC,os);
					}
					break;
				}
			}
			return rank;
		}

		void write (std::ostream& out)
		{
			size_t Omega =0;
			out << _n << " " << _m << " M\n";
			std::vector<int>::iterator occ_iter = colOcc.begin();
			int j =0;
			for (;occ_iter != colOcc.end(); ++occ_iter,++j) {
				if (*occ_iter > 0) {
					*occ_iter = 0;
					GridElement<Element>* tmp = AT[j];
					GridElement<Element>* tmp2;
					while (tmp != NULL) {
						++Omega;
						out << tmp->getJ()+1 <<" " << tmp->getI()+1 << " " << tmp->getX() << "\n";
						tmp2 = tmp->up;
						delete tmp;
						tmp = tmp2;
					}
				}
			}
			out << "0 0 0\n";
			std::cout << "Omega: " << Omega << "\n" << std::flush;
			rowOcc.clear();
			colOcc.clear();
			//rowGcd.clear();
			//colGcd.clear();
			A.clear();
			AT.clear();

		}

		~Grid()
		{

			std::vector<int>::iterator occ_iter = colOcc.begin();
			int j =0;
			for (;occ_iter != colOcc.end(); ++occ_iter,++j) {
				if (*occ_iter > 0) {
					GridElement<Element>* tmp = AT[j];
					GridElement<Element>* tmp2;
					while (tmp != NULL) {
						tmp2 = tmp->up;
						delete tmp;
						tmp = tmp2;
					}
				}
			}

			rowOcc.clear();
			colOcc.clear();
			//rowGcd.clear();
			//colGcd.clear();
			A.clear();
			AT.clear();

		}
	};

}

#endif //__LINBOX_matrix_grid_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
