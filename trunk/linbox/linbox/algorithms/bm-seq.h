/* -*- mode: C++; tab-width: 2; indent-tabs-mode: t; c-basic-offset: 2 -*- */

/* linbox/algorithms/bm-seq.h
 * Copyright (C) 2008 George Yuhasz
 *
 * Written by George Yuhasz gyuhasz@math.ncsu.edu
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



#ifndef __BM_SEQ_H
#define __BM_SEQ_H

//Preprocessor variables for the state of BM_iterators
#define DeltaExceeded  4
#define SequenceExceeded  3
#define GeneratorFound  2
#define GeneratorUnconfirmed  1 

#include <vector>
#include <list>
#include <set>

#include <linbox/util/timer.h>
#include <linbox/blackbox/dense.h>
#include <linbox/matrix/matrix-domain.h>



namespace LinBox {
	template<class _Field>
	class BM_Seq {
		
	public:
				
		typedef _Field Field;
		typedef DenseMatrix<Field> value_type;
		typedef typename std::list<value_type>::const_iterator const_iterator;
		typedef int size_type;
		
	private:

		Field& F;
		std::list<value_type > seq;
		size_type _size;
		size_t row, col;

	public:
		BM_Seq(Field& _F, size_t r, size_t c) : F(_F) {
			row = r;
			col = c;
			_size = 0;
		}
		BM_Seq(Field& _F, size_t r) : F(_F){
			row = r;
			col = r;
			_size = 0;
		}
		BM_Seq(int n, value_type& M) : F(M.field()),  seq(n, M), _size(n)
		{ row = M.rowdim(); col = M.coldim();}

		BM_Seq() {}

		BM_Seq(const BM_Seq<Field>& S) : F(S.F), seq(S.seq), _size(S._size), row(S.row), col(S.col){}

		BM_Seq & operator=(const BM_Seq<Field>& S){
			if(this != &S){
				(*this)._size = S._size;
				(*this).row = S.row;
				(*this).col = S.col;
				(*this).F = S.F;
				seq.clear();
				for(typename std::list<value_type>::const_iterator it = S.seq.begin(); it != S.seq.end(); it++)
								seq.push_back(*it);
				return *this;
			}
		}

		Field& field(){ return F; }

		size_t rowdim(){ return row; }

		size_t coldim(){ return col; }

		const_iterator begin() const { return seq.begin(); }

		const_iterator end() const { return seq.end(); }

		void push_back(const value_type &M){ 
			if(row==M.rowdim() && col==M.coldim()){
				seq.push_back(M);
				_size++;
			}
		}
		bool operator==(const BM_Seq<Field>& l){
			typename std::list<value_type>::const_iterator it, lit;
			bool test = false;
			if(_size==l._size && row==l.row && col==l.col){
				test = true;
				MatrixDomain<Field> MD(F);
				it = seq.begin();
				lit = l.seq.begin();
				if(_size==0){
					return test;
				}else{
					while(test && it!=seq.end()){
						test = MD.areEqual(*it,*lit);
						it++;
						lit++;
					}
				}	
			}
			return test;
		}
		bool operator!=(const BM_Seq<Field>& l){
			return !(*this == l);
		}

		size_type size(){ return _size; }

		class BM_iterator {
			public: 
				typedef std::list<typename BM_Seq<Field>::value_type> value_type;

			private:
				typedef typename BM_Seq<Field>::value_type matrix_type;
				Field& F;
				BM_Seq<Field>& seq;
				typename BM_Seq<Field>::size_type size;
				typename BM_Seq<Field>::size_type t;
				typename BM_Seq<Field>::const_iterator seqel;
				std::list<matrix_type> gen;
				std::vector<int> deg;
				int delta;
				int mu;
				int beta;
				int sigma;
				int gensize;
				size_t row, col;

			public:
				// This is an enumeration class that tells what state the berlekamp/massey algoithm iterator is in.
				// The four states are:
				// DeltaExceeded = 4
				// SequenceExceeded = 3
				// GeneratorFound = 2
				// GeneratorUnconfirmed = 1 
				class TerminationState{
					private:
						
						int state;
						friend class BM_iterator;
						TerminationState() : state(GeneratorUnconfirmed) {}
						TerminationState(int m) : state(m) {}
					public:
						TerminationState(const TerminationState& t) : state(t.state) {}
						TerminationState & operator=(const TerminationState & t){
							if(this != &t){
								(*this).state = t.state;
								return *this;
							}
						}
						bool IsGeneratorUnconfirmed() { return state==GeneratorUnconfirmed; }
						bool IsGeneratorFound() {return state==GeneratorFound; }
						bool IsSequenceExceeded() { return state==SequenceExceeded; }
						bool IsDeltaExceeded() { return state==DeltaExceeded; }
				};

			private:
				TerminationState _state;
			public:
				TerminationState state() const { return _state; }
				void setDelta(int d) { 
					delta=d;
					if((delta < 0 || beta < delta - sigma + mu +1) && _state.state!=3){
						if(sigma <= delta || delta < 0)
							_state.state = GeneratorUnconfirmed;
						else
							_state.state = DeltaExceeded;
					}else{
						if(sigma > delta)
							_state.state = DeltaExceeded;
						else
							_state.state = GeneratorFound;
					}
			       	}
				//Constructor
				explicit BM_iterator(BM_Seq<Field>& s, typename BM_Seq<Field>::size_type elinit=0) : F(s.field()), seq(s) {
					row = s.rowdim();
					col = s.coldim();
					size = seq.size();
					t = elinit;
					delta = -1;
					seqel = seq.begin();
					deg = std::vector<int>(row+col);
					for(size_t i = col; i < row+col; i++)
						deg[i] = 1;
					typename Field::Element one;
					F.init(one,1);
					matrix_type gen1(F,col,row+col);
					for(size_t i = 0; i<col; i++)
						gen1.setEntry(i,i,one);
					gen.push_back(gen1);
					gensize = 1;
					if(size==0 || t==size)
						_state.state = SequenceExceeded;
					sigma = 0;
					mu = 0;
					beta = 1;
				}
				
				//Copy constructor
				BM_iterator(const BM_Seq<Field>::BM_iterator & it) : F(it.F), seq(it.seq), size(it.size), t(it.t), seqel(it.seqel), gen(it.gen), deg(it.deg),  delta(it.delta), mu(it.mu), beta(it.beta), sigma(it.sigma), gensize(it.gensize), row(it.row), col(it.col), _state(it._state) {}

				//Assignment operator not overloaded since DenseMatrix class has overloaded assignment error
				//Overloaded assignment operator
				BM_iterator& operator=(const typename BM_Seq<Field>::BM_iterator& it){
					if(this != &it){
						(*this).F = it.F;
						(*this).row = it.row;
						(*this).col = it.col;
						(*this).seq = it.seq;
						(*this).size = it.size;
						(*this).t = it.t;
						(*this).seqel = it.seqel;
						(*this).deg = it.deg;
						(*this).gensize = it.gensize;
						(*this).delta = it.delta;
						(*this).mu = it.mu;
						(*this).sigma = it.sigma;
						(*this).beta = it.beta;
						(*this)._state = it._state;
						std::cout << "clearing generator" << std::endl;
						gen.clear();
						std::cout << "generator cleared" << std::endl;
						for(typename std::list<matrix_type>::const_iterator git = it.gen.begin(); git != it.gen.end(); git++)
							seq.push_back(*git);
						std::cout << "generator copied" << std::endl;
						return (*this);
					}
				}

				bool operator==(const BM_Seq<Field>::BM_iterator& it){
								TerminationState check = it.state();
								bool test1 = (seq==it.seq);
								bool test2 = (t==it.t);
								bool test3 = delta==it.delta;
								bool test4 = (_state.state == check.state && _state.IsSequenceExceeded());
					return (test1  && test2  && (test3 || test4));
				}
				bool operator!=(const BM_iterator& it){
					return !((*this) == it);
				}
			private:

				// Column Copy 
				template <class Matrix>
				void ColumnCopy(Matrix &M, const Matrix &A, size_t i){
					size_t rowd = A.rowdim();
					for(size_t j = 0; j<rowd; j++){
						M.setEntry(j,i,A.getEntry(j,i));
					}
				}
				// Column Swap 
				template <class Matrix>
				void ColumnSwap(Matrix &M, size_t i, size_t j){
					typename Matrix::Field F = M.field();
					typename Matrix::Element t;
					F.init(t,0);
					size_t rowd = M.rowdim();
					for(size_t k = 0; k < rowd; k++){
						F.assign(t,M.getEntry(k,i));
						M.setEntry(k,i,M.getEntry(k,j));
						M.setEntry(k,j,t);
					}
				}
				// Column Operation 
				template <class Matrix>
				void ColumnAdd(Matrix &M, size_t i, size_t j, typename Matrix::Element el){
					typename Matrix::Field F = M.field();
					typename Matrix::Element t;
					F.init(t,0);
					size_t rowd = M.rowdim();
					for (size_t k=0; k<rowd; k++){
						F.mul(t, M.getEntry(k,j), el);
						F.addin(t,M.getEntry(k,i));
						M.setEntry(k,i,t);
					}
				}
				template <class Matrix>
				Matrix  Algorithm3dot2(Matrix &D, std::vector<int> &d, int &mu, int &sigma, int &beta){
					typename Matrix::Field F = D.field();
					typename Matrix::Element one, pivel;
					F.init(one, 1);
					F.init(pivel,0);
					// Retrieve the row and column dimensions of the sequence and the dimension of the discrepancy
					size_t n = D.rowdim();
					size_t nm  = D.coldim();
					size_t m = nm-n;
					//Initialize tau to the identity matrix
					Matrix tau(F,nm,nm);
					for(size_t i = 0; i<nm; i++)
						tau.setEntry(i,i,one);
					//Create the set of generator columns
					std::set<size_t> gen;
					typedef std::set<size_t>::key_type index_type;
					for(index_type i=0; i<m; i++)
						gen.insert(i);
					for(index_type i = 0; i<n; i++){
							//Compute pi, the columns of D with nonzero entries in row i
							std::set<size_t> pi;
							pi.insert(m+i);
						for(typename std::set<size_t>::iterator genit = gen.begin(); genit != gen.end(); genit++){
							if(!F.isZero(D.getEntry(i,*genit)))
								pi.insert(*genit);
						}

						//Choose the pivot row with the smallest nominal degree
						index_type piv = m+i;
						for(std::set<size_t>::iterator itpi = pi.begin(); itpi != pi.end(); itpi++){
							size_t j = *itpi;
							if(d[j] <= d[piv]){
								if(d[j]==d[piv]){
									if(piv < m+i){
										if(j<piv)
											piv = j;
									}
								}else
									piv = j;
							}							
						}
						pi.erase(piv);
						F.assign(pivel,D.getEntry(i,piv));
						//Handle the case when piv=m+i, so no swap is done
						if(piv==m+i){
							for(std::set<size_t>::iterator itpi = pi.begin(); itpi != pi.end(); itpi++){
								typename Matrix::Element temp;
								F.init(temp,D.getEntry(i, *itpi));
								F.negin(temp);
								F.divin(temp,pivel);
								ColumnAdd(tau, *itpi, piv, temp);
								ColumnAdd(D, *itpi, piv, temp);
							}
						}else{
							//Remove column index m+i and handle it separately
							pi.erase(m+i);
							//Eliminate nonzero discrepancies in generator columns
							for(typename std::set<size_t>::iterator itpi = pi.begin(); itpi != pi.end(); itpi++){
								typename Matrix::Element temp;
								F.init(temp,D.getEntry(i, *itpi));
								F.negin(temp);
								F.divin(temp,pivel);
								ColumnAdd(tau, *itpi, piv, temp);
								ColumnAdd(D, *itpi, piv, temp);
							}
							typename Matrix::Element auxel;
							F.init(auxel,D.getEntry(i,m+i));
							//Perform a major change and update an initialized auxiliary column
							if(!F.isZero(auxel)){
								typename Matrix::Element temp;
								F.init(temp,D.getEntry(i, m+i));
								F.negin(temp);
								F.divin(temp,pivel);
								ColumnAdd(tau, m+i, piv, temp);
								ColumnAdd(D, m+i, piv, temp);
								ColumnSwap(tau,piv, m+i);
								ColumnSwap(D, piv, m+i);	}else{
								ColumnAdd(tau,m+i,piv,one);
								ColumnAdd(D,m+i,piv,one);
								gen.erase(piv);
							}
							int tempdeg = d[piv];
							d[piv] = d[m+i];
							d[m+i] = tempdeg;
							if(tempdeg < beta)
								beta = tempdeg;
							if(d[piv] > mu)
								mu = d[piv];
							sigma = sigma - tempdeg + d[piv];
						}
					}
					return tau;
				}
			public:
				BM_iterator& operator++(){
					//See if a matrix has been pushed on the sequence
					//if it has, then recompute the seqel since it may
					//have become corrupt.
					//Also reset the size to the correct size of the sequence
					if(size < seq.size()){
						seqel = seq.begin();
						for(int i = 0; i<t; i++)
							seqel++;
						size = seq.size();
					}
					//if the iterator points past the seq elements, do nothing
					if(t == size){
						return *this;
					}
					//Initialize the discrepancy
					matrix_type disc(F,row, row+col);
					//Create two iterators, one for seq, and one for gen
					typename BM_Seq<Field>::const_iterator cseqit;
					typename std::list<matrix_type>::iterator genit;
					//get a iterator to the seq element to be processed
					cseqit = seqel;
					//Create a matrix domain for addition and multiplication
					MatrixDomain<Field> MD(F);
					//Compute the discrepancy
					for(genit = gen.begin(); genit!=gen.end(); genit++){
						MD.axpyin(disc,*cseqit,*genit);
						cseqit--;
					}
					//Compute tau with Algorith3.2
					matrix_type tau(Algorithm3dot2(disc, deg, mu, sigma, beta));
					//Multiply tau into each matrix in the generator
					for(genit = gen.begin(); genit!=gen.end(); genit++){
						MD.mulin(*genit,tau);
					}
					//Increment the auxiliary degrees and beta
					for(size_t j = col; j <row+col; j++)
						deg[j]++;
					beta++;
					//Add a zero matrix to the end of the generator if needed.
					int tmax = deg[0];
					for(size_t j = 1; j<row+col; j++)
						if(tmax < deg[j])
							tmax = deg[j];
					if(tmax+1 > gensize){
						gen.push_back(matrix_type(F,col,row+col));
						gensize++;
					}
					//Mimic multiplication be z in the auxiliary columns
					typename std::list<matrix_type>::reverse_iterator g1,g2;
					g1 = gen.rbegin();
					g2 = gen.rbegin();
					g1++;
					while(g1!=gen.rend()){
						for(size_t k = col; k < row+col; k++){
							ColumnCopy(*g2,*g1,k);
						}
						g1++;
						g2++;
					}
					genit = gen.begin();
					matrix_type z1(F,col,row+col);
					for(size_t k = col; k < row+col; k++)
						ColumnCopy(*genit, z1,k);
					//Increment the t and seqel to the next element
					t++;
					seqel++;
					//Update the state
					if(delta < 0 || beta < delta - sigma + mu +1){
						if(t == size)
							_state.state = SequenceExceeded;
						else{
							if(sigma > delta && delta >= 0)
								_state.state = DeltaExceeded;
							else
								_state.state = GeneratorUnconfirmed;
						}
					}else{
						if(sigma > delta)
							_state.state = DeltaExceeded;
						else
							_state.state = GeneratorFound;
					}

					return *this;
				}
				BM_iterator operator++(int){
					//Create a copy of this
					BM_iterator temp(*this);
					
					//See if a matrix has been pushed on the sequence
					//if it has, then recompute the seqit since it may
					//have become corrupt.
					//Also reset the size to the correct size of the sequence
					if(size < seq.size()){
						seqel = seq.begin();
						for(int i = 0; i<t; i++)
							seqel++;
						size = seq.size();
					}
					//if the iterator points past the seq elements, do nothing
					if(t == size){
						return *this;
					}
					//Initialize the discrepancy
					matrix_type disc(F,row, row+col);
					//Create two iterators, one for seq, and one for gen
					typename BM_Seq<Field>::const_iterator cseqit;
					typename std::list<matrix_type>::iterator genit;
					//get an iterator to the seq element to be processed
					cseqit = seqel;
					//Create a matrix domain for addition and multiplication
					MatrixDomain<Field> MD(F);
					//Compute the discrepancy
					for(genit = gen.begin(); genit!=gen.end(); genit++){
						MD.axpyin(disc,*cseqit,*genit);
						cseqit--;
					}
					//Compute tau with Algorith3.2
					matrix_type tau(Algorithm3dot2(disc, deg, mu, sigma, beta));
					//Multiply tau into each matrix in the generator
					for(genit = gen.begin(); genit!=gen.end(); genit++){
						MD.mulin(*genit,tau);
					}
					//Increment the auxiliary degrees and beta
					for(size_t j = col; j <row+col; j++)
						deg[j]++;
					beta++;
					//Add a zero matrix to the end of the generator if needed.
					int tmax = deg[0];
					for(size_t j = 1; j<row+col; j++)
						if(tmax < deg[j])
							tmax = deg[j];
					if(tmax+1 > gensize){
						gen.push_back(matrix_type(F,col,row+col));
						gensize++;
					}
					//Mimic multiplication be z in the auxiliary columns
					typename std::list<matrix_type>::reverse_iterator g1,g2;
					g1 = gen.rbegin();
					g2 = gen.rbegin();
					g1++;
					while(g1!=gen.rend()){
						for(size_t k = col; k < row+col; k++){
							ColumnCopy(*g2,*g1,k);
						}
						g1++;
						g2++;
					}
					genit = gen.begin();
					matrix_type z1(F,col,row+col);
					for(size_t k = col; k < row+col; k++)
						ColumnCopy(*genit, z1,k);
					//Increment the t and seqel to the next element
					t++;
					seqel++;
					//Update the state
					if(delta < 0 || beta < delta - sigma + mu +1){
						if(t == size)
							_state.state = SequenceExceeded;
						else{
							if(sigma > delta && delta >= 0)
								_state.state = DeltaExceeded;
							else
								_state.state = GeneratorUnconfirmed;
						}
					}else{
						if(sigma > delta)
							_state.state = DeltaExceeded;
						else
							_state.state = GeneratorFound;
					}

					return temp;
				}
				//return a reference to the current generator, in its algorithmic reversed form
				value_type& operator*(){
					return gen;
				}
				//overload the pointer operator
				value_type* operator->(){
					return &gen;
				}
				//Return a vector representing the reversal, by nominal degree, of the current generator
				std::vector<matrix_type> GetGenerator(){
								std::vector<matrix_type> revgen(mu+1, matrix_type(F,col,col));
					for(size_t i = 0; i<col; i++){
						typename std::list<matrix_type>::iterator genit = gen.begin();
						for(int j = 0; j < deg[i]+1; j++){
							ColumnCopy(revgen[deg[i]-j], *genit,i);
						       	genit++;
						}
					}
					return revgen;
				}

				typename BM_Seq<Field>::size_type get_t() { return t; }
				int get_mu() { return mu; }
				int get_sigma() { return sigma; }
				int get_beta() { return beta; }
				int get_delta() { return delta; }
		};
		//return an initialized BM_iterator
		typename BM_Seq<Field>::BM_iterator BM_begin(){
			return typename BM_Seq<Field>::BM_iterator(*this);
		}
		//return an initialized BM_iterator that points to one past the end of the sequence
		typename BM_Seq<Field>::BM_iterator BM_end(){
			return typename BM_Seq<Field>::BM_iterator(*this, _size);
		}
		/**/
	};
}

#endif
