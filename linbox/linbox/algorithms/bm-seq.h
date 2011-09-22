/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
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
		typedef Protected::DenseMatrix<Field> value_type;
		typedef typename std::list<value_type>::const_iterator const_iterator;
		typedef int size_type;

	private:

		Field& _F;
		std::list<value_type > _seq;
		size_type _size;
		size_t _row, _col;

	public:
		BM_Seq(Field& F, size_t r, size_t c) : _F(F)
		{
			_row = r;
			_col = c;
			_size = 0;
		}
		BM_Seq(Field& F, size_t r) : _F(F)
		{
			_row = r;
			_col = r;
			_size = 0;
		}
		BM_Seq(int n, value_type& M) : _F(M.field()),  _seq(n, M), _size(n)
		{
			_row = M.rowdim();
			_col = M.coldim();
		}

		BM_Seq() {}

		BM_Seq(const BM_Seq<Field>& S) :
			_F(S._F), _seq(S._seq), _size(S._size), _row(S._row), _col(S._col)
		{}

		BM_Seq & operator=(const BM_Seq<Field>& S)
		{
			if(this != &S){
				(*this)._size = S._size;
				(*this)._row = S._row;
				(*this)._col = S._col;
				(*this)._F = S._F;
				_seq.clear();
				for(typename std::list<value_type>::const_iterator it = S._seq.begin(); it != S._seq.end(); it++)
					_seq.push_back(*it);
			}
			return *this;
		}

		Field& field()
		{
			return _F;
		}

		size_t rowdim()
		{
			return _row;
		}

		size_t coldim()
		{
			return _col;
		}

		const_iterator begin() const
		{
			return _seq.begin();
		}

		const_iterator end() const
		{
			return _seq.end();
		}

		void push_back(const value_type &M)
		{
			if(_row==M.rowdim() && _col==M.coldim()){
				_seq.push_back(M);
				_size++;
			}
		}

		bool operator==(const BM_Seq<Field>& l)
		{
			typename std::list<value_type>::const_iterator it, lit;
			bool test = false;
			if(_size==l._size && _row==l._row && _col==l._col){
				test = true;
				MatrixDomain<Field> MD(_F);
				it = _seq.begin();
				lit = l._seq.begin();
				if(_size==0){
					return test;
				}
				else{
					while(test && it!=_seq.end()){
						test = MD.areEqual(*it,*lit);
						it++;
						lit++;
					}
				}
			}
			return test;
		}

		bool operator!=(const BM_Seq<Field>& l)
		{
			return !(*this == l);
		}

		size_type size()
		{
			return _size;
		}

		class BM_iterator {
		public:
			typedef std::list<typename BM_Seq<Field>::value_type> value_type;

		private:
			typedef typename BM_Seq<Field>::value_type matrix_type;
			Field& _F;
			BM_Seq<Field>& _seq;
			typename BM_Seq<Field>::size_type _size;
			typename BM_Seq<Field>::size_type _t;
			typename BM_Seq<Field>::const_iterator _seqel;
			std::list<matrix_type> _gen;
			std::vector<int> _deg;
			int _delta;
			int _mu;
			int _beta;
			int _sigma;
			int _gensize;
			size_t _row, _col;

		public:
			// This is an enumeration class that tells what state the berlekamp/massey algoithm iterator is in.
			// The four states are:
			// DeltaExceeded = 4
			// SequenceExceeded = 3
			// GeneratorFound = 2
			// GeneratorUnconfirmed = 1
			class TerminationState{
			private:

				int _state;
				friend class BM_iterator;
				TerminationState() : _state(GeneratorUnconfirmed) {}
				TerminationState(int m) : _state(m) {}
			public:
				TerminationState(const TerminationState& t) : _state(t._state) {}
				TerminationState & operator=(const TerminationState & t){
					if(this != &t){
						(*this)._state = t._state;
					}
					return *this;
				}
				bool IsGeneratorUnconfirmed(){
					return _state==GeneratorUnconfirmed;
				}
				bool IsGeneratorFound()
				{
					return _state==GeneratorFound;
				}
				bool IsSequenceExceeded()
				{
					return _state==SequenceExceeded;
				}
				bool IsDeltaExceeded()
				{
					return _state==DeltaExceeded;
				}
			};

		private:
			TerminationState _state;
		public:
			TerminationState state() const
			{
				return _state;
			}
			void setDelta(int d)
			{
				_delta=d;
				if((_delta < 0 || _beta < _delta - _sigma + _mu +1) && _state._state!=3){
					if(_sigma <= _delta || _delta < 0)
						_state._state = GeneratorUnconfirmed;
					else
						_state._state = DeltaExceeded;
				}
				else{
					if(_sigma > _delta)
						_state._state = DeltaExceeded;
					else
						_state._state = GeneratorFound;
				}
			}
			//Constructor
			explicit BM_iterator(BM_Seq<Field>& s, typename BM_Seq<Field>::size_type elinit=0) :
				_F(s.field()), _seq(s)
			{
				_row = s.rowdim();
				_col = s.coldim();
				_size = _seq.size();
				_t = elinit;
				_delta = -1;
				_seqel = _seq.begin();
				_deg = std::vector<int>(_row+_col);
				for(size_t i = _col; i < _row+_col; i++)
					_deg[i] = 1;
				typename Field::Element one;
				_F.init(one,1);
				matrix_type gen1(_F,_col,_row+_col);
				for(size_t i = 0; i<_col; i++)
					gen1.setEntry(i,i,one);
				_gen.push_back(gen1);
				_gensize = 1;
				if(_size==0 || _t==_size)
					_state._state = SequenceExceeded;
				_sigma = 0;
				_mu = 0;
				_beta = 1;
			}

			//Copy constructor
			BM_iterator(const BM_Seq<Field>::BM_iterator & it) :
				_F(it._F), _seq(it._seq), _size(it._size), _t(it._t),
				_seqel(it._seqel), _gen(it._gen), _deg(it._deg),
				_delta(it._delta), _mu(it._mu), _beta(it._beta),
				_sigma(it._sigma), _gensize(it._gensize),
				_row(it._row), _col(it._col), _state(it._state) {}

			//Assignment operator not overloaded since Protected::DenseMatrix class has overloaded assignment error
			//Overloaded assignment operator
			BM_iterator& operator=(const typename BM_Seq<Field>::BM_iterator& it)
			{
				if(this != &it){
					(*this)._F       = it._F;
					(*this)._row     = it._row;
					(*this)._col     = it._col;
					(*this)._seq     = it._seq;
					(*this)._size    = it._size;
					(*this)._t       = it._t;
					(*this)._seqel   = it._seqel;
					(*this)._deg     = it._deg;
					(*this)._gensize = it._gensize;
					(*this)._delta   = it._delta;
					(*this)._mu      = it._mu;
					(*this)._sigma   = it._sigma;
					(*this)._beta    = it._beta;
					(*this)._state   = it._state;
					_gen.clear();
					for(typename std::list<matrix_type>::const_iterator git = it._gen.begin(); git != it._gen.end(); git++)
						_seq.push_back(*git);
				}
				return (*this);
			}

			bool operator==(const BM_Seq<Field>::BM_iterator& it)
			{
				TerminationState check = it.state();
				bool test1 = (_seq==it._seq);
				bool test2 = (_t==it._t);
				bool test3 = _delta==it._delta;
				bool test4 = (_state._state == check._state && _state.IsSequenceExceeded());
				return (test1  && test2  && (test3 || test4));
			}

			bool operator!=(const BM_iterator& it)
			{
				return !((*this) == it);
			}
		private:

			// Column Copy
			template <class Matrix>
			void ColumnCopy(Matrix &M, const Matrix &A, size_t i)
			{
				size_t rowd = A.rowdim();
				for(size_t j = 0; j<rowd; j++){
					M.setEntry(j,i,A.getEntry(j,i));
				}
			}
			// Column Swap
			template <class Matrix>
			void ColumnSwap(Matrix &M, size_t i, size_t j)
			{
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
			void ColumnAdd(Matrix &M, size_t i, size_t j, typename Matrix::Element el)
			{
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
			Matrix  Algorithm3dot2(Matrix &D, std::vector<int> &d, int &mu, int &sigma, int &beta)
			{
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
							}
							else
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
					}
					else{
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
							ColumnSwap(D, piv, m+i);
						}
						else{
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
			BM_iterator& operator++()
			{
				//See if a matrix has been pushed on the sequence
				//if it has, then recompute the seqel since it may
				//have become corrupt.
				//Also reset the size to the correct size of the sequence
				if(_size < _seq.size()){
					_seqel = _seq.begin();
					for(int i = 0; i<_t; i++)
						_seqel++;
					_size = _seq.size();
				}
				//if the iterator points past the seq elements, do nothing
				if(_t == _size){
					return *this;
				}
				//Initialize the discrepancy
				matrix_type disc(_F,_row, _row+_col);
				//Create two iterators, one for seq, and one for gen
				typename BM_Seq<Field>::const_iterator cseqit;
				typename std::list<matrix_type>::iterator genit;
				//get a iterator to the seq element to be processed
				cseqit = _seqel;
				//Create a matrix domain for addition and multiplication
				MatrixDomain<Field> MD(_F);
				//Compute the discrepancy
				for(genit = _gen.begin(); genit!=_gen.end(); genit++){
					MD.axpyin(disc,*cseqit,*genit);
					cseqit--;
				}
				//Compute tau with Algorith3.2
				matrix_type tau(Algorithm3dot2(disc, _deg, _mu, _sigma, _beta));
				//Multiply tau into each matrix in the generator
				for(genit = _gen.begin(); genit!=_gen.end(); genit++){
					MD.mulin(*genit,tau);
				}
				//Increment the auxiliary degrees and beta
				for(size_t j = _col; j <_row+_col; j++)
					_deg[j]++;
				_beta++;
				//Add a zero matrix to the end of the generator if needed.
				int tmax = _deg[0];
				for(size_t j = 1; j<_row+_col; j++)
					if(tmax < _deg[j])
						tmax = _deg[j];
				if(tmax+1 > _gensize){
					_gen.push_back(matrix_type(_F,_col,_row+_col));
					_gensize++;
				}
				//Mimic multiplication be z in the auxiliary columns
				typename std::list<matrix_type>::reverse_iterator g1,g2;
				g1 = _gen.rbegin();
				g2 = _gen.rbegin();
				g1++;
				while(g1!=_gen.rend()){
					for(size_t k = _col; k < _row+_col; k++){
						ColumnCopy(*g2,*g1,k);
					}
					g1++;
					g2++;
				}
				genit = _gen.begin();
				matrix_type z1(_F,_col,_row+_col);
				for(size_t k = _col; k < _row+_col; k++)
					ColumnCopy(*genit, z1,k);
				//Increment the t and seqel to the next element
				_t++;
				_seqel++;
				//Update the state
				if(_delta < 0 || _beta < _delta - _sigma + _mu +1){
					if(_t == _size)
						_state._state = SequenceExceeded;
					else{
						if(_sigma > _delta && _delta >= 0)
							_state._state = DeltaExceeded;
						else
							_state._state = GeneratorUnconfirmed;
					}
				}
				else{
					if(_sigma > _delta)
						_state._state = DeltaExceeded;
					else
						_state._state = GeneratorFound;
				}

				return *this;
			}
			BM_iterator operator++(int)
			{
				//Create a copy of this
				BM_iterator temp(*this);

				//See if a matrix has been pushed on the sequence
				//if it has, then recompute the seqel since it may
				//have become corrupt.
				//Also reset the size to the correct size of the sequence
				if(_size < _seq.size()){
					_seqel = _seq.begin();
					for(int i = 0; i<_t; i++)
						_seqel++;
					_size = _seq.size();
				}
				//if the iterator points past the seq elements, do nothing
				if(_t == _size){
					return *this;
				}
				//Initialize the discrepancy
				matrix_type disc(_F,_row, _row+_col);
				//Create two iterators, one for seq, and one for gen
				typename BM_Seq<Field>::const_iterator cseqit;
				typename std::list<matrix_type>::iterator genit;
				//get an iterator to the seq element to be processed
				cseqit = _seqel;
				//Create a matrix domain for addition and multiplication
				MatrixDomain<Field> MD(_F);
				//Compute the discrepancy
				for(genit = _gen.begin(); genit!=_gen.end(); genit++, cseqit--){
					MD.axpyin(disc,*cseqit,*genit);
				} // cost: k*n^3 (nxn matrix muladds where k is current generator length)
				  // is a reductive addition over independent muls.
				//Compute tau with Algorith3.2
				matrix_type tau(Algorithm3dot2(disc, _deg, _mu, _sigma, _beta));
				  // cost: n^3 for elim on n x about 2n
				//Multiply tau into each matrix in the generator
				for(genit = _gen.begin(); genit!=_gen.end(); genit++){
					MD.mulin(*genit,tau);
				} // cost: k*n^3 (nxn matrix muls where k is current generator length)
				  // is k independent muls with a shared mat tau.
				//Increment the auxiliary degrees and beta
				for(size_t j = _col; j <_row+_col; j++)
					_deg[j]++;
				_beta++;
				//Add a zero matrix to the end of the generator if needed.
				int tmax = _deg[0];
				for(size_t j = 1; j<_row+_col; j++)
					if(tmax < _deg[j])
						tmax = _deg[j];
				if(tmax+1 > _gensize){
					_gen.push_back(matrix_type(_F,_col,_row+_col));
					_gensize++;
				}
				//Mimic multiplication by z in the auxiliary columns
				typename std::list<matrix_type>::reverse_iterator g1,g2;
				g1 = _gen.rbegin();
				g2 = _gen.rbegin();
				g1++;
				while(g1!=_gen.rend()){
					for(size_t k = _col; k < _row+_col; k++){
						ColumnCopy(*g2,*g1,k);
					}
					g1++;
					g2++;
				}
				genit = _gen.begin();
				matrix_type z1(_F,_col,_row+_col);
				for(size_t k = _col; k < _row+_col; k++)
					ColumnCopy(*genit, z1,k);
				//Increment the t and seqel to the next element
				_t++;
				_seqel++;
				//Update the state
				if(_delta < 0 || _beta < _delta - _sigma + _mu +1){
					if(_t == _size)
						_state._state = SequenceExceeded;
					else{
						if(_sigma > _delta && _delta >= 0)
							_state._state = DeltaExceeded;
						else
							_state._state = GeneratorUnconfirmed;
					}
				}
				else{
					if(_sigma > _delta)
						_state._state = DeltaExceeded;
					else
						_state._state = GeneratorFound;
				}

				return temp;
			}
			//return a reference to the current generator, in its algorithmic reversed form
			value_type& operator*()
			{
				return _gen;
			}
			//overload the pointer operator
			value_type* operator->()
			{
				return &_gen;
			}
			//Return a vector representing the reversal, by nominal degree, of the current generator
			std::vector<matrix_type> GetGenerator()
			{
				std::vector<matrix_type> revgen(_mu+1, matrix_type(_F,_col,_col));
				for(size_t i = 0; i<_col; i++){
					typename std::list<matrix_type>::iterator genit = _gen.begin();
					for(int j = 0; j < _deg[i]+1; j++){
						ColumnCopy(revgen[_deg[i]-j], *genit,i);
						genit++;
					}
				}
				return revgen;
			}

			typename BM_Seq<Field>::size_type get_t()
			{
				return _t;
			}

			int get_mu()
			{
				return _mu;
			}

			int get_sigma()
			{
				return _sigma;
			}
			int get_beta()
			{
				return _beta;
			}
			int get_delta()
			{
				return _delta;
			}
		};
		//return an initialized BM_iterator
		typename BM_Seq<Field>::BM_iterator BM_begin()
		{
			return typename BM_Seq<Field>::BM_iterator(*this);
		}
		//return an initialized BM_iterator that points to one past the end of the sequence
		typename BM_Seq<Field>::BM_iterator BM_end()
		{
			return typename BM_Seq<Field>::BM_iterator(*this, _size);
		}
		/**/
	};
}

#endif
