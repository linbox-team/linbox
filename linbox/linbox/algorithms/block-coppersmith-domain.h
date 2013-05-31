/* linbox/algorithms/block-coppersmith-domain.h
 * Copyright (C) 2012 George Yuhasz
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301  USA
 * ========LICENCE========
 */



#ifndef __LINBOX_coppersmith_block_domain_H
#define __LINBOX_coppersmith_block_domain_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include "linbox/util/commentator.h"
#include "linbox/field/unparametric.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/blas-matrix.h"

#define DEFAULT_BLOCK_EARLY_TERM_THRESHOLD 10
//Preprocessor variables for the state of BM_iterators
#define DeltaExceeded  4
#define SequenceExceeded  3
#define GeneratorFound  2
#define GeneratorUnconfirmed  1

namespace LinBox
{

    /** Compute the linear generator of a sequence of matrices.
     *
     * This class encapsulates the functionality required for computing
     * the block minimal polynomial of a matrix.
     * @bib
     * Yuhasz thesis ...
     */
    template<class _Field, class _Sequence>
    class BlockCoppersmithDomain {

    public:
        typedef _Field                           Field;
        typedef typename Field::Element        Element;
        typedef _Sequence                     Sequence;
        typedef typename MatrixDomain<Field>::Matrix Coefficient;


    protected:
        Sequence                          *_container;
        const Field                           *_field;
        MatrixDomain<Field>                       _MD;
        unsigned long            EARLY_TERM_THRESHOLD;


    public:

        BlockCoppersmithDomain (const BlockCoppersmithDomain<Field,
Sequence> &Mat, unsigned long ett_default =
DEFAULT_BLOCK_EARLY_TERM_THRESHOLD) :
            _container(Mat._container), _field(Mat._field),
            _MD(Mat.field()),  EARLY_TERM_THRESHOLD (ett_default)
        {}
        BlockCoppersmithDomain (Sequence *D, unsigned long ett_default
= DEFAULT_BLOCK_EARLY_TERM_THRESHOLD) :
            _container(D), _field(&(D->field ())), _MD(D->field ()),
EARLY_TERM_THRESHOLD (ett_default)
        {}


        // field of the domain
        const Field &field    () const { return *_field; }
        const Field &getField    () const { return *_field; } // deprecated

        // sequence of the domain
        Sequence *getSequence () const
        { return _container; }

        // the principal function
        std::vector<size_t>  right_minpoly (std::vector<Coefficient> &P);

        // left minimal generating polynomial of the sequence
        // This _MAY_ get defined eventually.
        std::vector<size_t> & left_minpoly (std::vector<Coefficient> &P);
        // { /* transpose seq and transpose result*/ }

        // For compatibility with massey-domain (element sequence case).
        std::vector<size_t> & operator()(std::vector<Coefficient> &P)
        { return right_minpoly(P); }


    private:

        // bm-seq.h stuff can go here.
	class BM_Seq {

	public:

		typedef _Field Field;
		typedef MatrixDomain<Field> Domain;
		typedef typename MatrixDomain<Field>::Matrix value_type;
		typedef typename std::list<value_type>::const_iterator const_iterator;
		typedef int size_type;

	private:

		const Field* _field;
		Domain  _MD;
		std::list<value_type > _seq;
		size_type _size;
		size_t _row, _col;

	public:
		BM_Seq(const Field& F, size_t r, size_t c) : _field(&F), _MD(F)
		{
			_row = r;
			_col = c;
			_size = 0;
		}
		BM_Seq(const Field& F, size_t r) : _field(&F), _MD(F)
		{
			_row = r;
			_col = r;
			_size = 0;
		}
		BM_Seq(int n, value_type& M) : _field(&(M.field())),  _MD(M.field()), _seq(n, M), _size(n)
		{
			_row = M.rowdim();
			_col = M.coldim();
		}

		BM_Seq() {}

		BM_Seq(const BM_Seq& S) :
			_field(S._field), _MD(S._MD), _seq(S._seq), _size(S._size), _row(S._row), _col(S._col)
		{}

		BM_Seq & operator=(const BM_Seq& S)
		{
			if(this != &S){
				(*this)._size = S._size;
				(*this)._row = S._row;
				(*this)._col = S._col;
				(*this)._field = S._field;
				(*this)._MD = S._MD;
				_seq.clear();
				for(typename std::list<value_type>::const_iterator it = S._seq.begin(); it != S._seq.end(); it++)
					_seq.push_back(*it);
			}
			return *this;
		}

		const Field& field()
		{
			return *_field;
		}

		Domain& domain()
		{
			return _MD;
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

		bool operator==(const BM_Seq& l)
		{
			typename std::list<value_type>::const_iterator it, lit;
			bool test = false;
			if(_size==l._size && _row==l._row && _col==l._col){
				test = true;
				it = _seq.begin();
				lit = l._seq.begin();
				if(_size==0){
					return test;
				}
				else{
					while(test && it!=_seq.end()){
						test = _MD.areEqual(*it,*lit);
						it++;
						lit++;
					}
				}
			}
			return test;
		}

		bool operator!=(const BM_Seq& l)
		{
			return !(*this == l);
		}

		size_type size()
		{
			return _size;
		}

		class BM_iterator {
		public:
			typedef std::list<typename BM_Seq::value_type> value_type;

		private:
			typedef typename BM_Seq::value_type matrix_type;
			const Field* _field;
			BM_Seq& _seq;
			typename BM_Seq::size_type _size;
			typename BM_Seq::size_type _t;
			typename BM_Seq::const_iterator _seqel;
			std::list<matrix_type> _gen;
			std::vector<size_t> _deg;
			size_t _delta;
			size_t _mu;
			size_t _beta;
			size_t _sigma;
			size_t _gensize;
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
			}; // TerminationState

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
			explicit BM_iterator(BM_Seq& s, typename BM_Seq::size_type elinit=0) :
				_field(&(s.field())), _seq(s)
			{
				_row = s.rowdim();
				_col = s.coldim();
				_size = _seq.size();
				_t = elinit;
				_delta = (size_t)-1; // BB: is it meant ?
				_seqel = _seq.begin();
				_deg = std::vector<size_t>(_row+_col);
				for(size_t i = _col; i < _row+_col; i++)
					_deg[i] = 1;
				typename Field::Element one;
				field().init(one,1);
				matrix_type gen1(field(),_col,_row+_col);
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
			BM_iterator(const BM_Seq::BM_iterator & it) :
				_field(it._field), _seq(it._seq), _size(it._size), _t(it._t),
				_seqel(it._seqel), _gen(it._gen), _deg(it._deg),
				_delta(it._delta), _mu(it._mu), _beta(it._beta),
				_sigma(it._sigma), _gensize(it._gensize),
				_row(it._row), _col(it._col), _state(it._state) {}

			inline const Field & field() const { return *_field; }

			//Assignment operator not overloaded since BlasMatrix class has overloaded assignment error
			//Overloaded assignment operator
			BM_iterator& operator=(const typename BM_Seq::BM_iterator& it)
			{
				if(this != &it){
					(*this)._field       = it._field;
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

			bool operator==(const BM_Seq::BM_iterator& it)
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
				typename Domain::Submatrix MC(M,0,i,rowd,1);
				typename Domain::Submatrix AC(A,0,i,rowd,1);
				Domain& MD = _seq.domain();
				MD.copy(MC,AC);
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
			Matrix  Algorithm3dot2(Matrix &D, std::vector<size_t> &d, size_t &mu, size_t &sigma, size_t &beta)
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
						size_t tempdeg = d[piv];
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
				matrix_type disc(field(),_row, _row+_col);
				//Create two iterators, one for seq, and one for gen
				typename BM_Seq::const_iterator cseqit;
				typename std::list<matrix_type>::iterator genit;
				//get a iterator to the seq element to be processed
				cseqit = _seqel;
				//Create a matrix domain for addition and multiplication

				Domain& MD = _seq->domain();
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
					_gen.push_back(matrix_type(field(),_col,_row+_col));
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
				matrix_type z1(field(),_col,_row+_col);
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
				matrix_type disc(field(),_row, _row+_col);
				//Create two iterators, one for seq, and one for gen
				typename BM_Seq::const_iterator cseqit;
				typename std::list<matrix_type>::iterator genit;
				//get an iterator to the seq element to be processed
				cseqit = _seqel;
				//Create a matrix domain for addition and multiplication
				Domain& MD = _seq.domain();
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
					_gen.push_back(matrix_type(field(),_col,_row+_col));
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
				matrix_type z1(field(),_col,_row+_col);
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
				std::vector<matrix_type> revgen(_mu+1, matrix_type(field(),_col,_col));
				for(size_t i = 0; i<_col; i++){
					typename std::list<matrix_type>::iterator genit = _gen.begin();
					for(int j = 0; j < _deg[i]+1; j++){
						ColumnCopy(revgen[_deg[i]-j], *genit,i);
						genit++;
					}
				}
				return revgen;
			}

			typename BM_Seq::size_type get_t()
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
			std::vector<size_t> get_deg()
			{
				std::vector<size_t> gendegree(&_deg[0], &_deg[_col-1]);
				return gendegree;
			}
		}; //End of BM_iterator

		//return an initialized BM_iterator
		typename BM_Seq::BM_iterator BM_begin()
		{
			return typename BM_Seq::BM_iterator(*this);
		}
		//return an initialized BM_iterator that points to one past the end of the sequence
		typename BM_Seq::BM_iterator BM_end()
		{
			return typename BM_Seq::BM_iterator(*this, _size);
		}
		/**/
	};//End of BM_Seq

    }; //end of class BlockCoppersmithDomain


                       template<class _Field, class _Sequence>
    std::vector<size_t>  BlockCoppersmithDomain<_Field,
_Sequence>::right_minpoly (std::vector<Coefficient> &P)
    {
	    //Get the row and column dimensions
	    const size_t r = _container->rowdim();
	    const size_t c = _container->coldim();

	    typename Sequence::const_iterator contiter(_container->begin());
	    //Create the BM_Seq, that will use the Coppersmith Block Berlekamp Massey Algorithm to compute the minimal generator.
	    BM_Seq seq(field(),r,c);

	    //Push the first projection onto the BM_Seq

	    seq.push_back(*contiter);

	    //Create the BM_Seq iterator whose incrementation performs a step of the generator
	    typename BM_Seq::BM_iterator bmit(seq.BM_begin());
	    bmit.setDelta(EARLY_TERM_THRESHOLD);
	    typename BM_Seq::BM_iterator::TerminationState check = bmit.state();
	    while(!check.IsGeneratorFound() ){
		    bmit++;
		    check = bmit.state();
		    if(check.IsSequenceExceeded()){
			    ++contiter;
			    seq.push_back(*contiter);
		    }
	    }
	    P = bmit.GetGenerator();
	    std::vector<size_t> deg(bmit.get_deg());
	    return deg;
    }

} // end of namespace LinBox

#endif // __LINBOX_coppersmith_block_domain_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

