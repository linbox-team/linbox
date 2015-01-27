/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2013  Pascal Giorgi
 *                     Romain Lebreton
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *            Romain Lebreton <lebreton@lirmm.fr>
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
#ifndef __LINBOX_polynomial_matrix_H
#define __LINBOX_polynomial_matrix_H
#include <vector>
#include "linbox/vector/subvector.h"
#include <algorithm>
#include "linbox/matrix/dense-matrix.h"

#define COPY_BLOCKSIZE 32

namespace LinBox{

	enum PMType {polfirst, matfirst};
	enum PMStorage {plain, view, const_view};

	// Generic handler class for Polynomial Matrix
	template<size_t type, size_t storage, class Field>
	class PolynomialMatrix;

	// Class for Polynomial Matrix stored as a Matrix of Polynomials
	template<class _Field>
	class PolynomialMatrix<PMType::polfirst,PMStorage::plain,_Field> {
	public:
		typedef _Field Field;
		typedef typename Field::Element   Element;
		typedef BlasMatrix<Field>          Matrix;
		typedef typename std::vector<Element>::iterator  Iterator;
		typedef typename std::vector<Element>::const_iterator  ConstIterator;
		//typedef vector<Element>        Polynomial;
		typedef Subvector<Iterator,ConstIterator>   Polynomial;
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,_Field>  Self_t;
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,_Field> Other_t;

		PolynomialMatrix() {}

		// construct a polynomial matrix in f[x]^(m x n) of degree (s-1)
		PolynomialMatrix(const Field& f, size_t r, size_t c, size_t s) :
			_repview(r*c),_rep(r*c*s,f.zero), _row(r), _col(c), _size(s), _fld(&f) {
			for (size_t i=0;i<_row;i++)
				for (size_t j=0;j<_col;j++)
					_repview[i*_col+j]= Polynomial(_rep.begin()+(i*_col+j)*_size,_size);
		}

		// retrieve the matrix of degree k in the polynomial matrix
		Matrix     operator[](size_t k)const {
			Matrix A(field(), _row, _col);
			typename Matrix::Iterator it=A.Begin();
			for(size_t i=0;i<_row*_col;i++,it++)
				*it = get(i,k);
			return A;
		}
		// retrieve the polynomial at entry (i,j) in the matrix
		inline Polynomial&       operator()(size_t i, size_t j)      {return operator()(i*_col+j);}
		inline const Polynomial& operator()(size_t i, size_t j)const {return operator()(i*_col+j);}
		// retrieve the polynomial at the position i in the storage of the matrix
		inline Polynomial&       operator()(size_t i)      {return _repview[i];}
		inline const Polynomial& operator()(size_t i)const {return _repview[i];}

		// resize the polynomial length of the polynomial matrix
		void resize(size_t s){
			if (s>_size){
				_rep.resize(s*_row*_col);
				size_t k=s*_row*_col-1;
				for(size_t i=0;i<_row*_col;i++){
					size_t j=s;
					for(;j>=_size;j--,k--)
						_rep[k]=_fld->zero;
					for(;j>size_t(-1);j--,k--)
						_rep[k]=_rep[i*_size+j];
				}
			}
			else {
				size_t k=0;
				for(size_t i=0;i<_row*_col;i++)
					for (size_t j=0;j<s;j++,k++)
						_rep[k]=_rep[i*_size+j];
				_rep.resize(s*_row*_col);
			}
			_size=s;

			for (size_t i=0;i<_row;i++)
				for (size_t j=0;j<_col;j++)
					_repview[i*_col+j]= Polynomial(_rep.begin()+(i*_col+j)*_size,_size);
		}

		// copy elt from M[beg..end], _size must be >= j-i
		void copy(const Self_t& M, size_t beg, size_t end){
			//cout<<"copying.....polfirst to polfirst.....same field"<<endl;
			for (size_t k=0;k<_row*_col;k++){
				size_t j=0;
				for (size_t i=beg;i<=end;i++){
					ref(k,j)=M.get(k,i);
					j++;
				}
			}
		}
		template<typename OtherField>
		void copy(const PolynomialMatrix<PMType::polfirst,PMStorage::plain,OtherField> & M, size_t beg, size_t end){
			//cout<<"copying.....polfirst to polfirst.....other field"<<endl;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for (size_t k=0;k<_row*_col;k++){
				size_t j=0;
				for (size_t i=beg;i<=end;i++,j++){
					hom.image(ref(k,j),M.get(k,i));
				}
			}
		}

		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Polynomial of Matrices
		template <size_t storage>
		void copy(const PolynomialMatrix<PMType::matfirst,storage,Field>& M, size_t beg, size_t end){
			//cout<<"copying.....matfirst to polfirst.....same field"<<endl;
			const size_t ls = COPY_BLOCKSIZE;
			for (size_t i = beg; i <= end; i+=ls)
				for (size_t j = 0; j < _col * _row; j+=ls)
					// Rk: the two loop must be interchanged in some cases
					for (size_t _i = i; _i < std::min(end+1, i + ls); _i++)
						for (size_t _j = j; _j < std::min(_col * _row, j + ls);++_j)
							ref(_j,_i-beg)= M.get(_j,_i);
		}

		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Polynomial of Matrices with a different field
		template<size_t storage, typename OtherField>
		void copy(const PolynomialMatrix<PMType::matfirst,storage,OtherField> & M, size_t beg, size_t end){
			//cout<<"copying.....matfirst to polfirst.....other field"<<endl;
			const size_t ls = COPY_BLOCKSIZE;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for (size_t i = beg; i <= end; i+=ls)
				for (size_t j = 0; j < _col * _row; j+=ls)
					for (size_t _i = i; _i < std::min(end+1, i + ls); _i++) {
						for (size_t _j = j; _j < std::min(_col * _row, j + ls);++_j)
							hom.image(ref(_j,_i-beg),M.get(_j,_i) );
					}
		}

		template<typename Mat>
		void copy(const Mat& M){
			linbox_check (M.size()==_size && M.rowdim()== _row && M.coldim()== _col);
			copy(M,0,M.size()-1);
		}

		// rebind functor to change base field (e.g. apply modulo reduction)
		template<typename _Tp1>
		struct rebind {
			typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain, Field> Self_t;
			typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain, _Tp1> Other_t;

			void operator() (Other_t& Ap,
					 const Self_t&  A){
				Hom<Field, _Tp1> hom(A.field(), Ap.field()) ;
				for (size_t i = 0; i < A._row * A._col; i++)
					for (size_t j = 0; j < A._size; j++)
						hom.image (Ap.ref(i,j), A.get(i,j));
			}
		};

		// get access to the the k-th coeff  of the ith matrix entry
		inline Element& ref(size_t i, size_t k){
			return _rep[i*_size+k];
		}
		inline Element& ref(size_t i, size_t j, size_t k){
			return ref(i*_col+j,k);
		}
		inline Element get(size_t i, size_t k)const {
			return _rep[i*_size+k];
		}
		inline Element get(size_t i, size_t j, size_t k) const{
			return get(i*_col+j,k);
		}

		size_t rowdim() const {return _row;}
		size_t coldim() const {return _col;}
		size_t degree() const {return _size-1;}
		size_t size()   const {return _size;}
		const Field& field()  const {return *_fld;}

		std::ostream& write(std::ostream& os) const { return write(os,0,_size-1);}

                std::ostream& write(std::ostream& os, size_t deg_min, size_t deg_max) const {
                        integer c;
                        int wid,b;
                        field().cardinality (c);

                        if (c >0){
                                wid = (int) ceil (log ((double) c) / M_LN10);
                                b= ((int) ceil (log ((double) _size) / M_LN10));
                                wid*=10*(b*(b-1)/2.);
                        }
                        for (size_t i = 0; i< _row;++i) {
                                os << "  [ ";
                                for (size_t j = 0;j<_col;++j){
                                        os.width (wid);
                                        field().write(os,get(i+j*_col,deg_min));
					for (size_t k=deg_min+1;k<deg_max;++k){
                                                os<<"+";
                                                field().write(os,get(i+j*_col,k))<<"x^"<<k-deg_min;
                                        }
                                        if (deg_max-deg_min>0){
                                                os<<"+";
                                                field().write(os,get(i+j*_col,deg_max))<<"x^"<<deg_max-deg_min<<"  ";
                                        }
                                }
                                os << "]" << std::endl;
                        }
			return os;
                }

		Element* getWritePointer(){return &_rep[0];}
		const Element* getPointer() const {return &_rep[0];}

	private:
		std::vector<Polynomial> _repview;
		std::vector<Element>    _rep;
		size_t             _row;
		size_t             _col;
		size_t            _size;
		const Field*       _fld;
	};

	// Class for Polynomial Matrix stored as a Polynomial of Matrices
	template<class _Field>
	class PolynomialMatrix<PMType::matfirst,PMStorage::plain,_Field> {
	public:
		typedef _Field Field;
		typedef typename Field::Element   Element;
		typedef BlasMatrix<Field>         Matrix;
		typedef std::vector<Element>       Polynomial;
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,_Field>  Self_t;
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,_Field> Other_t;
		typedef Matrix value_type;

		PolynomialMatrix() {}

		PolynomialMatrix(const Field& f, size_t r, size_t c, size_t s) :
			 _rep(s,Matrix(f)), _row(r), _col(c), _size(s), _fld(&f) {
			for(size_t i=0;i<s;i++)
				_rep[i].init(f,r,c);
		}

		// retrieve the matrix of degree k in the polynomial matrix
		inline Matrix&       operator[](size_t k)      {return _rep[k];}
		inline const Matrix& operator[](size_t k)const {return _rep[k];}

		// retrieve the polynomial at entry (i,j) in the matrix
		inline Polynomial  operator()(size_t i, size_t j) const {
			Polynomial A(_size);
			for(size_t k=0;k<_size;k++)
				A[k]=_rep[k].refEntry(i,j);
			return A;
		}
		// retrieve the polynomial at psotion (i) in the storage of the matrix
		inline Polynomial  operator()(size_t i) const {
			Polynomial A(_size);
			for(size_t k=0;k<_size;k++)
				A[k]=*(_rep[k].Begin()+i);
			return A;
		}

		// resize the polynomial length of the polynomial matrix
		void resize(size_t s){
			_rep.resize(s,Matrix(*_fld));
			if (s>_size){
				for(size_t i=_size;i<s;i++)
					_rep[i].init(field(),_row,_col);
			}
			_size=s;
		}

		// copy elt from M[beg..end], _size must be >= j-i
		template <size_t storage>
		void copy(const PolynomialMatrix<PMType::matfirst,storage,Field>& M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....matfirst to matfirst.....same field"<<endl;
			for (size_t i=beg;i<=end;i++)
				_rep[start+i-beg]=M[i];

		}

		template<size_t storage, typename OtherField>
		void copy(const PolynomialMatrix<PMType::matfirst,storage,OtherField> & M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....matfirst to matfirst.....other field"<<endl;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for (size_t i=beg;i<=end;i++)
				for (size_t j=0;j<_row*_col;j++)
					hom.image(ref(j,start+i-beg) , M.get(j,i));
		}

		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Matrix of Polynomials
		void copy(const Other_t& M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....polfirst to matfirst.....same field"<<endl;
			const size_t ls = COPY_BLOCKSIZE; // Loop unrooling to be more cache-friendly (cache block transposition)
			for(size_t i = beg; i < end+1; i+=ls)
				for (size_t j = 0; j < _col * _row; j+=ls)
					for (size_t _i = i; _i < std::min(end+1, i + ls); _i++) {
						for (size_t _j = j; _j < std::min(_col * _row, j + ls); _j++)
							ref(_j,start+_i-beg) = M.get(_j,_i);
					}
		}

		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Matrix of Polynomials with a different field
		template<typename OtherField>
		void copy(const PolynomialMatrix<PMType::polfirst,PMStorage::plain,OtherField> & M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....polfirst to matfirst.....other field"<<endl;
			const size_t ls = COPY_BLOCKSIZE;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for(size_t i = beg; i < end+1; i+=ls)
				for (size_t j = 0; j < _col * _row; j+=ls)
					for (size_t _i = i; _i < std::min(end+1, i + ls); _i++) {
						for (size_t _j = j; _j < std::min(_col * _row, j + ls); _j++)
							hom.image(ref(_j,start+_i-beg) , M.get(_j,_i));
					}
		}

		template<typename Mat>
		void copy(const Mat& M){
			linbox_check (M.size()==_size && M.rowdim()== _row && M.coldim()== _col);
			copy(M,0,M.size()-1);
		}

		// rebind functor to change base field (e.g. apply modulo reduction)
		template<typename _Tp1>
		struct rebind {
			typedef PolynomialMatrix<PMType::matfirst, PMStorage::plain,Field>  Self_t;
			typedef PolynomialMatrix<PMType::matfirst, PMStorage::plain,_Tp1>  Other_t;

			template<size_t storage>
			void operator() (PolynomialMatrix<PMType::matfirst, storage,_Tp1>& Ap,
					 const PolynomialMatrix<PMType::matfirst, storage,Field>&  A){
				Hom<Field, _Tp1> hom(A.field(), Ap.field()) ;
				for (size_t j = 0; j < A._size; j++)
					for (size_t i = 0; i < A._row * A._col; i++)
						hom.image (Ap.ref(i,j), A.get(i,j));
			}
		};

		typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>           plain;
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::view,Field>             view;
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::const_view,Field> const_view;

		inline view       at(size_t i, size_t j)      {return view(*this,i,j);}
		inline const_view at(size_t i, size_t j)const {return const_view(*this,i,j);}


		// get access to the the k-th coeff  of the ith matrix entry
		inline Element& ref(size_t i, size_t k){
			return 	_rep[k].refRep()[i];
			//return *(_rep[k].Begin()+i);
		}
		inline Element& ref(size_t i, size_t j, size_t k){
			return ref(i*_col+j,k);
		}
		inline Element get(size_t i, size_t k)const {
			return *(_rep[k].Begin()+i);
		}
		inline Element get(size_t i, size_t j, size_t k) const{
			return get(i*_col+j,k);
		}

		inline size_t rowdim() const {return _row;}
		inline size_t coldim() const {return _col;}
		inline size_t degree() const {return _size-1;}
		inline size_t real_degree() const {
			MatrixDomain<Field> MD(field());
			size_t d= _size-1;
			while(d>0 && MD.isZero(_rep[d])) d--;
			return d;
		}
		inline size_t size()   const {return _size;}
		inline const Field& field()  const {return *_fld;}

		std::ostream& write(std::ostream& os) const { return write(os,0,_size-1);}

                std::ostream& write(std::ostream& os, size_t deg_min, size_t deg_max) const {
                        integer c;
                        int wid,b;
                        field().cardinality (c);

                        if (c >0)
                                wid = (int) ceil (log ((double) c) / M_LN10);
			else
				wid=(int) ceil (log ((double) _rep[0].getEntry(0,0)) / M_LN10);

			b= ((int) ceil (log ((double) _size) / M_LN10));
			wid*=10*(b*(b-1)/2.);

                        for (size_t i = 0; i< _row;++i) {
                                os << "  [ ";
                                for (size_t j = 0;j<_col;++j){
                                        os.width (wid);
                                        field().write(os,_rep[deg_min].getEntry(i,j));
                                        for (size_t k=deg_min+1;k<deg_max;++k){
                                                os<<"+";
                                                field().write(os,_rep[k].getEntry(i,j))<<"x^"<<k-deg_min;
                                        }
                                        if (deg_max-deg_min>0){
                                                os<<"+";
                                                field().write(os,_rep[deg_max].getEntry(i,j))<<"x^"<<deg_max-deg_min<<"  ";
                                        }
                                }
                                os << "]" << std::endl;
                        }
			return os;
                }

		// NEED FOR YUHASZ
		typedef typename std::vector<Matrix>::const_iterator const_iterator;
		const_iterator begin() const {return _rep.begin();}

	private:
		std::vector<Matrix>     _rep;
		size_t             _row;
		size_t             _col;
		size_t            _size;
		const Field*       _fld;
	};

	template<typename _Field, size_t T, size_t S>
	std::ostream& operator<<(std::ostream& os, const PolynomialMatrix<T,S,_Field>& P) {
		return P.write(os);
	}

	// Class to handle the view of a Polynomial Matrix according to some degree range
	// the view is read/write
	template<class _Field>
	class PolynomialMatrix<PMType::matfirst,PMStorage::view,_Field> {
	public:
		typedef _Field Field;
		typedef typename Field::Element   Element;
		typedef BlasMatrix<Field>          Matrix;
		typedef std::vector<Element>        Polynomial;

		PolynomialMatrix() {}

		// constructor of a view between i and j from a plain Polynomial Matrix
		PolynomialMatrix(PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>& M, size_t i,size_t j)
			: _ptr(&M), _size(j-i+1), _shift(i)
		{//cout<<i<<"-"<<j<<" -> "<<M.size()<<endl;
			linbox_check(i<M.size() && i<=j && j< M.size());}

		// constructor of a view between i and j from a view Polynomial Matrix
		PolynomialMatrix(PolynomialMatrix<PMType::matfirst,PMStorage::view,Field>& M, size_t i,size_t j)
			: _ptr(M._ptr), _size(j-i+1), _shift(i+M._shift)
		{//cout<<i<<"-"<<j<<" -> "<<M.size()<<endl;
			linbox_check(i<M.size() && i<=j && j< M.size());}


		// retrieve the matrix of degree k in the polynomial matrix
		inline Matrix&       operator[](size_t k)      {return (*_ptr)[k+_shift];}
		inline const Matrix& operator[](size_t k)const {return (*_ptr)[k+_shift];}

		// retrieve the polynomial at entry (i,j) in the matrix
		Polynomial    operator()(size_t i, size_t j){
			Polynomial A(_size, Matrix(_ptr->field()));
			for(size_t k=0;k<_size;k++)
				A[k]=(*_ptr)[k+_shift].refEntry(i,j);
			return A;
		}

		// get access to the the k-th coeff  of the ith matrix entry
		inline Element& ref(size_t i, size_t k){
			return 	_ptr->ref(i,k+_shift);
		}
		inline Element& ref(size_t i, size_t j, size_t k){
			return ref(i*coldim()+j,k);
		}
		inline Element get(size_t i, size_t k)const {
			return 	_ptr->get(i,k+_shift);
		}
		inline Element get(size_t i, size_t j, size_t k) const{
			return get(i*coldim()+j,k);
		}

		inline size_t rowdim() const {return _ptr->rowdim();}
		inline size_t coldim() const {return _ptr->coldim();}
		inline size_t degree() const {return _size-1;}
		inline size_t size()   const {return _size;}
		inline const Field& field()  const {return _ptr->field();}

		typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>           plain;
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::view,Field>             view;
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::const_view,Field> const_view;

		template<typename Mat>
		void copy(const Mat& M){
			this->copy(M,0,M.size()-1);
		}

		template<typename Mat>
		void copy(const Mat& M, size_t beg, size_t end){_ptr->copy(M,beg,end,_shift);}


		inline view       at(size_t i, size_t j)      {return view(*this,i,j);}
		inline const_view at(size_t i, size_t j)const {return const_view(*this,i,j);}

		std::ostream& write(std::ostream& os) const { return _ptr->write(os,_shift,_shift+_size-1);}

		friend class  PolynomialMatrix<PMType::matfirst,PMStorage::const_view,Field>;

	private:
		PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>* _ptr;
		size_t  _size;
		size_t _shift;
	};

	// Class to handle the view of a Polynomial Matrix according to some degree range
	// the view is read only
	template<class _Field>
	class PolynomialMatrix<PMType::matfirst,PMStorage::const_view,_Field> {
	public:
		typedef _Field Field;
		typedef typename Field::Element   Element;
		typedef BlasMatrix<Field>          Matrix;
		typedef std::vector<Element>        Polynomial;

		PolynomialMatrix() {}

		// constructor of a view between i and j from a plain Polynomial Matrix
		PolynomialMatrix(const PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>& M, size_t i,size_t j)
			: _ptr(&M), _size(j-i+1), _shift(i)
		{linbox_check(i<M.size() && i<=j && j< M.size());}

		// constructor of a view between i and j from a view Polynomial Matrix
		PolynomialMatrix(const PolynomialMatrix<PMType::matfirst,PMStorage::view,Field>& M, size_t i,size_t j)
			: _ptr(M._ptr), _size(j-i+1), _shift(i+M._shift)
		{linbox_check(i<M.size() && i<=j && j< M.size());}


		// constructor of a view between i and j from a view Polynomial Matrix
		PolynomialMatrix(const PolynomialMatrix<PMType::matfirst,PMStorage::const_view,Field>& M, size_t i,size_t j)
			: _ptr(M._ptr), _size(j-i+1), _shift(i+M._shift)
		{//cout<<"("<<i<<","<<j<<") -> "<<M.size()<<endl;
			linbox_check(i<M.size()&& i<=j && j<= M.size());}


		// retrieve the matrix of degree k in the polynomial matrix
		inline const Matrix& operator[](size_t k)const {return (*_ptr)[k+_shift];}

		// retrieve the polynomial at entry (i,j) in the matrix
		inline Polynomial    operator()(size_t i, size_t j){
			Polynomial A(_size, Matrix(_ptr->field()));
			for(size_t k=0;k<_size;k++)
				A[k]=(*_ptr)[k+_shift].refEntry(i,j);
			return A;
		}

		inline Element get(size_t i, size_t k) const {
			return 	_ptr->get(i,k+_shift);
		}
		inline Element get(size_t i, size_t j, size_t k) const{
			return get(i*coldim()+j,k);
		}


		inline size_t rowdim() const {return _ptr->rowdim();}
		inline size_t coldim() const {return _ptr->coldim();}
		inline size_t degree() const {return _size-1;}
		inline size_t size()   const {return _size;}
		inline const Field& field()  const {return _ptr->field();}

		typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>           plain;
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::const_view,Field> const_view;

		inline const_view at(size_t i, size_t j)const {return const_view(*this,i,j);}
		std::ostream& write(std::ostream& os) const { return _ptr->write(os,_shift,_shift+_size-1);}

	private:
		const PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>* _ptr;
		size_t  _size;
		size_t _shift;
	};




} //end of namespace LinBox

#endif // __LINBOX_polynomial_matrix_H
