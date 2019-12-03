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
#include "linbox/vector/vector.h"
#include "linbox/matrix/dense-matrix.h"

#include "linbox/field/hom.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include "givaro/modular.h"
#include <algorithm>




#define COPY_BLOCKSIZE 32

namespace LinBox{

	enum PMType {polfirst, matfirst};
	//enum PMStorage {plain, view, const_view};

	// Generic handler class for Polynomial Matrix
	template<class Field, enum PMType>
	class PolynomialMatrix;
    
    template<class MatPoly>
	class SubPolynomialMatrix;
    
	template<class _Field>
	class PolynomialMatrixBase<_Field> {
    public:
        typedef _Field                         Field;
		typedef typename Field::Element      Element;
		typedef std::vector<Element,AlignedAllocator<Element, Alignment::DEFAULT>> VECT;
        typedef BlasVector<Field>                    Data;        
        typedef BlasMatrix<Field>                           Matrix;
        typedef BlasSubmatrix<Matrix>                    SubMatrix;
        typedef BlasSubmatrix<const Matrix>         constSubMatrix;
        typedef BlasVector<Field>                       Polynomial;
        typedef BlasSubvector<Data>                  SubPolynomial;
        typedef const BlasSubvector<const Data> constSubPolynomial;        

		// construct a polynomial matrix in f[x]^(m x n) of degree (s-1)
		PolynomialMatrixBase(const Field& f, size_t r, size_t c, size_t s)
            : _rep(f,r*c*s), _row(r), _col(c), _size(s) {}



        size_t rowdim()  const {return _row;}
		size_t coldim()  const {return _col;}
		size_t degree()  const {return _size-1;}
		size_t size()    const {return _size;}
		size_t storage() const {return _store;}
        
		const Field& field()  const {return _rep.field();}        


    protected:
        Data     _rep;
        size_t   _row;
		size_t   _col;
		size_t  _size;
        size_t _store;
    };
    
	// Class for Polynomial Matrix stored as a Matrix of Polynomials
	template<class _Field>
	class PolynomialMatrix<_Field, PMType::polfirst>
        : public PolynomialMatrixBase<Field> {
    public:
		typedef PolynomialMatrix<_Field, PMType::polfirs>  Self_t;
        using typename PolynomialMatrixBase<Field>::Field;
        using typename PolynomialMatrixBase<Field>::Element;
        typedef typename PolynomialMatrixBase<Field>::Matrix                      Matrix;
        typedef typename PolynomialMatrixBase<Field>::SubPolynomial           Polynomial;
        typedef typename PolynomialMatrixBase<Field>::constSubPolynomial constPolynomial;
        typedef SubPolynomialMatrix<Selft>                                          view;
        typedef SubPolynomialMatrix<const Selft>                              const_view;
        
		// construct a polynomial matrix in f[x]^(m x n) of degree (s-1)
		//PolynomialMatrix(const Field& f, size_t r, size_t c, size_t s, size_t stor=0);
        using typename PolynomialMatrixBase<_Field>::PolynomialMatrixBase<_Field>;

        /*******************************
         * Data access functionnalities*
         ******************************/
        
        // retrieve the polynomial at entry (i,j) in the matrix
        Polynomial&       operator()(size_t i, size_t j)      {return operator()(i*_col+j);}
        constPolynomial&  operator()(size_t i, size_t j)const {return operator()(i*_col+j);}

        // retrieve the polynomial at the position i in the storage of the matrix
        Polynomial&       operator()(size_t i)      {return      Polynomial(_rep, i*_store,1,_size);}
        constPolynomial&  operator()(size_t i)const {return constPolynomial(_rep, i*_store,1,_size);}

		// get write access to the the k-th coeff  of the ith matrix entry
        Element& ref(size_t i, size_t k)      {return _rep[i*_store+k];}
        // get read access to the the k-th coeff  of the ith matrix entry
        Element  get(size_t i, size_t k) const {return _rep[i*_store+k];}
        // get write access to the the k-th coeff  of the entry (i,j) in matrix 
        Element& ref(size_t i, size_t j, size_t k)      { return _rep[(i*_col+j)*_store+k];}
        // get read access to the the k-th coeff  of the entry (i,j) in matrix 
        Element  get(size_t i, size_t j, size_t k) const { return _rep[(i*_col+j)*_store+k];}

		Element*       getPointer()       {return _rep.getPointer();}
		const Element* getPointer() const {return _rep.getConstPointer();}

        view       at(size_t i, size_t j)      {return view(*this,i,j);}
        const_view at(size_t i, size_t j)const {return const_view(*this,i,j);}

        // set the matrix A to the matrix of degree k in the polynomial matrix
		void setMatrix(Matrix& A, size_t k) const {
			typename Matrix::Iterator it=A.Begin();
			for(size_t i=0;i<_row*_col;i++,it++)
				*it = get(i,k);			
		}
		// retrieve the matrix of degree k in the polynomial matrix
		Matrix     operator[](size_t k) const {
			Matrix A(field(), _row, _col);
			setMatrix(A,k);
			return A;
		}
        
		// resize the polynomial length of the polynomial matrix
		void resize(size_t s){
			if (s==_store) return;
			if (s>_store){
				_rep.resize(s*_row*_col, field().zero);
				size_t k=s*_row*_col-1;
				for(size_t i=0;i<_row*_col;i++){
					size_t j=_size;
					for(;j>size_t(-1);j--,k--)
						_rep[k]=_rep[i*_store+j];
				}
			}
			else {
				size_t k=0;
				for(size_t i=0;i<_row*_col;i++)
					for (size_t j=0;j<s;j++,k++)
						_rep[k]=_rep[i*_store+j];
				_rep.resize(s*_row*_col);
			}
			_store=s;
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
		void copy(const PolynomialMatrix<PMType::polfirst,OtherField> & M, size_t beg, size_t end){
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
			//std::cout<<"copying.....matfirst to polfirst.....same field"<<std::endl;
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
			//std::cout<<"copying.....matfirst to polfirst.....other field"<<std::endl;
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
			typedef PolynomialMatrix<PMType::polfirst, Field> Self_t;
			typedef PolynomialMatrix<PMType::polfirst, _Tp1> Other_t;

			void operator() (Other_t& Ap,
                             const Self_t&  A){
				Hom<Field, _Tp1> hom(A.field(), Ap.field()) ;
				for (size_t i = 0; i < A._row * A._col; i++)
					for (size_t j = 0; j < A._size; j++)
						hom.image (Ap.ref(i,j), A.get(i,j));
			}
		};


		void dump(std::ostream& os) const {
			os<<_row<<" "<<_col<<" "<<_size<<std::endl;
			for(size_t i=0;i<_row*_col;i++)
				for(size_t k=0;k<_size;k++)
					os<<get(i,k)<<" ";

			os<<std::endl;
		}
		
		std::ostream& write(std::ostream& os) const { return write(os,0,_size-1);}

        std::ostream& write(std::ostream& os, size_t deg_min, size_t deg_max) const {
            integer c;
            int wid=-1,b;
            field().cardinality (c);

            if (c >0){
                wid = (int) ceil (log ((double) c) / M_LN10);
                b= ((int) ceil (log ((double) _size) / M_LN10));
                wid*=10*(b*(b-1)/2.);
            }
            std::cout<<"Matrix([";
            for (size_t i = 0; i< _row;++i) {
                os << "  [ ";
                for (size_t j = 0;j<_col;++j){
                    os.width (wid);
                    field().write(os,get(i,j,deg_min));
                    for (size_t k=deg_min+1;k<deg_max;++k){
                        os<<"+";
                        field().write(os,get(i,j,k))<<"*x^"<<k-deg_min;
                    }
                    if (deg_max-deg_min>0){
                        os<<"+";
                        field().write(os,get(i,j,deg_max))<<"*x^"<<deg_max-deg_min;
                    }
                    os<<(j<_col-1?",":"" );
                }
                os << (i<_row-1?"],":"]" )<< std::endl;
            }
			std::cout<<"]);";
		
            return os;
        }
	};

	// Class for Polynomial Matrix stored as a Polynomial of Matrices
	template<class _Field>
	class PolynomialMatrix<PMType::matfirst,_Field>
        : public PolynomialMatrixBase<Field> {
    public:
		typedef PolynomialMatrix<_Field, PMType::polfirst>  Self_t;
        using typename PolynomialMatrixBase<Field>::Field;
        using typename PolynomialMatrixBase<Field>::Element;

        typedef typename PolynomialMatrixBase<Field>::subMatrix                   Matrix;
        typedef typename PolynomialMatrixBase<Field>::constSubMatrix         constMatrix;
        typedef typename PolynomialMatrixBase<Field>::SubPolynomial           Polynomial;
        typedef typename PolynomialMatrixBase<Field>::constSubPolynomial constPolynomial; 
        typedef SubPolynomialMatrix<Selft>                                          view;
        typedef SubPolynomialMatrix<const Selft>                              const_view;
        
		// construct a polynomial matrix in f[x]^(m x n) of degree (s-1)
		//PolynomialMatrix(const Field& f, size_t r, size_t c, size_t s, size_t stor=0);
        using typename PolynomialMatrixBase<_Field>::PolynomialMatrixBase<_Field>;
       

        /*******************************
         * Data access functionnalities*
         ******************************/
        
        // retrieve the polynomial at entry (i,j) in the matrix
		Polynomial       operator()(size_t i, size_t j)      {return operator()(i*_col+j);}
		constPolynomial  operator()(size_t i, size_t j)const {return operator()(i*_col+j);}

        // retrieve the polynomial at the position i in the storage of the matrix
		Polynomial       operator()(size_t i)      {return      Polynomial(_rep, i,_row*_col,_size);}
		constPolynomial  operator()(size_t i)const {return constPolynomial(_rep, i,_row*_col,_size);}

		// get write access to the the k-th coeff  of the ith matrix entry
		Element& ref(size_t i, size_t k)      {return _rep[i+k*_row*_col];}
        // get read access to the the k-th coeff  of the ith matrix entry
		Element  get(size_t i, size_t k) const {return _rep[i+k*_row*_col];}
        // get write access to the the k-th coeff  of the entry (i,j) in matrix 
		Element& ref(size_t i, size_t j, size_t k)      { return _rep[i*_col+j+k*_row*_col];}
        // get read access to the the k-th coeff  of the entry (i,j) in matrix 
		Element  get(size_t i, size_t j, size_t k) const { return _rep[i*_col+j+k*_row*_col];}

		Element*       getPointer()       {return _rep.getPointer();}
		const Element* getPointer() const {return _rep.getConstPointer();}

		
		// retrieve the matrix of degree k in the polynomial matrix
		Matrix       operator[](size_t k)       {return      Matrix(field(), getPointer()+k*_row*_col, _row,_col,_col);}
		constMatrix  operator[](size_t k) const {return constMatrix(field(), getPointer()+k*_row*_col, _row,_col,_col);}

        view       at(size_t i, size_t j)       {return view(*this,i,j);}
        const_view at(size_t i, size_t j) const {return const_view(*this,i,j);}


		// resize the polynomial length of the polynomial matrix
		void resize(size_t s){
			_rep.resize(s,field().zero);
		}
		
		// copy elt from M[beg..end], _size must be >= j-i
		template <size_t storage>
		void copy(const PolynomialMatrix<PMType::matfirst,storage,Field>& M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....matfirst to matfirst.....same field"<<endl;
            FFLAS::fassign(field(), (end-beg+1)*_row*_col, M.getPointer()+beg*_row*_col,1, getPointer()+start,1);
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
		void copy(const PolynomialMatrix<PMType::polfirst,OtherField> & M, size_t beg, size_t end, size_t start=0){
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
			typedef PolynomialMatrix<PMType::matfirst,Field>  Self_t;
			typedef PolynomialMatrix<PMType::matfirst,_Tp1>  Other_t;

			template<size_t storage>
			void operator() (PolynomialMatrix<PMType::matfirst, storage,_Tp1>& Ap,
                             const PolynomialMatrix<PMType::matfirst, storage,Field>&  A){
				Hom<Field, _Tp1> hom(A.field(), Ap.field()) ;
				for (size_t j = 0; j < A._size; j++)
					for (size_t i = 0; i < A._row * A._col; i++)
						hom.image (Ap.ref(i,j), A.get(i,j));
			}
		};

       
        size_t real_degree() const {
			MatrixDomain<Field> MD(field());
			size_t d= _size-1;
			while(d>0 && MD.isZero(operator[d])) d--;
			return d;
		}

        
		void dump(std::ostream& os) const {
			os<<_row<<" "<<_col<<" "<<_size<<std::endl;
			for(size_t k=0;k<_size;k++)
				for(size_t i=0;i<_row*_col;i++)		
					os<<get(i,k)<<" ";
			os<<std::endl;
		}

		
		std::ostream& write(std::ostream& os) const { return write(os,0,real_degree());}

        std::ostream& write(std::ostream& os, size_t deg_min, size_t deg_max) const {
            integer c;
            int wid,b;
            field().cardinality (c);

            if (c >0)
                wid = (int) ceil (log ((double) c) / M_LN10);
			else
				wid=(int) ceil (log ((double) _rep[0].getEntry(0,0)) / M_LN10);

			b= ((int) ceil (log ((double) (deg_max-deg_min+1)) / M_LN10));
			wid*=10*(b*(b-1)/2.);
			std::cout<<"Matrix([";
			for (size_t i = 0; i< _row;++i) {
                os << "  [ ";
                for (size_t j = 0;j<_col;++j){
                    os.width (wid);
                    //field().write(os,_rep[deg_min].getEntry(i,j));
                    field().write(os,get(i,j,deg_min));
                    for (size_t k=deg_min+1;k<deg_max;++k){
                        if (!field().isZero(get(i,j,k))){
                            os<<"+";
                            field().write(os,get(i,j,k))<<"*x^"<<k-deg_min;
                        }
                    }
                    if (deg_max-deg_min>0){
                        if (!field().isZero(get(i,j,deg_max))){
                            os<<"+";
                            field().write(os,get(i,j,deg_max))<<"*x^"<<deg_max-deg_min;
                        }
                    }
					os<<(j<_col-1?",":"" );
                }
                os << (i<_row-1?"],":"]" )<< std::endl;
            }
			std::cout<<"]);";
			return os;
        }
        
		
		// NEED FOR YUHASZ
		typedef typename std::vector<Matrix>::const_iterator const_iterator;
		const_iterator begin() const {return _rep.begin();}
        
	};

	template<typename _Field, size_t T, size_t S>
	std::ostream& operator<<(std::ostream& os, const PolynomialMatrix<T,S,_Field>& P) {
		return P.write(os);
	}
    

	// Class to handle the view of a Polynomial Matrix according to some degree range
	template<class MatPoly>
	class SubPolynomialMatrix {
	public:
        typedef SubPolynomialMatrix<MatPoly> Self_t;
        typedef typename MaPoly::Matrix Matrix;
        typedef typename MaPoly::Polynomial Polynomial;
        typedef typename MatPoly::view view;
        typedef typename MatPoly::const_view const_view;
        
		SubPolynomialMatrix() {}

		// constructor of a view between i and j from a plain Polynomial Matrix
		SubPolynomialMatrix(MatPoly& M, size_t i,size_t j)
			: _ptr(&M), _size(j-i+1), _shift(i)
		{linbox_check(i<M.size() && i<=j && j< M.size());}

		// constructor of a view between i and j from a Sub Polynomial Matrix
		SubPolynomialMatrix(Self_t & M, size_t i,size_t j)
			: _ptr(M._ptr), _size(j-i+1), _shift(i+M._shift)
		{linbox_check(i<M.size() && i<=j && j< M.size());}

		// retrieve the matrix of degree k in the polynomial matrix
        Matrix  operator[](size_t k)const {return (*_ptr)[k+_shift];}

		// retrieve the polynomial at entry (i,j) in the matrix
        Polynomial    operator()(size_t i, size_t j){
            return Polynomial(??? );
        }
        
        //     return Polynomial(field()
		// 	Polynomial A(_size, Matrix(_ptr->field()));
		// 	for(size_t k=0;k<_size;k++)
		// 		A[k]=(*_ptr)[k+_shift].refEntry(i,j);
		// 	return A;
		// }

        Element get(size_t i, size_t k) const { return 	_ptr->get(i,k+_shift);}

        Element get(size_t i, size_t j, size_t k) const{ return get(i*coldim()+j,k);}

        size_t rowdim() const {return _ptr->rowdim();}
        size_t coldim() const {return _ptr->coldim();}
        size_t degree() const {return _size-1;}
        size_t size()   const {return _size;}
        const Field& field()  const {return _ptr->field();}


        const_view at(size_t i, size_t j) const {return const_view(*_ptr,i+_shift,j+_shift);}
        view       at(size_t i, size_t j)       {return       view(*_ptr,i+_shift,j+_shift);}
        
		std::ostream& write(std::ostream& os) const { return _ptr->write(os,_shift,_shift+_size-1);}

	private:
		MatPoly* _ptr;
		size_t  _size;
		size_t _shift;
	};





} //end of namespace LinBox

#endif // __LINBOX_polynomial_matrix_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
