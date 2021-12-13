/*
 * Copyright (C) 2020  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
#include "linbox/matrix/matrix-domain.h"

#include "linbox/field/hom.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include "givaro/modular.h"
#include <algorithm>

#ifdef TRACK_MEMORY_MATPOL
uint64_t max_memory=0, cur_memory=0;
#define ADD_MEM(x) {cur_memory+=x; max_memory=std::max(max_memory,cur_memory);}
#define DEL_MEM(x) {cur_memory-=x;}
#define STR_MEMINFO std::right<<"\033[31m [ MEM: cur="<<cur_memory/1000000.<<" Mo --- max="<<max_memory/1000000.<<" Mo \033[0m]"
#define PRINT_MEMINFO std::cerr<<"\033[31m[ MEM: cur="<<cur_memory/1000000.<<" Mo --- max="<<max_memory/1000000.<<" Mo ]\033[0m"<<std::endl;
#else
#define ADD_MEM(X) ;
#define DEL_MEM(X) ;     
#endif



#define COPY_BLOCKSIZE 32

namespace LinBox{

    
	enum PMType {polfirst, matfirst, matrowfirst};


    // matrix example
    /*  [ 1+2X+3X^2  4+5X+6X^2 ]
     *  [ 7+8X+9X^2  10+11X+12X^2 ]
     */
    
    /* polfirst storage : raw matrix of dimension (_row*_col) x _size
     *  [  1  2  3 ]
     *  [  4  5  6 ]  -> polynomial are stored contiguously (good for FFT)
     *  [  7  8  9 ]  -> storage is not BLAS compliant
     *  [ 10 11 12 ]
     *
     * matfirst storage : raw matrix of dimension _size x (_row*_col) (basically transposed of polfirst)
     *  [ 1 4 7 10 ]
     *  [ 2 5 8 11 ]  -> storage is BLAS compliant
     *  [ 3 6 9 12 ]  -> polynomial are not stored contiguously (gap of row*col, not good for FFT)
     *
     * matrowfirst storage: raw matrix of dimension _row x (_size*_col) 
     * [ 1  4 2  5 3  6 ]  -> storage is BLAS compliant (using stride)
     * [ 7 10 8 11 9 12 ]  -> polynomial are not stored contiguously (gap of row, might be good for FFT)
     *
     */
     
    // PG: the following does not make sense with characteristic 0 field !!! 
	template<typename Field> uint64_t element_storage(const Field& F)      {
        integer p;F.characteristic(p);
        if (p==0) throw LinboxError("PolynomialMatrix: use of element_storage with characterisrtic zero field is not allowed, aborting ...");
        return length(p);}
	template<> uint64_t element_storage(const Givaro::Modular<Givaro::Integer> &F) { integer p;F.characteristic(p); return length(p)+sizeof(Givaro::Integer);}
    
    
	// Generic handler class for Polynomial Matrix
	template<class Field, enum PMType>
	class PolynomialMatrix;
    
    template<class MatPoly>
	class SubPolynomialMatrix;
    
    
	// Class for Polynomial Matrix stored as a Matrix of Polynomials
	template<class _Field>
	class PolynomialMatrix<_Field, PMType::polfirst> {
    public:
        typedef _Field                         Field;
		typedef typename Field::Element      Element;
		typedef std::vector<Element,AlignedAllocator<Element, Alignment::DEFAULT>> VECT;
        typedef BlasVector<Field, VECT>                             Data;        

        typedef BlasMatrix<Field>                                   Matrix;
        typedef const BlasMatrix<Field>                        constMatrix;
        typedef BlasSubvector<Data>                             Polynomial;
        typedef const BlasSubvector<const Data>            constPolynomial;        

		typedef PolynomialMatrix<Field, PMType::polfirst>  Self_t;
        typedef SubPolynomialMatrix<Self_t>                  view;
        typedef SubPolynomialMatrix<const Self_t>      const_view;
        
		// construct a polynomial matrix in f[x]^(m x n) of degree (s-1)
		PolynomialMatrix(const Field& f, size_t r, size_t c, size_t s)
            : _rep(f,r*c*s), _row(r), _col(c), _size(s) {}


        /*******************************
         * Data access functionnalities*
         ******************************/
        
        // retrieve the polynomial at entry (i,j) in the matrix
        Polynomial       operator()(size_t i, size_t j)      {return operator()(i*_col+j);}
        constPolynomial  operator()(size_t i, size_t j)const {return operator()(i*_col+j);}

        // retrieve the polynomial at the position i in the storage of the matrix
        Polynomial       operator()(size_t i)      {return      Polynomial(_rep, i*_size,1,_size);}
        constPolynomial  operator()(size_t i)const {return constPolynomial(_rep, i*_size,1,_size);}

		// get write access to the the k-th coeff  of the ith matrix entry
        Element& ref(size_t i, size_t k)      {return _rep[i*_size+k];}
        // get read access to the the k-th coeff  of the ith matrix entry
        Element  get(size_t i, size_t k) const {return _rep[i*_size+k];}
        // get write access to the the k-th coeff  of the entry (i,j) in matrix 
        Element& ref(size_t i, size_t j, size_t k)      { return _rep[(i*_col+j)*_size+k];}
        // get read access to the the k-th coeff  of the entry (i,j) in matrix 
        Element  get(size_t i, size_t j, size_t k) const { return _rep[(i*_col+j)*_size+k];}

		Element*       getPointer()       {return _rep.getPointer();}
		const Element* getPointer() const {return _rep.getConstPointer();}

        view       at(size_t i, size_t j)      {return view(*this,i,j);}
        const_view at(size_t i, size_t j)const {return const_view(*this,i,j);}


        // initializee matrix entries at random
        void random(){
            _rep.random();
        }        
        template<class RandIter>
        void random(RandIter &I){
            _rep.random(I);
        }
        
        // copy the matrix of degree k into A
        template<typename DenseMatrix>
		DenseMatrix& getMatrix(DenseMatrix& A, size_t k) const {            
			auto it=A.Begin();
			for(size_t i=0;i<_row*_col;i++,it++)
				*it = get(i,k);
            return A;
		}

        // copy the matrix A into the matrix of degree k
        template<typename DenseMatrix>
		void setMatrix(const DenseMatrix& A, size_t k)  {
            //std::cout<<"SET MATRIX: "<<std::endl;
			auto it=A.Begin();
			for(size_t i=0;i<_row*_col;i++,it++){
                //std::cout<<*it<<" ";
                ref(i,k)=*it ;
            }
            //std::cout<<std::endl;
		}

        // retrieve the matrix of degree k in the polynomial matrix
		Matrix     operator[](size_t k) const {
			Matrix A(field(), _row, _col);
			getMatrix(A,k);
			return A;
		}
        
		// resize the polynomial length of the polynomial matrix
		void resize(size_t s){
			if (s==_size) return;
			if (s>_size){
				_rep.resize(s*_row*_col, field().zero);
				size_t k=s*_row*_col-1;
				for(size_t i=0;i<_row*_col;i++){
					size_t j=_size;
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
		}	
		
		
		// copy elt from M[beg..end], _size must be >= j-i
		void copy(const PolynomialMatrix<Field, PMType::polfirst>& M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....polfirst to polfirst.....same field"<<endl;
			for (size_t k=0;k<_row*_col;k++){
				size_t j=0;
				for (size_t i=beg;i<=end;i++){
					ref(k,start+j)=M.get(k,i);
					j++;
				}
			}
		}

        // copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Polynomial of Matrices
		void copy(const PolynomialMatrix<Field, PMType::matfirst>& M, size_t beg, size_t end, size_t start=0){
			//std::cout<<"copying.....matfirst to polfirst.....same field"<<std::endl;
			const size_t ls = COPY_BLOCKSIZE;
			for (size_t i = beg; i <= end; i+=ls)
				for (size_t j = 0; j < _col * _row; j+=ls)
					// Rk: the two loop must be interchanged in some cases
					for (size_t _i = i; _i < std::min(end+1, i + ls); _i++)
						for (size_t _j = j; _j < std::min(_col * _row, j + ls);++_j)
							ref(_j,start+_i-beg)= M.get(_j,_i);
		}

        // copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored with the hybrid format matrowfirst
		void copy(const PolynomialMatrix<Field, PMType::matrowfirst>& M, size_t beg, size_t end, size_t start=0){
			//std::cout<<"copying.....matrowfirst to polfirst.....same field"<<std::endl;
			const size_t ls = COPY_BLOCKSIZE;
            for(size_t k=0;k<_row;k++){
                for (size_t i = beg; i <= end; i+=ls)
                    for (size_t j = 0; j < _col; j+=ls)        
                        for (size_t _i = i; _i < std::min(end+1, i + ls); _i++)
                            for (size_t _j = j; _j < std::min(_col, j + ls);++_j)
                                ref(k,_j,start+_i-beg)= M.get(k,_j,_i);
            }
        }

        
		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Polynomial of Matrices with a different field
		template<typename OtherField>
		void copy(const PolynomialMatrix<OtherField, PMType::matfirst> & M, size_t beg, size_t end, size_t start=0){
			//std::cout<<"copying.....matfirst to polfirst.....other field"<<std::endl;
			const size_t ls = COPY_BLOCKSIZE;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for (size_t i = beg; i <= end; i+=ls)
				for (size_t j = 0; j < _col * _row; j+=ls)
					for (size_t _i = i; _i < std::min(end+1, i + ls); _i++) {
						for (size_t _j = j; _j < std::min(_col * _row, j + ls);++_j)
							hom.image(ref(_j,start+_i-beg),M.get(_j,_i) );
					}
		}



        template<typename OtherField>
		void copy(const PolynomialMatrix<OtherField, PMType::polfirst> & M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....polfirst to polfirst.....other field"<<endl;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for (size_t k=0;k<_row*_col;k++){
				size_t j=0;
				for (size_t i=beg;i<=end;i++,j++){
					hom.image(ref(k,start+j),M.get(k,i));
				}
			}
		}

        // copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored with the hybrid format matrowfirst
		template<typename OtherField>
		void copy(const PolynomialMatrix<OtherField, PMType::matrowfirst> & M, size_t beg, size_t end, size_t start=0){
			//std::cout<<"copying.....matrowfirst to polfirst.....same field"<<std::endl;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			const size_t ls = COPY_BLOCKSIZE;
            for(size_t k=0;k<_row;k++){
                for (size_t i = beg; i <= end; i+=ls)
                    for (size_t j = 0; j < _col; j+=ls)        
                        for (size_t _i = i; _i < std::min(end+1, i + ls); _i++)
                            for (size_t _j = j; _j < std::min(_col, j + ls);++_j)
                                hom.image(ref(k,_j,start+_i-beg), M.get(k,_j,_i));
            }
        }

	

		template<typename Mat>
		void copy(const Mat& M){
			linbox_check (M.size()==_size && M.rowdim()== _row && M.coldim()== _col);
			copy(M,0,M.size()-1);
		}


        template<class MatPoly>
        void copy(const SubPolynomialMatrix<MatPoly>& M, size_t beg, size_t end, size_t start=0){
            copy(*M._ptr, beg+M._shift, end+M._shift,start);
        }

        
		// rebind functor to change base field (e.g. apply modulo reduction)
		template<typename _Tp1>
		struct rebind {
			typedef PolynomialMatrix<Field, PMType::polfirst> Self_t;
			typedef PolynomialMatrix<_Tp1, PMType::polfirst> Other_t;

			void operator() (Other_t& Ap,
                             const Self_t&  A){
				Hom<Field, _Tp1> hom(A.field(), Ap.field()) ;
				for (size_t i = 0; i < A._row * A._col; i++)
					for (size_t j = 0; j < A._size; j++)
						hom.image (Ap.ref(i,j), A.get(i,j));
			}
		};


    
        size_t realmeminfo()const {			
			return _rep.size()*element_storage(field());
		}

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
            os<<"Matrix([";
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
                os << (i<_row-1?"],":"]" );//<< std::endl;
            }
			os<<"]);";
		
            return os;
        }

        size_t rowdim()  const {return _row;}
		size_t coldim()  const {return _col;}
		size_t degree()  const {return _size-1;}
		size_t size()    const {return _size;}
        size_t poly_stride() const {return 1;}
		const Field& field()  const {return _rep.field();}        


    private:
        Data     _rep;
        size_t   _row;
		size_t   _col;
		size_t  _size;
	};

    // Class for Polynomial Matrix stored as a Polynomial of Matrices
    template<class _Field>
	class PolynomialMatrix<_Field, PMType::matfirst> {
    public:
        typedef _Field                         Field;
		typedef typename Field::Element      Element;
		typedef std::vector<Element,AlignedAllocator<Element, Alignment::DEFAULT>> VECT;
        typedef BlasVector<Field, VECT>                             Data;        

        typedef BlasSubmatrix<BlasMatrix<Field>>                    Matrix;
        typedef BlasSubmatrix<const BlasMatrix<Field>>         constMatrix;
        typedef BlasSubvector<Data>                             Polynomial;
        typedef const BlasSubvector<const Data>            constPolynomial;        

		typedef PolynomialMatrix<_Field, PMType::matfirst>  Self_t;
        typedef SubPolynomialMatrix<Self_t>                   view;
        typedef SubPolynomialMatrix<const Self_t>       const_view;
        
		// construct a polynomial matrix in f[x]^(m x n) of degree (s-1)
        PolynomialMatrix(const Field& f, size_t r, size_t c, size_t s)
            : _rep(f,r*c*s), _row(r), _col(c), _size(s) {}


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

        // initializee matrix entries at random
        void random(){
            _rep.random();
        }        
        template<class RandIter>
        void random(RandIter &I){
            _rep.random(I);
        }

        // copy the matrix of degree k into A
        template<typename DenseMatrix>
		DenseMatrix& getMatrix(DenseMatrix& A, size_t k) const {
			return A= Matrix(field(), getPointer()+k*_row*_col, _row,_col,_col);            
        }

        // copy the matrix A into the matrix of degree
        template<typename DenseMatrix>
		void setMatrix(const DenseMatrix& A, size_t k)  {
            if (A.getPointer() != getPointer()+k*_row*_col){
                //std::cout<<"NEW MAtPol setMatrix: "<<A.getPointer()<<"<>"<<getPointer()+k*_row*_col<<std::endl;
                auto it=A.Begin();
                for(size_t i=0;i<_row*_col;i++,it++)
                    ref(i,k)=*it ;
            }
            else {
                //std::cout<<"NEW MAtPol setMatrix: does nothing\n";
            }
		} 


        Matrix       operator[](size_t k)          {return Matrix(field(), getPointer()+k*_row*_col, _row,_col,_col);}
		constMatrix  operator[](size_t k) const {return constMatrix(field(), getPointer()+k*_row*_col, _row,_col,_col);}

        view       at(size_t i, size_t j)       {return view(*this,i,j);}
        const_view at(size_t i, size_t j) const {return const_view(*this,i,j);}


		// resize the polynomial length of the polynomial matrix
		void resize(size_t s){
			_rep.resize(s*_row*_col,field().zero);
            _size=s;
		}
		
		// copy elt from M[beg..end], _size must be >= j-i
		void copy(const PolynomialMatrix<Field, PMType::matfirst>& M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....matfirst to matfirst.....same field"<<endl;
            FFLAS::fassign(field(), (end-beg+1)*_row*_col, M.getPointer()+beg*_row*_col,1, getPointer()+start,1);
		}

        // copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored with the hybrid format matrowfirst
		void copy(const PolynomialMatrix<Field, PMType::matrowfirst>& M, size_t beg, size_t end, size_t start=0){
			//std::cout<<"copying.....matrowfirst to polfirst.....same field"<<std::endl;
			const size_t ls = COPY_BLOCKSIZE;
            for (size_t i = beg; i <= end; i+=ls)
                for (size_t j = 0; j < _row; j+=ls)        
                    for (size_t _i = i; _i < std::min(end+1, i + ls); _i++)
                        for (size_t _j = j; _j < std::min(_row, j + ls);++_j)
                            for(size_t k=0;k<_col;k++)
                                ref(k,_j,start+_i-beg)= M.get(k,_j,_i);            
        }


        
		template<typename OtherField>
		void copy(const PolynomialMatrix<OtherField, PMType::matfirst> & M, size_t beg, size_t end, size_t start=0){
			//cout<<"copying.....matfirst to matfirst.....other field"<<endl;
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for (size_t i=beg;i<=end;i++)
				for (size_t j=0;j<_row*_col;j++)
					hom.image(ref(j,start+i-beg) , M.get(j,i));
		}

		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Matrix of Polynomials
		void copy(const PolynomialMatrix<Field, PMType::polfirst>& M, size_t beg, size_t end, size_t start=0){
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
		void copy(const PolynomialMatrix<OtherField, PMType::polfirst> & M, size_t beg, size_t end, size_t start=0){
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


        template<class MatPoly>
        void copy(const SubPolynomialMatrix<MatPoly>& M, size_t beg, size_t end, size_t start=0){
            copy(*M._ptr, beg+M._shift, end+M._shift,start);
        }

        
        
		// rebind functor to change base field (e.g. apply modulo reduction)
		template<typename _Tp1>
		struct rebind {
			typedef PolynomialMatrix<Field, PMType::matfirst>  Self_t;
			typedef PolynomialMatrix<_Tp1, PMType::matfirst>  Other_t;

			void operator() (PolynomialMatrix<_Tp1, PMType::matfirst>& Ap,
                             const PolynomialMatrix<Field, PMType::matfirst>&  A){
				Hom<Field, _Tp1> hom(A.field(), Ap.field()) ;
				for (size_t j = 0; j < A._size; j++)
					for (size_t i = 0; i < A._row * A._col; i++)
						hom.image (Ap.ref(i,j), A.get(i,j));
			}
		};

       
        size_t real_degree() const {
			MatrixDomain<Field> MD(field());
			size_t d= _size-1;
			while(d>0 && MD.isZero(operator[](d))) d--;
			return d;
		}

        size_t realmeminfo()const {			
			return _rep.size()*element_storage(field());
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
				wid=(int) ceil (log ((double) get(0,0,0)) / M_LN10);

			b= ((int) ceil (log ((double) (deg_max-deg_min+1)) / M_LN10));
			wid*=10*(b*(b-1)/2.);
			os<<"Matrix([";
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
                os << (i<_row-1?"],":"]" );//<< std::endl;
            }
			os<<"]);";
			return os;
        }        	

        size_t rowdim()  const {return _row;}
		size_t coldim()  const {return _col;}
		size_t degree()  const {return _size-1;}
		size_t size()    const {return _size;}
        size_t poly_stride() const {return _row*_col;}        
		const Field& field()  const {return _rep.field();}        

    private:
        Data     _rep;
        size_t   _row;
		size_t   _col;
		size_t  _size;
	};


    // Class for Polynomial Matrix stored with an hybrid variant that store each row as a Polynomial of row 
    template<class _Field>
	class PolynomialMatrix<_Field, PMType::matrowfirst> {
    public:
        typedef _Field                         Field;
		typedef typename Field::Element      Element;
		typedef std::vector<Element,AlignedAllocator<Element, Alignment::DEFAULT>> VECT;
        typedef BlasVector<Field, VECT>                               Data;        

        typedef BlasSubmatrix<BlasMatrix<Field>>                    Matrix;
        typedef BlasSubmatrix<const BlasMatrix<Field>>         constMatrix;
        typedef BlasSubvector<Data>                             Polynomial;
        typedef const BlasSubvector<const Data>            constPolynomial;        

		typedef PolynomialMatrix<_Field, PMType::matrowfirst>  Self_t;
        typedef SubPolynomialMatrix<Self_t>                      view;
        typedef SubPolynomialMatrix<const Self_t>          const_view;
        
		// construct a polynomial matrix in f[x]^(m x n) of degree (s-1)
        PolynomialMatrix(const Field& f, size_t r, size_t c, size_t s)
            : _rep(f,r*c*s), _row(r), _col(c), _size(s), _stride(c*s) {}


        /*******************************
         * Data access functionnalities*
         ******************************/
        
        // retrieve the polynomial at entry (i,j) in the matrix
		Polynomial       operator()(size_t i, size_t j)      {return operator()(i*_stride+j);}
		constPolynomial  operator()(size_t i, size_t j)const {return operator()(i*_stride+j);}

        // retrieve the polynomial at the position i in the storage of the matrix
		Polynomial       operator()(size_t i)      {return      Polynomial(_rep, i,_row,_size);}
		constPolynomial  operator()(size_t i)const {return constPolynomial(_rep, i,_row,_size);}
        
		// get write access to the the k-th coeff  of the ith matrix entry
		Element& ref(size_t i, size_t k)      {return ref(i/_col, i%_col, k);}
        // get read access to the the k-th coeff  of the ith matrix entry
		Element  get(size_t i, size_t k) const {return get(i/_col, i%_col, k);}
        // get write access to the the k-th coeff  of the entry (i,j) in matrix 
		Element& ref(size_t i, size_t j, size_t k)      { return _rep[i*_stride+j+k*_col];}
        // get read access to the the k-th coeff  of the entry (i,j) in matrix 
		Element  get(size_t i, size_t j, size_t k) const{ return _rep[i*_stride+j+k*_col];}

		Element*       getPointer()       {return _rep.getPointer();}
		const Element* getPointer() const {return _rep.getConstPointer();}
        
        // initializee matrix entries at random
        void random(){
            _rep.random();
        }        
        template<class RandIter>
        void random(RandIter &I){
            _rep.random(I);
        }
        
		// retrieve the matrix of degree k in the polynomial matrix
		//Matrix       operator[](size_t k)       {return      Matrix(field(), getPointer()+k*_row*_col, _row,_col,_col);}

        // copy the matrix of degree k into A
        template<typename DenseMatrix>
		DenseMatrix& getMatrix(DenseMatrix& A, size_t k) const {
			return A= Matrix(field(), getPointer()+k*_col, _row,_col,_stride);            
        }

        // copy the matrix A into the matrix of degree
        template<typename DenseMatrix>
		void setMatrix(const DenseMatrix& A, size_t k)  {
            if (A.getPointer() != getPointer()+k*_col){
                //std::cout<<"NEW MAtPol setMatrix: "<<A.getPointer()<<"<>"<<getPointer()+k*_row*_col<<std::endl;
                auto it=A.Begin();
                for(size_t i=0;i<_row*_col;i++,it++)
                    ref(i,k)=*it ;
            }
            else {
                //std::cout<<"NEW MAtPol setMatrix: does nothing\n";
            }
		} 

        //Matrix&&  operator[](size_t k)          {Matrix tmp(field(), getPointer()+k*_row*_col, _row,_col,_col); return std::move(tmp);}
        Matrix       operator[](size_t k)          {return Matrix(field(), getPointer()+k*_col, _row,_col,_stride);}
		constMatrix  operator[](size_t k) const {return constMatrix(field(), getPointer()+k*_col, _row,_col,_stride);}

        view       at(size_t i, size_t j)       {return view(*this,i,j);}
        const_view at(size_t i, size_t j) const {return const_view(*this,i,j);}




		// resize the polynomial length of the polynomial matrix
		void resize(size_t s){
			if (s==_size) return;
			if (s>_size ){
				_rep.resize(s*_row*_col, field().zero);

                if (s*_col > _stride) { // move data if needed
                    for(size_t i=0;i<_row;i++){
                        size_t k1=(i+1)*_stride-1;
                        size_t k2=(i+1)*s*_col-1;
                        for(size_t j=0;j<_size*_col;k1--,k2--)
                            _rep[k2]=_rep[k1];
                    }
                }
                // put zeros at the right place
                size_t diff=(s-_size)*_col;
                for(size_t i=0;i<_row;i++){
                    size_t k1=(i+1)*_stride;
                    FFLAS::fzero(field(), diff, getPointer()+k1, 1);
                }
                // correct the stride when necessary
                if (s*_col > _stride) _stride=s*_col;
			}
			// if (s< _size) do nothing since we have the right stride
			_size=s;
		}	

        /*
         * PG -> NEXT LINES NEED TO BE CHECKED AND TESTED
         *
         */

        // copy elt from M[beg..end], _size must be >= end_beg+1
		void copy(const PolynomialMatrix<Field, PMType::matfirst>& M, size_t beg, size_t end, size_t start=0){
            // 1st version, write access compliant
            for (size_t i=0;i<_row;i++)
                FFLAS::fassign(field(), (end-beg+1),_col, M.getPointer()+beg*_row*_col+i*_col, _row*_col, getPointer()+start+i*_stride,_col);
            // 2nd version, read access compliant
            // for (size_t i=beg;i<end;i++)
            //     FFLAS::fassign(field(), _row, _col, M.getPointer()+beg*_row*_col, _col, getPointer()+start+(beg-i),_stride);
            
		}        
		
		// copy elt from M[beg..end], _size must be >= end-beg+1
		void copy(const PolynomialMatrix<Field, PMType::matrowfirst>& M, size_t beg, size_t end, size_t start=0){
            FFLAS::fassign(field(), _row, (end-beg+1)*_col, M.getPointer()+beg*_col, M._stride, getPointer()+start,_stride);
		}
        
		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Matrix of Polynomials
		void copy(const PolynomialMatrix<Field, PMType::polfirst>& M, size_t beg, size_t end, size_t start=0){
			const size_t ls = COPY_BLOCKSIZE; // Loop unrooling to be more cache-friendly (cache block transposition)
            for(size_t k=0;k<_row;k++) {
                for(size_t i = beg; i < end+1; i+=ls)
                    for (size_t j = 0; j < _col; j+=ls)
                        for (size_t _i = i; _i < std::min(end+1, i + ls); _i++) {
                            for (size_t _j = j; _j < std::min(_col, j + ls); _j++)
                                ref(k,_j,start+_i-beg) = M.get(k,_j,_i);
                        }
            }
        }

		// copy elt from M[beg..end], _size must be >= end-beg+1
		// M is stored as a Matrix of Polynomials with a different field
        // PG -> could be specialized for each storage
		template<typename OtherField, PMType Storage>
		void copy(const PolynomialMatrix<OtherField, Storage> & M, size_t beg, size_t end, size_t start=0){
			Hom<OtherField,Field> hom(M.field(),field()) ;
			for (size_t i=beg;i<=end;i++)
				for (size_t j=0;j<_row*_col;j++)
					hom.image(ref(j,start+i-beg) , M.get(j,i));
		}



		template<typename Mat>
		void copy(const Mat& M){
			linbox_check (M.size()==_size && M.rowdim()== _row && M.coldim()== _col);
			copy(M,0,M.size()-1);
		}


        template<class MatPoly>
        void copy(const SubPolynomialMatrix<MatPoly>& M, size_t beg, size_t end, size_t start=0){
            copy(*M._ptr, beg+M._shift, end+M._shift,start);
        }

        
        
		// rebind functor to change base field (e.g. apply modulo reduction)
		template<typename _Tp1>
		struct rebind {
            typedef PolynomialMatrix<Field, PMType::matrowfirst>  Self_t;
            typedef PolynomialMatrix<_Tp1, PMType::matrowfirst>  Other_t;
            
            void operator() (PolynomialMatrix<_Tp1, PMType::matrowfirst>& Ap,
                             const PolynomialMatrix<Field, PMType::matrowfirst>&  A){
                Hom<Field, _Tp1> hom(A.field(), Ap.field()) ;
                for (size_t j = 0; j < A._size; j++)
                    for (size_t i = 0; i < A._row * A._col; i++)
                        hom.image (Ap.ref(i,j), A.get(i,j));
            }
        };


        size_t real_degree() const {
			MatrixDomain<Field> MD(field());
			size_t d= _size-1;
			while(d>0 && MD.isZero(operator[](d))) d--;
			return d;
		}

		size_t realmeminfo()const {			
			return _rep.size()*element_storage(field());
		}

        
		void dump(std::ostream& os) const {
			os<<_row<<" "<<_col<<" "<<_size<<std::endl;
			for(size_t k=0;k<_row;k++)
				for(size_t i=0;i<_size*_col;i++)		
					os<<getPointer()[k*_stride+i]<<" ";
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
				wid=(int) ceil (log ((double) get(0,0,0)) / M_LN10);

			b= ((int) ceil (log ((double) (deg_max-deg_min+1)) / M_LN10));
			wid*=10*(b*(b-1)/2.);
			os<<"Matrix([";
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
                os << (i<_row-1?"],":"]" );//<< std::endl;
            }
			os<<"]);";
			return os;
        }        	

        size_t rowdim()  const {return _row;}
		size_t coldim()  const {return _col;}
		size_t degree()  const {return _size-1;}
		size_t size()    const {return _size;}
        size_t poly_stride() const {return _row;}        
		const Field& field()  const {return _rep.field();}        

    private:
        Data     _rep;
        size_t   _row;
		size_t   _col;
		size_t  _size;
        size_t _stride;
	};


    
	template<typename _Field, LinBox::PMType T>
	std::ostream& operator<<(std::ostream& os, const PolynomialMatrix<_Field,T>& P) {
		return P.write(os);
	}
    


    template<typename MatPoly>
    struct SubPolyConst{
        typedef typename MatPoly::Matrix         Matrix;
        typedef typename MatPoly::constMatrix    constMatrix;
        typedef typename MatPoly::view             view;
        typedef typename MatPoly::const_view const_view;                
    };

    template<typename MatPoly>
    struct SubPolyConst <const MatPoly>{
        typedef typename MatPoly::constMatrix    constMatrix;
        typedef typename MatPoly::constMatrix         Matrix;
        typedef typename MatPoly::const_view        view;
        typedef typename MatPoly::const_view const_view;                
    };

    
    
	// Class to handle the view of a Polynomial Matrix according to some degree range
	template<class MatPoly>
	class SubPolynomialMatrix {
	public:
        typedef SubPolynomialMatrix<MatPoly>     Self_t;
        typedef typename MatPoly::Field           Field;
        typedef typename MatPoly::Element       Element;
        typedef typename MatPoly::Polynomial Polynomial;        
        typedef typename SubPolyConst<MatPoly>::Matrix            Matrix;
        typedef typename SubPolyConst<MatPoly>::constMatrix  constMatrix;
        typedef typename SubPolyConst<MatPoly>::view             view;
        typedef typename SubPolyConst<MatPoly>::const_view const_view;
        
		SubPolynomialMatrix() {}

		// constructor of a view between i and j from a plain Polynomial Matrix
		SubPolynomialMatrix(MatPoly& M, size_t i,size_t j)
			: _ptr(&M), _size(j-i+1), _shift(i)
		{
            linbox_check(i<M.size() && i<=j && j< M.size());
            //if (i>=M.size() || (i>j)) {_size=0;}
        }

		// constructor of a view between i and j from a Sub Polynomial Matrix
		SubPolynomialMatrix(Self_t & M, size_t i,size_t j)
			: _ptr(M._ptr), _size(j-i+1), _shift(i+M._shift)
		{
            linbox_check(i<M.size() && i<=j && j< M.size());
            //if (i>=M.size() || (i>j)) {_size=0;}
        }

        // copy the matrix of degree k into A
        template<typename DenseMatrix>
		DenseMatrix& getMatrix(DenseMatrix& A, size_t k) const {
            return _ptr->getMatrix(A,k+_shift);
        }

        // copy the matrix A into the matrix of degree 
        template<typename DenseMatrix>
		void setMatrix(const DenseMatrix& A, size_t k)  {
            return _ptr->setMatrix(A,k+_shift);
        }

        // only shrink down the size, increase is not possible
		void resize(size_t s){
            if (s<=_size) _size=s;
            else{
                std::cerr<<"LinBox -> SubPolynomialMatrix: resizing with larger size than available, exiting ..."<<std::endl;
                throw LinboxError("SubPolynomialMatrix: resizing with larger size than available, exiting ...");
            }
        }

        template <typename Mat>
        void copy(const Mat & M, size_t beg, size_t end){
            _ptr->copy(M,beg,end,_shift);
		}
        
		template<typename Mat>
        void copy(const Mat& M){
            copy(M,0,M.size()-1);
        }

        
		// retrieve the matrix of degree k in the polynomial matrix
        Matrix  operator[](size_t k) const {return _ptr->operator[](k+_shift);}

		// retrieve the polynomial at entry (i,j) in the matrix
        Polynomial    operator()(size_t i, size_t j){
            return Polynomial(_ptr->operator()(i,j),i, _ptr->poly_stride(), j-i+1);
        }

        Element get(size_t i, size_t k) const { return 	_ptr->get(i,k+_shift);}

        Element get(size_t i, size_t j, size_t k) const{ return get(i*coldim()+j,k);}

        size_t rowdim() const {return _ptr->rowdim();}
        size_t coldim() const {return _ptr->coldim();}
        size_t degree() const {return _size-1;}
        size_t size()   const {return _size;}
        const Field& field()  const {return _ptr->field();}


        const_view at(size_t i, size_t j) const {return const_view(*_ptr,i+_shift,j+_shift);}
        view       at(size_t i, size_t j)       {return       view(*_ptr,i+_shift,j+_shift);}
        
		std::ostream& write(std::ostream& os) const { if (_size!=0) return _ptr->write(os,_shift,_shift+_size-1); else return os;}

        

        // give friendship access to the main PolynomialMatrix classes
        template<class Field, enum PMType>
        friend class PolynomialMatrix;

	protected:
		MatPoly* _ptr;
		size_t  _size;
		size_t _shift;
	};

	template<typename MatPoly>      
	std::ostream& operator<<(std::ostream& os, const SubPolynomialMatrix<MatPoly>& P) {
		return P.write(os);
	}





} //end of namespace LinBox

#endif // __LINBOX_polynomial_matrix_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
