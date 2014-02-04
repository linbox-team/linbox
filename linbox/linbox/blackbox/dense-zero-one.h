/* linbox/blackbox/dense-zero-one.h
 * Copyright (c) LinBox
 * written by Nick Messina <nmessina@cis.udel.edu>
 *
 * -----------------------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_dense_zero_one_H
#define __LINBOX_dense_zero_one_H



#include <iostream>
//#include "linbox/field/hom.h"
#include "linbox/linbox-config.h"
#include <vector>
#include "linbox/blackbox/blackbox-interface.h"


//const size_t SIZE = sizeof(LinBox::uint64_t) * 8;
#define SIZE  (sizeof(uint64_t) * 8)


namespace LinBox
{

	template <class _MatrixDomain>
	class DenseZeroOne : public  BlackboxInterface {

	public:

		// The following two typedefs are not useful in LinBox yet, as matrix domains aren't yet defined
		typedef _MatrixDomain MatrixDomain;
		typedef typename MatrixDomain::Scalar Scalar;
		typedef typename MatrixDomain::Block Block;
		typedef MatrixDomain Field;

		typedef uint64_t PackedUnit;
		typedef DenseZeroOne<MatrixDomain> Self_t;

		/* Notes:
		 * _rows gives the number of rows of PackedUnits
		 * _cols gives the number of bits in each row */

		DenseZeroOne(MatrixDomain MD, size_t m, size_t n) :
			_MD(MD), _rows(m), _cols(n)
		{
			_valPerWord = SIZE;
			_mask = 1;
			_rep.resize( _rows * ( ((_cols-1)/_valPerWord) + 1));
		}

		DenseZeroOne(const Self_t & other) :
			_MD(other._MD), _rows(other._rows), _cols(other._cols), _valPerWord(other._valPerWord), _rep(other._rep), _mask(other._mask)
		{}

		/* Constructor for making the blackbox out of a vector of PackedUnits.
		 * Note that if the vector given is too short for given size(m by n),
		 * the blackbox is enlarged and padded by 0's added after the given vector.
		 * If the vector given is too large, the blackbox truncates this vector
		 * to the appropriate number of PackedUnits.  */
		template <class Vector>
		DenseZeroOne(MatrixDomain MD, size_t m, size_t n, Vector & v) :
			_MD(MD), _rows(m), _cols(n)
		{
			_valPerWord = SIZE;
			_mask = 1;
			_rep = v;
			_rep.resize( _rows * ( ((_cols-1)/_valPerWord) + 1));
		}


		Scalar & getEntry(Scalar & x, size_t i, size_t j) const
		{
			//Get all index values needed.
			size_t unitNum = i*(_rep.size()/_rows) + (j/_valPerWord);  //Starts from 0
			size_t bitPlace = (j%_valPerWord) + 1;                    //Starts from 1

			_MD.init(x, (_rep[unitNum] >> (_valPerWord - bitPlace) ) & _mask );
			return x;
		}

		void setEntry(size_t i, size_t j, const Scalar & x)
		{
			//Get all index values needed.
			size_t unitNum = i*(_rep.size()/_rows) + (j/_valPerWord);  //Starts from 0
			size_t bitPlace = (j%_valPerWord) + 1;                    //Starts from 1

			PackedUnit y = 0;
			if (_MD.isOne(x))
				y = 1;

			_rep[unitNum] &= ~(_mask << (_valPerWord - bitPlace));
			_rep[unitNum] |= (y << (_valPerWord - bitPlace));
		}

		//Apply Variants

		// matrix-vector product
		template <class OutVector, class InVector>
		OutVector & apply(OutVector & y, const InVector & x) const
		{
			Scalar sum = 0;
			std::vector<PackedUnit>::const_iterator p = _rep.begin();
			typename InVector::const_iterator xp = x.begin();

			for(size_t yi = 0; yi != _rows; ++yi){
				//	for(; xp != x.end(); ++p)
				for(; p != _rep.begin() + ((yi+1) * (_rep.size()/_rows)); ++p){
					for(PackedUnit mask = _mask << (_valPerWord -1); (mask!=0) && (xp!=x.end()); mask>>=1){
						if(mask & (*p))
							_MD.addin(sum, (*xp));
						++xp;
					}
				}
				y[yi] = sum;
				sum = 0;
				xp = x.begin();
			}

			return y;
		}


		//vector-matrix product
		template <class OutVector, class InVector>
		OutVector & applyTranspose(OutVector & y, const InVector & x) const
		{
			Scalar sum = 0;
			std::vector<PackedUnit>::const_iterator p = _rep.begin();
			typename InVector::const_iterator xp = x.begin();
			PackedUnit mask = _mask << (_valPerWord - 1);

			for(size_t yi = 0; yi != _cols; ++yi){
				for(; xp != x.end(); ++xp, p+=(_rep.size()/_rows)){
					if(mask & (*p))
						_MD.addin(sum, (*xp));
				}
				y[yi] = sum;
				sum = 0;
				xp = x.begin();
				p = _rep.begin() + ((yi+1)/_valPerWord);
				mask >>= 1;
				if(mask == 0)
					mask = _mask << (_valPerWord - 1);
			}

			return y;
		}


		//Block apply
		Block & apply(Block & Y, Block & X) const
		{
			//Let i,j represent current spot in BB, used for Block constructions
			std::vector<PackedUnit>::const_iterator p = _rep.begin();
			size_t j = 0;
			_MD.zero(Y);
			Block x(X.field(),1, X.coldim()), y(Y.field(),1, Y.coldim());

			for(size_t i = 0; i != _rows; ++i){
				//! @bug submatrix ?
				Y.subBlock(y,i,0,1,Y.coldim());

				for(; p != _rep.begin() + ((i+1) * (_rep.size()/_rows)); ++p){
					for(PackedUnit mask = _mask << (_valPerWord - 1); (mask!=0) && (j != _cols); mask>>=1, ++j){
						if(mask & (*p))
							_MD.addin(y, X.subBlock(x,j,0,1,X.coldim()));
					}
				}
				j = 0;
			}

			return Y;
		}


		//Block applyTranspose
		Block & applyTranspose(Block & Y, Block & X) const
		{
			//Traverse the BB by going through each column.
			//When BB[i,j] = 1, add the ith column of X
			//to the jth column of Y

			std::vector<PackedUnit>::const_iterator p = _rep.begin();
			PackedUnit mask = _mask << (_valPerWord - 1);
			_MD.zero(Y);
			Block x(X.rowdim(),1), y(Y.rowdim(),1);

			for(size_t j = 0; j != _cols; ++j){
				Y.subBlock (y,0,j,Y.rowdim(),1);
				for(size_t i = 0; i != _rows; ++i, p+=(_rep.size()/_rows)){
					if(mask & (*p))
						_MD.addin(y, X.subBlock(x,0,i,X.rowdim(),1));
				}
				p = _rep.begin() + ((j+1)/_valPerWord);
				mask >>= 1;
				if(mask == 0)
					mask = _mask << (_valPerWord - 1);
			}

			return Y;
		}


		//Unpacking apply
		Block & unpackingApply(Block & Out, Block & In, const size_t U = 1024)
		{
			_MD.zero(Out);
			Block A, B, C;
			size_t iend, jend, kend;

			for(size_t i = 0; i < _rows; i += U){
				if(i + U >= _rows) iend = _rows - i;
				else iend = U;

				for(size_t j = 0; j < In.coldim(); j += U){
					if (j + U >= In.coldim()) jend = In.coldim() - j;
					else jend = U;

					Out.subBlock(C,i,j,iend,jend);
					for(size_t k = 0; k < _cols; k += U){
						if(k + U >= _cols) kend = _cols - k;
						else kend = U;

						expandBlock(A,i,k,iend,kend);
						In.subBlock(B,k,j,kend,jend);
						_MD.axpyin(C, A, B);
					}
				}
			}

			return Out;
		}


		//Unpacking applyTranspose
		Block & unpackingApplyTranspose(Block & Out, Block & In, const size_t U = 1024)
		{
			_MD.zero(Out);
			Block A, B, C;
			size_t iend, jend, kend;

			for(size_t i = 0; i < In.rowdim(); i += U){
				if(i + U >= In.rowdim()) iend = In.rowdim() - i;
				else iend = U;

				for(size_t j = 0; j < _cols; j += U){
					if (j + U >= _cols) jend = _cols - j;
					else jend = U;

					Out.subBlock(C,i,j,iend,jend);
					for(size_t k = 0; k < _rows; k += U){
						if(k + U >= _rows) kend = _rows - k;
						else kend = U;

						expandBlock(A,k,j,kend,jend);
						In.subBlock(B,i,k,iend,kend);
						_MD.axpyin(C, B, A);
					}
				}
			}

			return Out;
		}


		size_t rowdim() const {return _rows;}

		size_t coldim() const {return _cols;}

		MatrixDomain domain() const {return _MD;}

		Field field() const {return _MD;}

		/* Note that if the user calls the following function outside of the range of the
		 * blackbox, the block returned is padded with zeros for out of range values of the blackbox.
		 */
		Block & expandBlock(Block & Out, const size_t i, const size_t j, const size_t numRows, const size_t numCols)
		{
			//Check size, and resize if necessary
			if ( (Out.rowdim() != numRows) || (Out.coldim() != numCols) )
				Out.resize(numRows, numCols);

			//Index values
			size_t unitNum = i*(_rep.size()/_rows) + (j/_valPerWord);  //Starts from 0
			size_t bitPlace = (j%_valPerWord) + 1;                    //Starts from 1

			std::vector<PackedUnit>::const_iterator p;
			size_t outRow = 0;
			size_t outCol = 0;
			Scalar x;
			size_t bp = bitPlace;

			for(;outRow != numRows; ++outRow){
				p = _rep.begin() + unitNum + outRow*(_rep.size()/_rows);
				if(i + outRow >= _rows){
					for(; outCol != numCols; ++outCol)
						Out.setEntry(outRow, outCol, _MD._zro);
				}

				else{
					for(; outCol != numCols; ++p){
						for(; (bp != _valPerWord + 1) && (outCol != numCols); ++bp){
							if (j + outCol >= _cols) x = _MD._zro;
							else _MD.init(x, ((*p) >> (_valPerWord - bp)) & _mask);
							Out.setEntry(outRow, outCol, x);
							++outCol;
						}
						bp = 1;
					}
				}

				bp = bitPlace;
				outCol = 0;
			}

			return Out;
		}

		std::ostream & write(std::ostream & os, Block & B)
		{
			Scalar x;
			for(size_t i = 0; i!= B.rowdim(); ++i){
				for(size_t j = 0; j!= B.coldim(); ++j)
					os << B.getEntry(x,i,j) << ' ';
				os << std::endl;
			}

			return os;
		}

		std::ostream & write(std::ostream & os)
		{
			Scalar x;
			for(size_t i = 0; i != _rows; ++i){
				for(size_t j = 0; j!= _cols; ++j){
					os << getEntry(x,i,j) << ' ';
				}
				os << std::endl;
			}
			return os;
		}


		/** @todo Functions that should be written:
		 * read()
		 * rebind()
		 */

	private:

		//Field which blackbox is over
		MatrixDomain _MD;

		//Number of rows and columns of the matrix
		size_t _rows, _cols;

		//Number of bits in a PackedUnit
		size_t _valPerWord;

		//Representation of the matrix
		std::vector<PackedUnit> _rep;

		//Mask used
		PackedUnit _mask;



	}; // class Dense01Blackbox

} // namespace LinBox

#endif // __LINBOX_dense_zero_one_H



// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

