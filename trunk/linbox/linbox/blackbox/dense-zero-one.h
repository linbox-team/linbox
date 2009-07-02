/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/dense-zero-one.h
 * 
 * Written by Nick Messina <nmessina@mail.eecis.udel.edu>
 *
 * -----------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __DENSE_ZERO_ONE_H
#define __DENSE_ZERO_ONE_H



#include <iostream>
//#include <linbox/field/hom.h>
#include "linbox/linbox-config.h"
#include <vector>
#include <linbox/blackbox/blackbox-interface.h>


const size_t SIZE = sizeof(LinBox::uint64) * 8;


namespace LinBox
{

	template <class _Field>
	class DenseZeroOne : public  BlackboxInterface 
	{

	    public:

		typedef _Field Field;
		typedef typename Field::Element Element;
		typedef LinBox::uint64 PackedUnit;
		typedef DenseZeroOne<_Field> Self_t;
	
		/* Notes:
		 * _rows gives the number of rows of PackedUnits
		 * _cols gives the number of bits in each row */

		DenseZeroOne(Field F, size_t m, size_t n) 
		: _F(F), _rows(m), _cols(n)
		{
			_valPerWord = SIZE;
			_mask = 1;
			_rep.resize( _rows * ( ((_cols-1)/_valPerWord) + 1));
		}
	
		DenseZeroOne(const Self_t & other)
		: _F(other._F), _rows(other._rows), _cols(other._cols), _valPerWord(other._valPerWord), _rep(other._rep), _mask(other._mask)
		{}

		/* Constructor for making the blackbox out of a vector of PackedUnits.
		 * Note that if the vector given is too short for given size(m by n),
		 * the blackbox is enlarged and padded by 0's added after the given vector.
		 * If the vector given is too large, the blackbox truncates this vector
		 * to the appropriate number of PackedUnits.  */
		template <class Vector>
		DenseZeroOne(Field F, size_t m, size_t n, Vector & v)
		: _F(F), _rows(m), _cols(n)
		{
			_valPerWord = SIZE;
			_mask = 1;
			_rep = v;
			_rep.resize( _rows * ( ((_cols-1)/_valPerWord) + 1));
		}

	
		Element & getEntry(Element & x, size_t i, size_t j) const
		{
			//Get all index values needed.
			size_t unitNum = i*(_rep.size()/_rows) + (j/_valPerWord);  //Starts from 0
			size_t bitPlace = ((j+1)%_valPerWord);                    //Starts from 1

			_F.init(x, (_rep[unitNum] >> (_valPerWord - bitPlace) ) & _mask );
			return x;
		}
	
		void setEntry(size_t i, size_t j, PackedUnit & x)
		{
			//Get all index values needed.
			size_t unitNum = i*(_rep.size()/_rows) + (j/_valPerWord);  //Starts from 0
			size_t bitPlace = ((j+1)%_valPerWord);                    //Starts from 1

			_rep[unitNum] &= ~(_mask << (_valPerWord - bitPlace));
			_rep[unitNum] |= (x << (_valPerWord - bitPlace)); 
		}

	//Apply Variants
	
		//Speedup:
		//instead of using getEntry, walk through each PackedUnit
	
		//apply, y := Ax, matrix-vector product
		//should there be both handled and handleless apply?

		//handleless apply
		template <class OutVector, class InVector>
		OutVector & apply(OutVector & y, const InVector & x) const 
		{
		 	Element sum = 0;
			Element a;
			for(size_t i = 0; i != _rows; ++i){
				for(size_t j = 0; j != _cols; ++j){
					_F.axpyin(sum, getEntry(a,i,j), x[j]); 
				}
				y[i] = sum;
				sum = 0;
			}
		
			return y;
		}


		//applyTranspose, y := xA, vector-matrix product
		//should there be both handled and handleless apply?

		//handleless applyTranspose
		template <class OutVector, class InVector>
		OutVector & applyTranspose(OutVector & y, const InVector & x) const 
		{
		 	Element sum = 0;
			Element a;
			for(size_t j = 0; j != _cols; ++j){
				for(size_t i = 0; i != _rows; ++i){
					if(!_F.isZero(getEntry(a,i,j)))
						_F.addin(sum, x[i]);	
				}
				y[j] = sum;
				sum = 0;
			}

			return y;
		}		
		
		
		size_t rowdim() const {return _rows;}

		size_t coldim() const {return _cols;}

		Field field() const {return _F;}


		std::ostream & write(std::ostream & os)
		{		
			Element x;
			for(size_t i = 0; i != _rows; ++i){
				for(size_t j = 0; j!= _cols; ++j){
					os << getEntry(x,i,j) << ' ';
				}
				os << std::endl;
			}
			return os;
		}
	
	
	/* Functions that should be written: 
	 * read()
	 * rebind()
   	 */ 

	    private:
		
		//Field which blackbox is over
		Field _F;
		
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

#endif // __DENSE_ZERO_ONE_H

