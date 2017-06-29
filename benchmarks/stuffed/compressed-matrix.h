#ifndef __compressed_matrix_h__
#define __compressed_matrix_h__

#include <functional>
#include <algorithm>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include <numeric>
#include <limits.h>
#include <cstddef>

#include "compressed-word.h"
#include "compressed-unit.h"
#include "compressed-field.h"

template <typename Field_>
class CompressedMatrix
{
	public: // typedefs
		using Field     = Field_;
		using Unit      = typename Field::Unit;
		using Word      = typename Field::Word;
		using Base      = typename Field::Base;
		using Element   = typename Field::Element;
		using Matrix    = CompressedMatrix<Field>;
		using Submatrix = Matrix;
	
	public: // data
		uint64_t rows;
		uint64_t cols;
		uint64_t words;
		uint64_t stride;
		uint64_t loff;
		uint64_t roff;
		Field& f;
		bool alloc;
		Unit * rep;
		
	public: // static helpers
		static uint64_t get_words (uint64_t loff, uint64_t sc, uint64_t nc);
		static uint64_t get_loff  (uint64_t loff, uint64_t sc);
		static uint64_t get_roff  (uint64_t loff, uint64_t sc, uint64_t nc);
	
	public: // constructors
		CompressedMatrix  (Field &f_);
		CompressedMatrix  (Field &f_, uint64_t nr, uint64_t nc);
		CompressedMatrix  (const CompressedMatrix &other);
		~CompressedMatrix ();
		CompressedMatrix& operator = (const CompressedMatrix &other);
		CompressedMatrix  (const CompressedMatrix &other, uint64_t sr, uint64_t sc, uint64_t nr, uint64_t nc);
	
	public: // information
		uint64_t rowdim () const;
		uint64_t coldim () const;
		uint64_t getwords  () const;
		uint64_t getstride () const;
		uint64_t getloff   () const;
		uint64_t getroff   () const;
		Field& field ();
	
	public: // operations
		CompressedMatrix& submatrix (const CompressedMatrix &other, uint64_t sr, uint64_t sc, uint64_t nr, uint64_t nc);
		void clear();
		
		Element getEntry (uint64_t ri, uint64_t ci) const;
		void setEntry (unsigned ri, uint64_t ci, const Element eij);
		const Unit * getRow (uint64_t i) const;
		
	public: // debug
		void random();
		std::ostream& write (std::ostream &os = std::cout) const;
};

#include "compressed-matrix.inl"

#endif