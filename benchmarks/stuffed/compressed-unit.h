#ifndef __compressed_unit_h__
#define __compressed_unit_h__

#include <functional>
#include <algorithm>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include <numeric>
#include <limits.h>
#include <cstddef>

#include "compressed-word.h"

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
class CompressedUnit
{
	public: // typedefs
		using Word = CompressedWord<T,N,R+L>;
		using Base = typename Word::Base;
	
	public: // static helpers
		static uint64_t get_r();
		static uint64_t get_l();
		static uint64_t get_d();
		static Base r_mask(); // R bit mask
		static Base l_mask(); // L bit mask
		static Base r_base_mask(); // all Rs in base mask
		static Base l_base_mask(); // all Ls in base mask
		static Base base_mask(); // 1 at first bit of each entry
		static Base base_entry_mask(); // R+L bit mask
		static Base base_entries_mask(); // all R+Ls in base mask
		static Word word_entries_mask(); // all R+Ls in word mask
	
	public: // data
		std::array<Word,D> w;
	
	public: // constructors
		CompressedUnit () = default;
		CompressedUnit (const CompressedUnit &other) = default;
		~CompressedUnit () = default;
		
	
	public: // operators
		CompressedUnit& operator =   (const CompressedUnit &rhs) = default;
	
		CompressedUnit  operator &   (const Word &rhs) const;
		CompressedUnit& operator &=  (const Word &rhs);
		CompressedUnit  operator &   (const Base rhs) const;
		CompressedUnit& operator &=  (const Base rhs);
		
		CompressedUnit  operator |   (const CompressedUnit &rhs) const;
		CompressedUnit& operator |=  (const CompressedUnit &rhs);
		
		CompressedUnit  operator <<  (const uint64_t rhs) const;
		CompressedUnit& operator <<= (const uint64_t rhs);
		CompressedUnit  operator >>  (const uint64_t rhs) const;
		CompressedUnit& operator >>= (const uint64_t rhs);
		
		CompressedUnit  operator ~   () const;
		
		bool operator == (const CompressedUnit &rhs) const;
		bool operator != (const CompressedUnit &rhs) const;
	
	public: // operations
		void clear  ();
		bool isZero () const;
						
		CompressedUnit& negin  ();				
						
		CompressedUnit  nand   (const Word &rhs) const;
		CompressedUnit& nandin (const Word &rhs);
		CompressedUnit  nand   (const Base rhs) const;
		CompressedUnit& nandin (const Base rhs);
	
		CompressedUnit  lshift   (const uint64_t rhs) const; // element-wise
		CompressedUnit  rshift   (const uint64_t rhs) const;
		CompressedUnit& lshiftin (const uint64_t rhs);
		CompressedUnit& rshiftin (const uint64_t rhs);
						
		CompressedUnit  select     (uint64_t idx, uint64_t l) const;
		CompressedUnit& selectin   (uint64_t idx, uint64_t l);
		CompressedUnit  deselect   (uint64_t idx, uint64_t l) const;
		CompressedUnit& deselectin (uint64_t idx, uint64_t l);
		
		void swapEntry (uint64_t i, uint64_t j);
		void swapWord  (uint64_t i, uint64_t j);
		
		Base getEntry  (uint64_t word, uint64_t i) const;
		Base getEntry  (uint64_t word, uint64_t base, uint64_t i) const;
		
		void setEntry (uint64_t word, uint64_t i, const Base b);
	
	public: // debug
		void random ();
		std::ostream& write (std::ostream &os) const;
	
};

#include "compressed-unit.inl"

#endif