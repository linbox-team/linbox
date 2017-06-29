#ifndef __compressed_word_h__
#define __compressed_word_h__

#include <functional>
#include <algorithm>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include <numeric>
#include <limits.h>
#include <cstddef>

template <typename T, uint64_t N, uint64_t R>
class CompressedWord
{
	public: // typedefs
		using Base = T;
	
	public: // static helpers
		static uint64_t entries_per_base(); // CHAR_BIT*sizeof(N) / R
		static uint64_t entries_per_word(); // N * entries_per_word()
		static uint64_t num_bases();        // N
		static uint64_t entry_length();     // R
		static uint64_t base_length();      // CHAR_BIT * entries_per_base()
		static Base entry_mask();               // R-bit mask
		static Base base_select_mask();         // 1 at first bit of each entry
		static Base entries_mask();             // base_length()-bit mask
		static std::pair<uint64_t, uint64_t> make_index(uint64_t idx);
	
	public: // data
		std::array<Base,N> b;
		
	public: // constructors
		CompressedWord  () = default;
		CompressedWord  (const CompressedWord &other) = default;
		~CompressedWord () = default;
		CompressedWord  (const CompressedWord &w1, const CompressedWord &w2, const CompressedWord &m);
		CompressedWord  (uint64_t s, uint64_t l);
	
	public: // operators (operate base-by-base)
		CompressedWord& operator =  (const CompressedWord &rhs) = default;
	
		CompressedWord  operator +  (const CompressedWord &rhs) const;
		CompressedWord  operator -  (const CompressedWord &rhs) const;
		CompressedWord  operator *  (const CompressedWord &rhs) const;
		CompressedWord  operator /  (const CompressedWord &rhs) const;
		CompressedWord  operator ^  (const CompressedWord &rhs) const;
		CompressedWord  operator &  (const CompressedWord &rhs) const;
		CompressedWord  operator |  (const CompressedWord &rhs) const;
		CompressedWord& operator += (const CompressedWord &rhs);
		CompressedWord& operator -= (const CompressedWord &rhs);
		CompressedWord& operator *= (const CompressedWord &rhs);
		CompressedWord& operator /= (const CompressedWord &rhs);
		CompressedWord& operator ^= (const CompressedWord &rhs);
		CompressedWord& operator &= (const CompressedWord &rhs);
		CompressedWord& operator |= (const CompressedWord &rhs);
		
		CompressedWord  operator *   (const Base rhs) const;
		CompressedWord  operator /   (const Base rhs) const;
		CompressedWord  operator <<  (const Base rhs) const; // rhs = number of elements
		CompressedWord  operator >>  (const Base rhs) const;
		CompressedWord  operator &   (const Base rhs) const;
		CompressedWord& operator *=  (const Base rhs);
		CompressedWord& operator /=  (const Base rhs);
		CompressedWord& operator <<= (const Base rhs);
		CompressedWord& operator >>= (const Base rhs);
		CompressedWord& operator &=  (const Base rhs);
		
		CompressedWord  operator ~ () const;
		
		bool operator == (const CompressedWord &rhs) const;
		bool operator != (const CompressedWord &rhs) const;
	
	public: // operations
		void clear  ();
		bool isZero () const;
		
		uint64_t count (const Base e) const;
		
		CompressedWord& negin  ();
		
		CompressedWord  axpy   (const CompressedWord &rhs, const Base a) const {
			CompressedWord res{};
			for (uint64_t i = 0; i < N; ++i) {
				res.b[i] = b[i] + rhs.b[i] * a;
			}
			return res;
		}
		
		CompressedWord  nand   (const CompressedWord &rhs) const;
		CompressedWord& nandin (const CompressedWord &rhs);
		CompressedWord  nand   (const Base rhs) const;
		CompressedWord& nandin (const Base rhs);
	
		CompressedWord  lshift   (const uint64_t rhs) const; // element-wise, rhs = number of elements
		CompressedWord  rshift   (const uint64_t rhs) const;
		CompressedWord& lshiftin (const uint64_t rhs);
		CompressedWord& rshiftin (const uint64_t rhs);
		
		Base  getEntry (uint64_t idx) const;
		Base& getEntry (Base &base, uint64_t idx) const;
		Base  getEntry (uint64_t bidx, uint64_t idx) const;
		Base& getEntry (Base &base, uint64_t bidx, uint64_t idx) const;
		
		void setEntry (uint64_t idx, const Base base);
		void setEntry (uint64_t bidx, uint64_t idx, const Base base);
		
		CompressedWord  select     (uint64_t idx, uint64_t l) const;
		CompressedWord& selectin   (uint64_t idx, uint64_t l);
		CompressedWord  deselect   (uint64_t idx, uint64_t l) const;
		CompressedWord& deselectin (uint64_t idx, uint64_t l);
		
		void swapBase  (uint64_t i, uint64_t j);
		void swapEntry (uint64_t i, uint64_t j);
	
	public: // debug
		void random ();
		std::ostream& write (std::ostream &os) const;

};

#include "compressed-word.inl"

#endif