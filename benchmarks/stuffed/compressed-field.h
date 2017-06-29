#ifndef __compressed_field_h__
#define __compressed_field_h__

#include <functional>
#include <algorithm>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include <numeric>
#include <limits.h>
#include <cstddef>

#include <x86intrin.h>

#include "compressed-word.h"
#include "compressed-unit.h"

// P is the prime of the field
// T is storage unit (uint64_t, etc), unsigned
// N is number of storage units used
// R is number of bits required to hold prime
// L is number of overflow bits (for delayed normalization)
// D is the depth of the Unit (1 if packed, > 1 if Sliced)
template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
class CompressedField
{
	public: // typedefs
		using Unit = CompressedUnit<T,N,R,L,D>;
		using Word = typename Unit::Word;
		using Base = typename Unit::Base;
		using Element = Base;
		
	public: // data
		Element one = static_cast<Element>(1);
		Element mOne = static_cast<Element>(P-1);
		Element zero = static_cast<Element>(0);
		Unit    uZero{};
	
	public: // static helpers
		static uint64_t normalize_freq(); // if 0, no normalization needed
		static uint64_t addin_freq(); // if 0, no normalization needed
		static uint64_t axpy_freq(); // if 0, no normalization needed
		static uint64_t m4rm_freq(); // if 0, no normalization needed
		static uint64_t four_russians_depth() { return D*R; }
	
	public: // helpers
		uint64_t characteristic() const { return P; }
		uint64_t cardinality()    const { return P; }
		
		Element rand_entry() { return static_cast<Element>(rand() % P); }
	
	public: // constructor
		CompressedField() = default;
		~CompressedField() = default;

	public: // operations
		Unit  add    (const Unit &lhs, const Unit &rhs) const;
		Unit& addin  (Unit &lhs, const Unit &rhs) const { return lhs = add(lhs, rhs); }
		Unit  sub    (const Unit &lhs, const Unit &rhs) const;
		Unit& subin  (Unit &lhs, const Unit &rhs) const { return lhs = sub(lhs, rhs); }
		Unit  axpy   (const Unit &lhs, const Unit &rhs, const Element a) const;
		Unit& axpyin (Unit &lhs, const Unit &rhs, const Element a) const { return lhs = axpy(lhs, rhs, a); }
		
		Unit  neg    (const Unit &lhs) const;
		Unit& negin  (Unit &lhs) const { return lhs = neg(lhs); }
		Unit  smul   (const Unit &lhs, const Element a) const;
		Unit& smulin (Unit &lhs, const Element a) const { return lhs = smul(lhs, a); }
		
		Element  getEntry  (const Unit &lhs, uint64_t idx) const;
		Unit&    setEntry  (Unit &lhs, const Element e, uint64_t idx) const;
		Unit&    initEntry (Unit &lhs, const Element e, uint64_t idx) const { return setEntry(lhs, e % P, idx); }
		
		Unit normalize  (const Unit &lhs) const;
		Unit pnormalize (const Unit &lhs) const;
		
		Base getEntries (const Unit * lhs, uint64_t idx, uint64_t t, uint64_t d) const; // t entries starting at idx for d-th use in four russians (can span bases and/or words)
		Base getEntries (const Unit &lhs, uint64_t widx, uint64_t idx, uint64_t t, uint64_t d) const; // must be in 1 base (and one word)
		
	public: // debug
		Unit random () const;
		std::ostream& write (std::ostream &os, const Unit &lhs) const;
};

#include "compressed-field.inl"

#endif
