#ifndef __compressed_word_inl__
#define __compressed_word_inl__

template <typename T, uint64_t N, uint64_t R>
inline
uint64_t 
CompressedWord<T,N,R>::entries_per_base()
{
	return static_cast<uint64_t>(CHAR_BIT*sizeof(T)) / R;
}

template <typename T, uint64_t N, uint64_t R>
inline
uint64_t 
CompressedWord<T,N,R>::entries_per_word()
{
	return CompressedWord<T,N,R>::entries_per_base() * CompressedWord<T,N,R>::num_bases();
}

template <typename T, uint64_t N, uint64_t R>
inline
uint64_t 
CompressedWord<T,N,R>::num_bases()
{
	return N;
}

template <typename T, uint64_t N, uint64_t R>
inline
uint64_t 
CompressedWord<T,N,R>::entry_length()
{
	return R;
}

template <typename T, uint64_t N, uint64_t R>
inline
uint64_t 
CompressedWord<T,N,R>::base_length()
{
	return CHAR_BIT * CompressedWord<T,N,R>::entries_per_base();
}

template <typename T, uint64_t N, uint64_t R>
inline
typename CompressedWord<T,N,R>::Base 
CompressedWord<T,N,R>::entry_mask()
{
	return (static_cast<typename CompressedWord<T,N,R>::Base>(-1) >> (CHAR_BIT*sizeof(T) - R));
}

template <typename T, uint64_t N, uint64_t R>
inline
typename CompressedWord<T,N,R>::Base 
CompressedWord<T,N,R>::base_select_mask()
{
	const auto b1 = static_cast<typename CompressedWord<T,N,R>::Base>(1);
	CompressedWord<T,N,R>::Base res{};
	for (uint64_t i = 0; i < base_length(); i += CompressedWord<T,N,R>::entry_length()) {
		res |= (b1 << i);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
typename CompressedWord<T,N,R>::Base 
CompressedWord<T,N,R>::entries_mask()
{
	auto res = ~static_cast<typename CompressedWord<T,N,R>::Base>(0);
	res = res >> (CHAR_BIT*sizeof(T) - CompressedWord<T,N,R>::base_length());
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
std::pair<uint64_t, uint64_t> 
CompressedWord<T,N,R>::make_index(uint64_t idx)
{
	return std::make_pair(idx / CompressedWord<T,N,R>::entries_per_base(), R*(idx % CompressedWord<T,N,R>::entries_per_base()));
}

template <typename T, uint64_t N, uint64_t R>
CompressedWord<T,N,R>::CompressedWord(const CompressedWord<T,N,R> &w1, 
										const CompressedWord<T,N,R> &w2, 
										const CompressedWord<T,N,R> &m)
	: b{}
{
	for (uint64_t i = 0; i < N; ++i) {
		b[i] = w2.b[i] ^ ((w1.b[i] ^ w2.b[i]) & m.b[i]);
	}
}

template <typename T, uint64_t N, uint64_t R>
CompressedWord<T,N,R>::CompressedWord(uint64_t s, uint64_t l)
	: b{}
{
	const auto idxs = CompressedWord<T,N,R>::make_index(s);
	const auto idxe = CompressedWord<T,N,R>::make_index(s+l);
	
	if (idxs.first == idxe.first) {
		//one word
		b[idxs.first] = ((static_cast<typename CompressedWord<T,N,R>::Base>(1) << (R*l)) - 1) << idxs.first;
	}
	else {
		//first
		b[idxs.first] = ~((static_cast<typename CompressedWord<T,N,R>::Base>(1) << idxs.second) - 1);
		
		//middle
		for (uint64_t i = idxs.first + 1; i < idxe.first; ++i) {
			b[i] = ~static_cast<typename CompressedWord<T,N,R>::Base>(0);
		}
		
		//last
		b[idxe.first] = (static_cast<typename CompressedWord<T,N,R>::Base>(1) << idxe.second) - 1;
	}
}
	
template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator +  (const CompressedWord<T,N,R> &rhs) const
{
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), std::plus<typename CompressedWord<T,N,R>::Base>());
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] + rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator -  (const CompressedWord<T,N,R> &rhs) const
{
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), std::minus<typename CompressedWord<T,N,R>::Base>());
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] - rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator *  (const CompressedWord<T,N,R> &rhs) const
{
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), std::multiplies<typename CompressedWord<T,N,R>::Base>());
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] * rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator /  (const CompressedWord<T,N,R> &rhs) const
{
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), std::divides<typename CompressedWord<T,N,R>::Base>());
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] / rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator ^  (const CompressedWord<T,N,R> &rhs) const
{
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), std::bit_xor<typename CompressedWord<T,N,R>::Base>());
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] ^ rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator &  (const CompressedWord<T,N,R> &rhs) const
{
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), std::bit_and<typename CompressedWord<T,N,R>::Base>());
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] & rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R> 
CompressedWord<T,N,R>::operator |  (const CompressedWord<T,N,R> &rhs) const
{
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), std::bit_or<typename CompressedWord<T,N,R>::Base>());
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] | rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator += (const CompressedWord<T,N,R> &rhs)
{
	return *this = *this + rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator -= (const CompressedWord<T,N,R> &rhs)
{
	return *this = *this - rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator *= (const CompressedWord<T,N,R> &rhs)
{
	return *this = *this * rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator /= (const CompressedWord<T,N,R> &rhs)
{
	return *this = *this / rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator ^= (const CompressedWord<T,N,R> &rhs)
{
	return *this = *this ^ rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator &= (const CompressedWord<T,N,R> &rhs)
{
	return *this = *this & rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator |= (const CompressedWord<T,N,R> &rhs)
{
	return *this = *this | rhs;
}
		
template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator *   (const typename CompressedWord<T,N,R>::Base rhs) const
{
	//auto l = [rhs](auto a) { return a * rhs; };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), res.begin(), l);
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] * rhs;
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator /   (const typename CompressedWord<T,N,R>::Base rhs) const
{
	//auto l = [rhs](auto a) { return a / rhs; };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), res.begin(), l);
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] / rhs;
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator <<  (const typename CompressedWord<T,N,R>::Base rhs) const
{
	//auto l = [rhs,R](auto a) { return a << (R*rhs); };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), res.begin(), l);
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] << (R*rhs);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator >>  (const typename CompressedWord<T,N,R>::Base rhs) const
{
	//auto l = [rhs,R](auto a) { return a >> (R*rhs); };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), res.begin(), l);
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] >> (R*rhs);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::operator &   (const typename CompressedWord<T,N,R>::Base rhs) const
{
	//auto l = [rhs](auto a) { return a & rhs; };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.begin(), res.begin(), l);
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] & rhs;
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator *=  (const typename CompressedWord<T,N,R>::Base rhs)
{
	return *this = *this * rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator /=  (const typename CompressedWord<T,N,R>::Base rhs)
{
	return *this = *this / rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator <<= (const typename CompressedWord<T,N,R>::Base rhs)
{
	return *this = *this << rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator >>= (const typename CompressedWord<T,N,R>::Base rhs)
{
	return *this = *this >> rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::operator &=  (const typename CompressedWord<T,N,R>::Base rhs)
{
	return *this = *this & rhs;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>
CompressedWord<T,N,R>::operator ~  () const
{
	//auto l = [](auto a) { return ~a; };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), res.b.begin(), l)
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = ~b[i];
	}
	return res;
}
		
template <typename T, uint64_t N, uint64_t R>
inline
bool 
CompressedWord<T,N,R>::operator == (const CompressedWord<T,N,R> &rhs) const
{
	//return std::inner_product(b.begin(), b.end(), rhs.b.begin(), true, std::logical_and<typename CompressedWord<T,N,R>::Base>(), std::equal_to<typename CompressedWord<T,N,R>::Base>());
	bool res = true;
	for (uint64_t i = 0; i < N; ++i) {
		if (b[i] != rhs.b[i]) {
			res = false;
		}
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
bool 
CompressedWord<T,N,R>::operator != (const CompressedWord<T,N,R> &rhs) const
{
	return !(*this == rhs);
}
	
template <typename T, uint64_t N, uint64_t R>
inline
void 
CompressedWord<T,N,R>::clear  ()
{
	b.fill(static_cast<typename CompressedWord<T,N,R>::Base>(0));
}

template <typename T, uint64_t N, uint64_t R>
inline
bool 
CompressedWord<T,N,R>::isZero () const
{
	//auto l = [](auto a, auto res) { return ((a & CompressedWord<T,N,R>::entries_mask()) == 0) && res; };
	//return std::accumulate(b.begin(), b.end(), true, l);
	bool res = true;
	for (uint64_t i = 0; i < N; ++i) {
		if (b[i] != 0) {
			res = false;
		}
	}
	return res;
}
		
template <typename T, uint64_t N, uint64_t R>
inline
uint64_t 
CompressedWord<T,N,R>::count (const typename CompressedWord<T,N,R>::Base e) const
{
	uint64_t res = 0;
	for (uint64_t i = 0; i < N; ++i) {
		for (uint64_t j = 0; j < CompressedWord<T,N,R>::base_length(); j += R) {
			if (((b[i] >> j) & CompressedWord<T,N,R>::entry_mask()) == e) {
				++res;
			}
		}
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::negin ()
{
	*this = ~(*this);
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R> 
CompressedWord<T,N,R>::nand (const CompressedWord<T,N,R> &rhs) const
{
	//auto l = [](auto a, auto b) { return a & ~b; };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), l);
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] & ~rhs.b[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::nandin (const CompressedWord<T,N,R> &rhs)
{
	return *this = nand(rhs);
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R> 
CompressedWord<T,N,R>::nand (const typename CompressedWord<T,N,R>::Base rhs) const
{
	//auto l = [rhs](auto a) { return a & ~rhs; };
	CompressedWord<T,N,R> res{};
	//std::transform(b.begin(), b.end(), rhs.b.begin(), res.b.begin(), l);
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = b[i] & ~rhs;
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::nandin (const typename CompressedWord<T,N,R>::Base rhs)
{
	return *this = nand(rhs);
}
	
template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::lshift   (const uint64_t rhs) const // element-wise
{
	const auto off = CompressedWord<T,N,R>::make_index(rhs);
	const uint64_t rs = (-off.second) % (CHAR_BIT*sizeof(T));
	CompressedWord<T,N,R> res{};
	res.b[off.first] = (b[0] & CompressedWord<T,N,R>::entries_mask()) << off.second;
	for (uint64_t i = off.first+1, j = 0; i < N; ++i, ++j) {
		res.b[i] = (b[j] >> rs) | ((b[j+1] & CompressedWord<T,N,R>::entries_mask()) << off.second);
	} 
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::rshift   (const uint64_t rhs) const
{
	const auto off = CompressedWord<T,N,R>::make_index(rhs);
	const uint64_t ls = (-off.second) % (CHAR_BIT*sizeof(T));
	CompressedWord<T,N,R> res{};
	for (uint64_t i = 0, j = off.first; j < N-1; ++i, ++j) {
		res.b[i] = (b[j] >> off.second) | ((b[j+1] & CompressedWord<T,N,R>::entries_mask()) << ls);
	}
	res.b[off.first] = b[N-1] >> off.second; 
	return res;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::lshiftin (const uint64_t rhs)
{
	return *this = lshift(rhs);
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::rshiftin (const uint64_t rhs)
{
	return *this = rshift(rhs);
}
		
template <typename T, uint64_t N, uint64_t R>
inline
typename CompressedWord<T,N,R>::Base  
CompressedWord<T,N,R>::getEntry (uint64_t idx) const
{
	const auto off = CompressedWord<T,N,R>::make_index(idx);
	return getEntry(off.first, off.second);
}

template <typename T, uint64_t N, uint64_t R>
inline
typename CompressedWord<T,N,R>::Base& 
CompressedWord<T,N,R>::getEntry (typename CompressedWord<T,N,R>::Base &base, uint64_t idx) const
{
	return base = getEntry(idx);
}

template <typename T, uint64_t N, uint64_t R>
inline
typename CompressedWord<T,N,R>::Base  
CompressedWord<T,N,R>::getEntry (uint64_t bidx, uint64_t idx) const
{
	return (b[bidx] >> idx) & CompressedWord<T,N,R>::entry_mask();
}

template <typename T, uint64_t N, uint64_t R>
inline
typename CompressedWord<T,N,R>::Base& 
CompressedWord<T,N,R>::getEntry (typename CompressedWord<T,N,R>::Base &base, uint64_t bidx, uint64_t idx) const
{
	return base = getEntry(bidx, idx);
}
		
template <typename T, uint64_t N, uint64_t R>
inline
void 
CompressedWord<T,N,R>::setEntry (uint64_t idx, const typename CompressedWord<T,N,R>::Base base)
{
	const auto off = CompressedWord<T,N,R>::make_index(idx);
	setEntry(off.first, off.second, base);
}

template <typename T, uint64_t N, uint64_t R>
inline
void 
CompressedWord<T,N,R>::setEntry (uint64_t bidx, uint64_t idx, const typename CompressedWord<T,N,R>::Base base)
{
	b[bidx] &= ~(CompressedWord<T,N,R>::entry_mask() << idx);
	b[bidx] |=  (base << idx);
}
		
template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::select     (uint64_t idx, uint64_t l) const
{
	const CompressedWord<T,N,R> m(idx,l);
	return *this & m;
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::selectin   (uint64_t idx, uint64_t l)
{
	return *this = select(idx, l);
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>  
CompressedWord<T,N,R>::deselect   (uint64_t idx, uint64_t l) const
{
	const CompressedWord<T,N,R> m(idx,l);
	return nand(m);
}

template <typename T, uint64_t N, uint64_t R>
inline
CompressedWord<T,N,R>& 
CompressedWord<T,N,R>::deselectin (uint64_t idx, uint64_t l)
{
	return *this = deselect(idx, l);
}
		
template <typename T, uint64_t N, uint64_t R>
inline
void 
CompressedWord<T,N,R>::swapBase  (uint64_t i, uint64_t j)
{
	std::swap(b[i], b[j]);
}

template <typename T, uint64_t N, uint64_t R>
inline
void 
CompressedWord<T,N,R>::swapEntry (uint64_t i, uint64_t j)
{
	const auto offi = CompressedWord<T,N,R>::make_index(i);
	const auto offj = CompressedWord<T,N,R>::make_index(j);
	
	const auto ei = ((b[offi.first] >> offi.second) & CompressedWord<T,N,R>::entry_mask());
	const auto ej = ((b[offj.first] >> offj.second) & CompressedWord<T,N,R>::entry_mask());
	
	b[offi.first]  &= ~(CompressedWord<T,N,R>::entry_mask() << offi.second);
	b[offj.second] &= ~(CompressedWord<T,N,R>::entry_mask() << offj.second);
	
	b[offi.first] |= ej << offi.second;
	b[offj.first] |= ei << offj.second;
}
	
template <typename T, uint64_t N, uint64_t R>
inline
void 
CompressedWord<T,N,R>::random ()
{
	for (uint64_t i = 0; i < CompressedWord<T,N,R>::entries_per_base(); ++i) {
		setEntry(i, static_cast<typename CompressedWord<T,N,R>::Base>(rand() % (static_cast<uint64_t>(1) << R)));
	}
}

template <typename T, uint64_t N, uint64_t R>
inline
std::ostream& 
CompressedWord<T,N,R>::write (std::ostream &os) const
{
	for (uint64_t i = 0; i < CompressedWord<T,N,R>::entries_per_base(); ++i) {
		os << getEntry(i);
	}
	return os;
}

#endif