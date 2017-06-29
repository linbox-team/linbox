#ifndef __compressed_unit_inl__
#define __compressed_unit_inl__

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
uint64_t 
CompressedUnit<T,N,R,L,D>::get_r()
{
	return R;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
uint64_t 
CompressedUnit<T,N,R,L,D>::get_l()
{
	return L;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
uint64_t 
CompressedUnit<T,N,R,L,D>::get_d()
{
	return D;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::r_mask() // R bit mask
{
	return (static_cast<typename CompressedUnit<T,N,R,L,D>::Base>(1) << R) - 1;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::l_mask() // L bit mask
{
	return (static_cast<typename CompressedUnit<T,N,R,L,D>::Base>(1) << L) - 1;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::r_base_mask() // all Rs in base mask
{
	typename CompressedUnit<T,N,R,L,D>::Base res{};
	for (uint64_t i = 0; i < CompressedUnit<T,N,R,L,D>::Word::entries_per_base(); ++i) {
		res |= (CompressedUnit<T,N,R,L,D>::r_mask()) << (i*(R+L));
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::l_base_mask() // all Ls in base mask
{
	typename CompressedUnit<T,N,R,L,D>::Base res{};
	for (uint64_t i = 0; i < CompressedUnit<T,N,R,L,D>::Word::entries_per_base(); ++i) {
		res |= (CompressedUnit<T,N,R,L,D>::l_mask()) << (i*(R+L));
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::base_mask() // 1 at first bit of each entry
{
	typename CompressedUnit<T,N,R,L,D>::Base res{};
	for (uint64_t i = 0; i < CompressedUnit<T,N,R,L,D>::Word::entries_per_base(); ++i) {
		res |= (static_cast<typename CompressedUnit<T,N,R,L,D>::Base>(1)) << (i*(R+L));
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::base_entry_mask() // R+L bit mask
{ 
	return CompressedUnit<T,N,R,L,D>::Word::entry_mask();
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::base_entries_mask() // all R+Ls in base mask
{
	typename CompressedUnit<T,N,R,L,D>::Base res{};
	for (uint64_t i = 0; i < CompressedUnit<T,N,R,L,D>::Word::entries_per_base(); ++i) {
		res |= (CompressedUnit<T,N,R,L,D>::Word::entry_mask()) << (i*(R+L));
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Word 
CompressedUnit<T,N,R,L,D>::word_entries_mask() // all R+Ls in word mask
{
	typename CompressedUnit<T,N,R,L,D>::Word res{};
	for (uint64_t i = 0; i < N; ++i) {
		res.b[i] = CompressedUnit<T,N,R,L,D>::base_entries_mask();
	}
	return res;
}
	
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::operator &   (const typename CompressedUnit<T,N,R,L,D>::Word &rhs) const
{
	CompressedUnit<T,N,R,L,D> res{};
	//auto l = [rhs](auto a) { return a & rhs; };
	//std::transform(w.begin(), w.end(), res.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i] & rhs.w[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::operator &=  (const typename CompressedUnit<T,N,R,L,D>::Word &rhs)
{
	return *this = *this & rhs;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::operator &   (const typename CompressedUnit<T,N,R,L,D>::Base rhs) const
{
	CompressedUnit<T,N,R,L,D> res{};
	//auto l = [rhs](auto a) { return a & rhs; };
	//std::transform(w.begin(), w.end(), res.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i] & rhs;
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::operator &=  (const typename CompressedUnit<T,N,R,L,D>::Base rhs)
{
	return *this = *this & rhs;
}
		
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::operator |   (const CompressedUnit<T,N,R,L,D> &rhs) const
{
	CompressedUnit<T,N,R,L,D> res{};
	//std::transform(w.begin(), w.end(), rhs.w.begin(), res.w.begin(), std::bit_or<typename CompressedUnit<T,N,R,L,D>::Word>());
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i] | rhs.w[i];
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::operator |=  (const CompressedUnit<T,N,R,L,D> &rhs)
{
	return *this = *this | rhs;
}
		
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::operator <<  (const uint64_t rhs) const
{
	CompressedUnit<T,N,R,L,D> res{};
	//auto l = [rhs](auto a) { return a << rhs; };
	//std::transform(w.begin(), w.end(), res.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i] << rhs;
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::operator <<= (const uint64_t rhs)
{
	return *this = *this << rhs;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::operator >>  (const uint64_t rhs) const
{
	CompressedUnit<T,N,R,L,D> res{};
	//auto l = [rhs](auto a) { return a >> rhs; };
	//std::transform(w.begin(), w.end(), res.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i] >> rhs;
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::operator >>= (const uint64_t rhs)
{
	return *this = *this >> rhs;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline		
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::operator ~   () const
{
	CompressedUnit<T,N,R,L,D> res{};
	//auto l = [](auto a) { return ~a; };
	//std::transform(w.begin(), w.end(), res.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = ~w[i];
	}
	return res;
}
		
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
bool 
CompressedUnit<T,N,R,L,D>::operator == (const CompressedUnit<T,N,R,L,D> &rhs) const
{
	//return std::inner_product(w.begin(), w.end(), rhs.w.begin(), true, std::bit_and<typename CompressedUnit<T,N,R,L,D>::Word>(), std::equal_to<typename CompressedUnit<T,N,R,L,D>::Word>());
	bool res = true;
	for (uint64_t i = 0; i < D; ++i) {
		if (w[i] != rhs.w[i]) {
			res = false;
		}
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
bool 
CompressedUnit<T,N,R,L,D>::operator != (const CompressedUnit<T,N,R,L,D> &rhs) const
{
	return !(*this == rhs);
}
	
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
void 
CompressedUnit<T,N,R,L,D>::clear  ()
{
	//auto l = [](auto &a) { a.clear(); };
	//std::for_each(w.begin(), w.end(), l);
	for (uint64_t i = 0; i < D; ++i) {
		w[i].clear();
	}
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
bool 
CompressedUnit<T,N,R,L,D>::isZero () const
{
	//auto l = [](auto a, auto b) { return a.isZero() && b; };
	//return std::accumulate(w.begin(), w.end(), true, l);
	
	bool res = true;
	for (uint64_t i = 0; i < D; ++i) {
		if (w[i].isZero() == false) {
			res = false;
		}
	}
	return res;
}
						
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>&
CompressedUnit<T,N,R,L,D>::negin  ()
{
	return *this = ~(*this);
}				
						
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::nand   (const typename CompressedUnit<T,N,R,L,D>::Word &rhs) const
{
	//auto l = [](auto a, auto b) { return a.nand(b); };
	CompressedUnit<T,N,R,L,D> res{};
	//std::transform(w.begin(), w.end(), rhs.w.begin(), res.w.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i].nand(rhs);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::nandin (const typename CompressedUnit<T,N,R,L,D>::Word &rhs)
{
	return *this = nand(rhs);
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::nand   (const typename CompressedUnit<T,N,R,L,D>::Base rhs) const
{
	//auto l = [](auto a, auto b) { return a.nand(b); };
	CompressedUnit<T,N,R,L,D> res{};
	//std::transform(w.begin(), w.end(), rhs.w.begin(), res.w.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i].nand(rhs);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::nandin (const typename CompressedUnit<T,N,R,L,D>::Base rhs)
{
	return *this = nand(rhs);
}
	
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::lshift   (const uint64_t rhs) const // element-wise
{
	//auto l = [rhs](auto a) { return a.lshift(rhs); };
	CompressedUnit<T,N,R,L,D> res{};
	//std::transform(w.begin(), w.end(), rhs.w.begin(), res.w.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i].lshift(rhs);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::rshift   (const uint64_t rhs) const
{
	//auto l = [rhs](auto a) { return a.rshift(rhs); };
	CompressedUnit<T,N,R,L,D> res{};
	//std::transform(w.begin(), w.end(), rhs.w.begin(), res.w.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i].rshift(rhs);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::lshiftin (const uint64_t rhs)
{
	return *this = lshift(rhs);
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::rshiftin (const uint64_t rhs)
{
	return *this = rshift(rhs);
}
						
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::select     (uint64_t idx, uint64_t l) const
{
	//auto l = [idx, l](auto a) { return a.select(idx, l); };
	CompressedUnit<T,N,R,L,D> res{};
	//std::transform(w.begin(), w.end(), rhs.w.begin(), res.w.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i].select(idx,l);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::selectin   (uint64_t idx, uint64_t l)
{
	return *this = select(idx, l);
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>  
CompressedUnit<T,N,R,L,D>::deselect   (uint64_t idx, uint64_t l) const
{
	//auto l = [idx, l](auto a) { return a.deselect(idx, l); };
	CompressedUnit<T,N,R,L,D> res{};
	//std::transform(w.begin(), w.end(), rhs.w.begin(), res.w.begin(), l);
	for (uint64_t i = 0; i < D; ++i) {
		res.w[i] = w[i].deselect(idx,l);
	}
	return res;
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
CompressedUnit<T,N,R,L,D>& 
CompressedUnit<T,N,R,L,D>::deselectin (uint64_t idx, uint64_t l)
{
	return *this = deselect(idx, l);
}
		
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
void 
CompressedUnit<T,N,R,L,D>::swapEntry (uint64_t i, uint64_t j)
{
	//auto l = [i, j](auto &a) { a.swapEntry(i, j); };
	//std::for_each(w.begin(), w.end(), l);
	for (uint64_t i = 0; i < D; ++i) {
		w[i].swapEntry(i,j);
	}
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
void 
CompressedUnit<T,N,R,L,D>::swapWord (uint64_t i, uint64_t j)
{
	std::swap(w[i], w[j]);
}
		
template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::getEntry  (uint64_t word, uint64_t i) const
{
	return w[word].getEntry(i);
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedUnit<T,N,R,L,D>::Base 
CompressedUnit<T,N,R,L,D>::getEntry  (uint64_t word, uint64_t base, uint64_t i) const
{
	return w[word].getEntry(base,i);
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
void
CompressedUnit<T,N,R,L,D>::setEntry  (uint64_t word, uint64_t i, const typename CompressedUnit<T,N,R,L,D>::Base b)
{
	w[word].setEntry(i, b);
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
void 
CompressedUnit<T,N,R,L,D>::random ()
{
	//auto l = [](auto &a) { a.random(); };
	//std::for_each(w.begin(), w.end(), l);
	for (uint64_t i = 0; i < D; ++i) {
		w[i].random();
	}
}

template <typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
std::ostream& 
CompressedUnit<T,N,R,L,D>::write (std::ostream &os) const
{	
	for (uint64_t j = 0; j < N; ++j) {
		for (uint64_t k = 0; k < CompressedUnit<T,N,R,L,D>::Word::entries_per_base(); ++k) {
			os << "(";
			for (uint64_t i = 0; i < D; ++i) {
				os << w[i].getEntry(j,k);
			}
			os << ")";
		}
	}
	return os;
}

#endif