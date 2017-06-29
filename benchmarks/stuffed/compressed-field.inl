#ifndef __compressed_field_inl__
#define __compressed_field_inl__

namespace {
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_normalize_freq
	{
		uint64_t operator() () const { return ((static_cast<uint64_t>(1) << (R+L)) - 1) / ((P-1)*(P-1)); }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t D>
	struct CF_normalize_freq<P,T,N,static_cast<uint64_t>(1),0u,D>
	{
		uint64_t operator() () const { return 0; }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_addin_freq
	{
		uint64_t operator() () const { return ((static_cast<uint64_t>(1) << (R+L)) - 1) / (P-1); }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t D>
	struct CF_addin_freq<P,T,N,static_cast<uint64_t>(1),0u,D>
	{
		uint64_t operator() () const { return 0; }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_axpy_freq
	{
		uint64_t operator() () const { return 1 + (((static_cast<uint64_t>(1) << (R+L)) - 1) - (P-1)) / ((P-1)*(P-1)); }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_m4rm_freq
	{
		uint64_t operator() () const { return 1 + (((static_cast<uint64_t>(1) << (R+L)) - 1) - (P-1)) / (((1u << R) -1)*(P-1)); }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t D>
	struct CF_m4rm_freq<P,T,N,static_cast<uint64_t>(1),0u,D>
	{
		uint64_t operator() () const { return 0; }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t D>
	struct CF_axpy_freq<P,T,N,static_cast<uint64_t>(1),0u,D>
	{
		uint64_t operator() () const { return 0; }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_add
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		Unit operator ()(const Unit &lhs, const Unit &rhs) const { 
			Unit res{};
			res.w[0] = lhs.w[0] + rhs.w[0];
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_add<3,T,N,1,0,2>
	{
		using Unit = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Unit;
		Unit operator()(const Unit &lhs, const Unit &rhs) const {
			const auto temp0 = lhs.w[0] ^ rhs.w[1];
  			const auto temp1 = lhs.w[1] ^ rhs.w[0];
			const auto s = temp0 ^ lhs.w[1];
			const auto t = temp1 ^ rhs.w[1];
			const auto r0 = temp0 & temp1;
			const auto r1 = s | t;
			Unit res{};
			res.w[0] = r0;
			res.w[1] = r1;
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_add<2u,T,N,1,0,1>
	{
		using Unit = typename CompressedField<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Unit;
		Unit operator ()(const Unit &lhs, const Unit &rhs) const { 
			Unit res{};
			res.w[0] = lhs.w[0] ^ rhs.w[0];
			return res;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_sub
	{
		using Unit = typename CompressedField<P,T,N,R,L,static_cast<uint64_t>(1)>::Unit;
		Unit operator ()(const Unit &lhs, const Unit &rhs) const { 
			Unit res{};
			res.w[0] = lhs.w[0] - rhs.w[0];
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_sub<3u,T,N,static_cast<uint64_t>(1),0u,2u>
	{
		using Unit = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Unit;
		Unit operator()(const Unit &lhs, const Unit &rhs) const {
			const auto temp = lhs.w[0] ^ rhs.w[0];
			const auto x = lhs.w[1] ^ rhs.w[1];
			const auto r0 = temp | x;
			const auto z = temp ^ rhs.w[1];
			const auto w = rhs.w[0] ^ lhs.w[1];
			const auto r1 = z & w;
			Unit res{};
			res.w[0] = r0;
			res.w[1] = r1;
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_sub<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>
	{
		using Unit = typename CompressedField<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Unit;
		Unit operator ()(const Unit &lhs, const Unit &rhs) const { return ::CF_add<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>()(lhs,rhs); }
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_axpy
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		using Element = typename CompressedField<P,T,N,R,L,D>::Element;
		Unit operator ()(const Unit &lhs, const Unit &rhs, const Element a) const { 
			Unit res{};
			res.w[0] = lhs.w[0].axpy(rhs.w[0],a);
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_axpy<3u,T,N,1u,0u,2u>
	{
		using Unit = typename CompressedField<3u,T,N,1u,0u,2u>::Unit;
		using Element = typename CompressedField<3u,T,N,1u,0u,2u>::Element;
		Unit operator()(const Unit &lhs, const Unit &rhs, const Element a) const {
			Unit res{};
			switch (a) {
				case 0:
					res = lhs;
					break;
				case 1:
					res = CF_add<3u,T,N,static_cast<uint64_t>(1),0u,2u>()(lhs,rhs);
					break;
				case 2:
					res = CF_sub<3u,T,N,static_cast<uint64_t>(1),0u,2u>()(lhs,rhs);
			}
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_axpy<2u,T,N,1u,0u,1u>
	{
		using Unit = typename CompressedField<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Unit;
		using Element = typename CompressedField<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Element;
		Unit operator ()(const Unit &lhs, const Unit &rhs, const Element a) const { 
			Unit res{};
			if (a == 0) {
				res = lhs;
			}
			else {
				res = CF_add<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>()(lhs,rhs);
			}
			return res;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_neg
	{
		using Unit = typename CompressedField<P,T,N,R,L,static_cast<uint64_t>(1)>::Unit;
		Unit operator ()(const Unit &lhs) const { 
			Unit res{};
			res.w[0] = res.w[0] - lhs.w[0];
		}
			
	};
	
	template <typename T, uint64_t N>
	struct CF_neg<3u,T,N,static_cast<uint64_t>(1),0u,2u>
	{
		using Unit = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Unit;
		Unit operator()(const Unit &lhs) const {
			Unit res{};
			res.w[0] = lhs.w[0];
			res.w[1] = lhs.w[0] ^ lhs.w[1];
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_neg<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>
	{
		using Unit = typename CompressedField<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Unit;
		Unit operator ()(const Unit &lhs) const { 
			return lhs;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_smul
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		using Element = typename CompressedField<P,T,N,R,L,D>::Element;
		Unit operator ()(const Unit &lhs, const Element a) const { 
			Unit res{};
			res.w[0] = lhs.w[0] * a;
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_smul<3u,T,N,static_cast<uint64_t>(1),0u,2u>
	{
		using Unit = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Unit;
		using Element = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Element;
		Unit operator()(const Unit &lhs, const Element a) const {
			Unit res{};
			switch (a) {
				case 1:
					res = lhs;
					break;
				case 2:
					res = neg(lhs);
			}
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_smul<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>
	{
		using Unit = typename CompressedField<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Unit;
		using Element = typename CompressedField<2u,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Element;
		Unit operator()(const Unit &lhs, const Element a) const {
			Unit res{};
			if (a == 1) {
				res = lhs;
			}
			return res;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_get_entry
	{
		using Unit =    typename CompressedField<P,T,N,R,L,D>::Unit;
		using Element = typename CompressedField<P,T,N,R,L,D>::Element;
		Element operator ()(const Unit &lhs, uint64_t idx) const { 
			return lhs.getEntry(0,idx);
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_get_entry<3u,T,N,static_cast<uint64_t>(1),0u,2u>
	{
		using Unit = typename    CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Unit;
		using Element = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Element;
		Element operator()(const Unit &lhs, uint64_t idx) const {
			return lhs.getEntry(0,idx) + lhs.getEntry(1,idx);
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_set_entry
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		using Element = typename CompressedField<P,T,N,R,L,D>::Element;
		Unit& operator ()(Unit &lhs, const Element e, uint64_t idx) const { 
			lhs.setEntry(0,idx,e);
			return lhs;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_set_entry<3u,T,N,static_cast<uint64_t>(1),0u,2u>
	{
		using Unit = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Unit;
		using Element = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Element;
		Unit& operator()(Unit &lhs, const Element e, uint64_t idx) const {
			lhs.setEntry(0,idx,static_cast<typename Unit::Base>(e != 0));
			lhs.setEntry(1,idx,static_cast<typename Unit::Base>(e == 2));
			return lhs;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_normalize
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		Unit operator ()(const Unit &lhs) const {
			Unit res{};
			//for (uint64_t j = 0; j < Unit::Word::num_bases(); ++j) {
			for (uint64_t i = 0; i < Unit::Word::entries_per_base(); ++i) {
				for (uint64_t j = 0; j < Unit::Word::num_bases(); ++j) {
					res.w[0].b[j] |= ((((lhs.w[0].b[j]) >> (i*(R+L))) & (Unit::Word::entry_mask())) % P) << (i*(R+L));
				}
			}
			return res;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t D>
	struct CF_normalize<P,T,N,static_cast<uint64_t>(1),0u,D>
	{
		using Unit = typename CompressedField<P,T,N,static_cast<uint64_t>(1),0u,D>::Unit;
		Unit operator()(const Unit &lhs) const {
			return lhs;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_pnormalize
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		Unit operator ()(const Unit &lhs) const {
			return CF_normalize<P,T,N,R,L,D>(lhs);
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t D>
	struct CF_pnormalize<P,T,N,static_cast<uint64_t>(1),0u,D>
	{
		using Unit = typename CompressedField<P,T,N,static_cast<uint64_t>(1),0u,D>::Unit;
		Unit operator()(const Unit &lhs) const {
			return lhs;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_get_entries
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		using Base = typename CompressedField<P,T,N,R,L,D>::Base;
		Base operator ()(const Unit * lhs, uint64_t idx, uint64_t t, uint64_t d) const {
			const auto m = (Unit::Word::base_select_mask() << d);
			const auto off = Unit::Word::make_index(idx);
			Base res{};
			if (idx + t > Unit::Word::entries_per_word()) {
				//take from two words (last base of lhs, first base of lhs+1)
				const auto s = (lhs[0].w[0].b[off.first] >> off.second) & Unit::Word::entries_mask();
				const auto r = s | (lhs[1].w[0].b[0] << ((R+L)*Unit::Word::entries_per_base() - off.second));
				const auto bits = _pext_u64(static_cast<uint64_t>(r), static_cast<uint64_t>(m));
				res = bits;
			}
			else {
				//take from one
				if ((off.second + (R+L)*t) > (R+L)*Unit::Word::entries_per_base()) {
					//take from two bases in one word
					const auto s = (lhs[0].w[0].b[off.first] >> off.second) & Unit::Word::entries_mask();
					const auto r = s | (lhs[0].w[0].b[off.first + 1] << ((R+L)*Unit::Word::entries_per_base() - off.second));
					const auto bits = _pext_u64(static_cast<uint64_t>(r), static_cast<uint64_t>(m));
					res = bits;
				}
				else {
					//take from one base
					const auto bits = _pext_u64(static_cast<uint64_t>((lhs[0].w[0].b[off.first] >> off.second) & Unit::Word::entries_mask() ), static_cast<uint64_t>(m)) ;
					res = bits;
				}
			}
			return res & ((static_cast<Base>(1) << t) - 1);
		}
	};
	
	template <uint64_t P, typename T, uint64_t N>
	struct CF_get_entries<P,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>
	{
		using Unit = typename CompressedField<P,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Unit;
		using Base = typename CompressedField<P,T,N,static_cast<uint64_t>(1),0u,static_cast<uint64_t>(1)>::Base;
		Base operator()(const Unit * lhs, uint64_t idx, uint64_t t, uint64_t d) const {
			const auto off = Unit::Word::make_index(idx);
			Base res{};
			if (idx + t > Unit::Word::entries_per_word()) {
				//take from two (last base of lhs, first base of lhs+1)
				auto bits1 = (lhs[0].w[d].b[off.first] >> off.second) & (Unit::Word::entries_mask());
				auto bits2 = (lhs[1].w[d].b[0] << ((1)*Unit::Word::entries_per_base() - off.second));
				res = (bits1 | bits2) & ((static_cast<Base>(1) << t) - 1);
			}
			else {
				//take from one
				if (off.second + t > Unit::Word::entries_per_base()) {
					//take from two bases in one word
					auto bits1 = (lhs[0].w[d].b[off.first] >> off.second) & (Unit::Word::entries_mask());
					auto bits2 = (lhs[0].w[d].b[off.first+1] << ((1)*Unit::Word::entries_per_base() - off.second));
					res = (bits1 | bits2) & ((static_cast<Base>(1) << t) - 1);
				}
				else {
					//take from one base
					auto bits = (lhs[0].w[d].b[off.first] >> off.second) & ((static_cast<Base>(1) << t) - 1);
					res = bits;
				}
			}
			return res;
		}
	};
	
	template <typename T, uint64_t N>
	struct CF_get_entries<3u,T,N,static_cast<uint64_t>(1),0u,2u>
	{
		using Unit = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Unit;
		using Base = typename CompressedField<3u,T,N,static_cast<uint64_t>(1),0u,2u>::Base;
		Base operator()(const Unit * lhs, uint64_t idx, uint64_t t, uint64_t d) const {
			const auto off = Unit::Word::make_index(idx);
			Base res{};
			if (idx + t > Unit::Word::entries_per_word()) {
				//take from two (last base of lhs, first base of lhs+1)
				if (d == 0) {
					auto bits1 = ((lhs[0].w[0].b[off.first] ^ lhs[0].w[1].b[off.first]) >> off.second) & (Unit::Word::entries_mask());
					auto bits2 = ((lhs[1].w[0].b[0] ^ lhs[1].w[1].b[0]) << ((1)*Unit::Word::entries_per_base() - off.second));
					res = (bits1 | bits2) & ((static_cast<Base>(1) << t) - 1);
				}
				else {
					auto bits1 = (lhs[0].w[d].b[off.first] >> off.second) & (Unit::Word::entries_mask());
					auto bits2 = (lhs[1].w[d].b[0] << ((1)*Unit::Word::entries_per_base() - off.second));
					res = (bits1 | bits2) & ((static_cast<Base>(1) << t) - 1);
				}
			}
			else {
				//take from one
				if (off.second + t > Unit::Word::entries_per_base()) {
					//take from two bases in one word
					if (d == 0) {
						auto bits1 = ((lhs[0].w[0].b[off.first] ^ lhs[0].w[1].b[off.first]) >> off.second) & (Unit::Word::entries_mask());
						auto bits2 = ((lhs[0].w[0].b[off.first+1] ^ lhs[0].w[1].b[off.first+1]) << ((1)*Unit::Word::entries_per_base() - off.second));
						res = (bits1 | bits2) & ((static_cast<Base>(1) << t) - 1);
					}
					else {
						auto bits1 = (lhs[0].w[d].b[off.first] >> off.second) & (Unit::Word::entries_mask());
						auto bits2 = (lhs[0].w[d].b[off.first+1] << ((1)*Unit::Word::entries_per_base() - off.second));
						res = (bits1 | bits2) & ((static_cast<Base>(1) << t) - 1);
					}
				}
				else {
					//take from one base
					if (d == 0) {
						auto bits = ((lhs[0].w[0].b[off.first] ^ lhs[0].w[1].b[off.first])>> off.second) & ((static_cast<Base>(1) << t) - 1);
						res = bits;
					}
					else {
						auto bits = (lhs[0].w[d].b[off.first] >> off.second) & ((static_cast<Base>(1) << t) - 1);
						res = bits;
					}
				}
			}
			return res;
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
	struct CF_get_entries_w
	{
		using Unit = typename CompressedField<P,T,N,R,L,D>::Unit;
		using Base = typename CompressedField<P,T,N,R,L,D>::Base;
		Base operator ()(const Unit &lhs, uint64_t widx, uint64_t idx, uint64_t t, uint64_t d) const {
			auto bits = _pext_u64(static_cast<uint64_t>(lhs.w[0].b[widx]), static_cast<uint64_t>(Unit::Word::base_select_mask() << d));
			return static_cast<Base>(bits >> ((R+L)*idx)) & ((static_cast<Base>(1) << t) - 1);
		}
	};
	
	template <uint64_t P, typename T, uint64_t N, uint64_t D>
	struct CF_get_entries_w<P,T,N,static_cast<uint64_t>(1),0u,D>
	{
		using Unit = typename CompressedField<P,T,N,static_cast<uint64_t>(1),0u,D>::Unit;
		using Base = typename CompressedField<P,T,N,static_cast<uint64_t>(1),0u,D>::Base;
		
		Base operator ()(const Unit &lhs, uint64_t widx, uint64_t idx, uint64_t t, uint64_t d) const {
			return (lhs.w[d].b[widx] >> idx) & ((static_cast<Base>(1) << t) - 1);
		}
	};
	
};

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
uint64_t 
CompressedField<P,T,N,R,L,D>::normalize_freq()
{
	return ::CF_normalize_freq<P,T,N,R,L,D>()();
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
uint64_t 
CompressedField<P,T,N,R,L,D>::addin_freq()
{
	return ::CF_addin_freq<P,T,N,R,L,D>()();
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
uint64_t 
CompressedField<P,T,N,R,L,D>::axpy_freq()
{
	return ::CF_axpy_freq<P,T,N,R,L,D>()();
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
uint64_t 
CompressedField<P,T,N,R,L,D>::m4rm_freq()
{
	return ::CF_m4rm_freq<P,T,N,R,L,D>()();
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit  
CompressedField<P,T,N,R,L,D>::add    (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs, const typename CompressedField<P,T,N,R,L,D>::Unit &rhs) const
{
	return ::CF_add<P,T,N,R,L,D>()(lhs,rhs);
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit  
CompressedField<P,T,N,R,L,D>::sub    (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs, const typename CompressedField<P,T,N,R,L,D>::Unit &rhs) const
{
	return ::CF_sub<P,T,N,R,L,D>()(lhs,rhs);
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit  
CompressedField<P,T,N,R,L,D>::axpy   (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs, const typename CompressedField<P,T,N,R,L,D>::Unit &rhs, const typename CompressedField<P,T,N,R,L,D>::Element a) const
{
	return ::CF_axpy<P,T,N,R,L,D>()(lhs,rhs,a);
}	
		
template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit  
CompressedField<P,T,N,R,L,D>::neg    (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs) const
{
	return ::CF_neg<P,T,N,R,L,D>()(lhs);
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit  
CompressedField<P,T,N,R,L,D>::smul   (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs, const typename CompressedField<P,T,N,R,L,D>::Element a) const
{
	return ::CF_smul<P,T,N,R,L,D>()(lhs, a);
}
		
template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Element  
CompressedField<P,T,N,R,L,D>::getEntry  (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs, uint64_t idx) const
{
	return ::CF_get_entry<P,T,N,R,L,D>()(lhs,idx);
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit&   
CompressedField<P,T,N,R,L,D>::setEntry  (typename CompressedField<P,T,N,R,L,D>::Unit &lhs, const typename CompressedField<P,T,N,R,L,D>::Element e, uint64_t idx) const
{
	return ::CF_set_entry<P,T,N,R,L,D>()(lhs,e,idx);
}

		
template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit 
CompressedField<P,T,N,R,L,D>::normalize  (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs) const
{
	return ::CF_normalize<P,T,N,R,L,D>()(lhs);
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit 
CompressedField<P,T,N,R,L,D>::pnormalize (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs) const
{
	return ::CF_pnormalize<P,T,N,R,L,D>()(lhs);
}
		
template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Base 
CompressedField<P,T,N,R,L,D>::getEntries (const typename CompressedField<P,T,N,R,L,D>::Unit * lhs, uint64_t idx, uint64_t t, uint64_t d) const // t entries starting at idx for d-th use in four russians
{
	return ::CF_get_entries<P,T,N,R,L,D>()(lhs,idx,t,d);
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Base 
CompressedField<P,T,N,R,L,D>::getEntries (const typename CompressedField<P,T,N,R,L,D>::Unit &lhs, uint64_t widx, uint64_t idx, uint64_t t, uint64_t d) const
{
	return ::CF_get_entries_w<P,T,N,R,L,D>()(lhs,widx,idx,t,d);
}
		
template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
typename CompressedField<P,T,N,R,L,D>::Unit 
CompressedField<P,T,N,R,L,D>::random () const
{
	typename CompressedField<P,T,N,R,L,D>::Unit res{};
	for (uint64_t i = 0; i < CompressedField<P,T,N,R,L,D>::Word::entries_per_word(); ++i) {
		initEntry(res, rand(), i);
	}
	return res;
}

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
inline
std::ostream& 
CompressedField<P,T,N,R,L,D>::write (std::ostream &os, const typename CompressedField<P,T,N,R,L,D>::Unit &lhs) const
{
	for (uint64_t i = 0; i < CompressedField<P,T,N,R,L,D>::Word::entries_per_word(); ++i) {
		os << getEntry(lhs,i);
	}
	return os;
}

#endif