#ifndef __compressed_matrix_domain_inl__
#define __compressed_matrix_domain_inl__

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::mul_classical  (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A, const typename CompressedMatrixDomain<Field_>::Matrix &B) const
{
	for (uint64_t i = 0; i < A.rowdim(); ++i) {
		for (uint64_t j = 0; j < A.coldim(); ++j) {
			typename CompressedMatrixDomain<Field_>::Element cij = 0;
			for (uint64_t k = 0; k < B.coldim(); ++k) {
				cij += A.getEntry(i,j) * B.getEntry(j,k);
				cij %= f.cardinality();
			}
			C.setEntry(i, j, cij);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::mul_compressed (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A, const typename CompressedMatrixDomain<Field_>::Matrix &B) const
{
	for (uint64_t i = 0; i < A.rowdim(); ++i) {
		uint64_t rcount = 1;
		for (uint64_t j = 0; j < A.coldim(); ++j) {
			auto e = A.getEntry(i,j);
			//if (e != f.zero) {
				for (uint64_t k = 0; k < C.getwords(); ++k) {
					C.rep[i*C.getstride() + k] = f.axpy(C.rep[i*C.getstride() + k], B.rep[j*B.getstride() + k], e);
				}
				++rcount;
			//}
			if ((Field_::axpy_freq() != 0) and ((rcount >= Field_::axpy_freq()) or (j == A.coldim() - 1))) {
				for (uint64_t k = 0; k < C.getwords(); ++k) {
					C.rep[i*C.getstride() + k] = f.normalize(C.rep[i*C.getstride() + k]);
				}
				rcount = 1;
			}
		}
	}
}

/*
template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::mul_compressed (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A, const typename CompressedMatrixDomain<Field_>::Matrix &B) const
{
	for (uint64_t i = 0; i < A.rowdim(); ++i) {
		uint64_t rcount = 0;
		for (uint64_t j = 0; j < A.coldim(); ++j) {
			auto e = A.getEntry(i,j);
			if (e != f.zero) {
				if ((Field_::axpy_freq() != 0) and (rcount + e*f.cardinality() > ((1u << Field_::Word::entry_length()) - 1))) {
					for (uint64_t k = 0; k < C.getwords(); ++k) {
						C.rep[i*C.getstride() + k] = f.normalize(C.rep[i*C.getstride() + k]);
					}
					rcount = f.cardinality() - 1;
				}
				for (uint64_t k = 0; k < C.getwords(); ++k) {
					C.rep[i*C.getstride() + k] = f.axpy(C.rep[i*C.getstride() + k], B.rep[j*B.getstride() + k], e);
				}
				rcount += e*(f.cardinality() - 1);
			}
			if ((Field_::axpy_freq() != 0) and (j == A.coldim() - 1)) {
				for (uint64_t k = 0; k < C.getwords(); ++k) {
					C.rep[i*C.getstride() + k] = f.normalize(C.rep[i*C.getstride() + k]);
				}
			}
		}
	}
}
*/

template <typename Field_>
template <uint64_t NT, uint64_t K>
inline 
void 
CompressedMatrixDomain<Field_>::mul_four_russians (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A, const typename CompressedMatrixDomain<Field_>::Matrix &B) const
{
	using Unit = typename CompressedMatrixDomain<Field_>::Unit;
	Unit *table;
	table = new Unit[C.getwords() * (static_cast<uint64_t>(1) << K)];
		
	uint64_t rowcount = 0;
	
	for (uint64_t i = 0; i < A.getwords(); ++i) {
		for (uint64_t j = 0; j < Unit::Word::num_bases(); ++j) {
			for (uint64_t k = 0; k < Unit::Word::entries_per_base(); k += K) {
				build_four_russians_table<2>(table, B.getRow(rowcount), C.getwords(), B.getstride(),K);
				rowcount += K;
				for (uint64_t l = 0; l < A.rowdim(); ++l) {
					uint64_t rcount = 2;
					for (uint64_t d = 0; d < Field_::four_russians_depth(); ++d) {
						const auto bits = f.getEntries(A.rep[l*A.getstride() + i], j, k + K, K, d);
						for (uint64_t n = 0; n < C.getwords(); ++n) {
							C.rep[l*C.getstride() + n] = f.axpy(C.rep[l*C.getstride() + n], table[(bits)*C.getwords() + n], static_cast<uint64_t>(1) << d);
						}
						++rcount;
						if (Field_::normalize_freq() != 0 and rcount >= Field_::normalize_freq()) {
							for (uint64_t n = 0; n < C.getwords(); ++n) {
								C.rep[l*C.getstride() + n] = f.normalize(C.rep[l*C.getstride() + n]);
							}
							rcount = 2;
						}
					}
				}
			}
		}
	}
	delete[] table;
}

template <typename Field_>
inline
typename Field_::Base
get_bits(const Field_ &f, typename Field_::Unit *ptr, uint64_t K, uint64_t boff, uint64_t d)
{
	typename Field_::Base bits{};
	if (K > Field_::Word::entries_per_base()) {
		uint64_t k = 0;
		uint64_t off2 = boff;
		for (; k < K - Field_::Word::entries_per_base(); k += Field_::Word::entries_per_base()) {
			bits |= (f.getEntries(ptr, off2, Field_::Word::entries_per_base(), d) << k);
			off2 += Field_::Word::entries_per_base();
			if (off2 >= Field_::Word::entries_per_word()) {
				++ptr;
				off2 -= Field_::Word::entries_per_word();
			}
		}
		bits |= (f.getEntries(ptr, off2, K-k, d) << k);
	}
	else {
		bits = f.getEntries(ptr, boff, K, d);
	}
	return bits;
}

template <typename Field_>
inline 
void 
CompressedMatrixDomain<Field_>::mul_four_russians (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A, const typename CompressedMatrixDomain<Field_>::Matrix &B, uint64_t K) const
{
	using Unit = typename CompressedMatrixDomain<Field_>::Unit;
	Unit *table = new Unit[C.getwords() * (static_cast<uint64_t>(1) << K)];
	if (table == nullptr) {
		std::cout << "uh oh " << K << std::endl;
		return;
	}
	
	uint64_t s = 0;
	uint64_t uoff = 0;
	uint64_t boff = A.loff;
	
	//std::cout << "beginning four russians" << std::endl;
	
	uint64_t rcount = 1;
	while (s < A.coldim()) {
		//std::cout << s << ", " << K  << std::endl;
		build_four_russians_table<2>(table, B.getRow(s), C.getwords(), B.getstride(), K);	
		++rcount;
		for (uint64_t j = 0; j < A.rowdim(); ++j) {
			for (uint64_t d = 0; d < Field_::four_russians_depth(); ++d) {
				const auto bits = get_bits(f,A.rep+j*A.getstride() + uoff, K, boff, d);// & ((1u << K) - 1);
				for (uint64_t n = 0; n < C.getwords(); ++n) {
					C.rep[j*C.getstride() + n] = f.axpy(C.rep[j*C.getstride() + n], table[bits*C.getwords() + n], static_cast<uint64_t>(1) << d);
				}	
			}
			if ((Field_::axpy_freq() != 0) and (rcount >= Field_::m4rm_freq() or (s+K == A.coldim()))) {
				for (uint64_t n = 0; n < C.getwords(); ++n) {
					C.rep[j*C.getstride() + n] = f.normalize(C.rep[j*C.getstride() + n]);
				}
			}
		}
		if (rcount >= Field_::m4rm_freq()) {
			rcount = 1;
		}
		boff += K;
		if (boff >= Field_::Word::entries_per_word()) {
			boff -= Field_::Word::entries_per_word();
			++uoff;
		}
		s += K;
		if (s + K > A.coldim()) {
			K = A.coldim() - s;
		}
	}
		
	delete[] table;
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::addin  (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			C.rep[i*C.getstride()+j] = f.add(C.rep[i*C.getstride()+j], A.rep[i*A.getstride()+j]);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::subin  (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			C.rep[i*C.getstride()+j] = f.sub(C.rep[i*C.getstride()+j], A.rep[i*A.getstride()+j]);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::copy   (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			C.rep[i*C.getstride()+j] = A.rep[i*A.getstride()+j];
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::neg    (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			C.rep[i*C.getstride()+j] = f.neg(C.rep[i*C.getstride()+j], A.rep[i*A.getstride()+j]);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::axpyin (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Matrix &A, const typename CompressedMatrixDomain<Field_>::Element a) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			f.axpyin(C.rep[i*C.getstride()+j], A.rep[i*A.getstride()+j], a);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::negin  (typename CompressedMatrixDomain<Field_>::Matrix &C) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			f.negin(C.rep[i*C.getstride()+j]);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::smulin (typename CompressedMatrixDomain<Field_>::Matrix &C, const typename CompressedMatrixDomain<Field_>::Element c) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			f.smulin(C.rep[i*C.getstride()+j], c);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::zero   (typename CompressedMatrixDomain<Field_>::Matrix &C) const
{
	C.clear();
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::normalize   (typename CompressedMatrixDomain<Field_>::Matrix &C) const
{
	for (uint64_t i = 0; i < C.rowdim(); ++i) {
		for (uint64_t j = 0; j < C.getwords(); ++j) {
			C.rep[i*C.getstride()+j] = f.normalize(C.rep[i*C.getstride()+j]);
		}
	}
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::random (typename CompressedMatrixDomain<Field_>::Matrix &C) const
{
	C.random();
}

template <typename Field_>
inline
void 
CompressedMatrixDomain<Field_>::print  (typename CompressedMatrixDomain<Field_>::Matrix &C, std::ostream &os) const
{
	C.write(os);
}

template <typename Field_>
template <uint64_t S> 
inline
void 
CompressedMatrixDomain<Field_>::build_four_russians_table (typename CompressedMatrixDomain<Field_>::Unit *dest, const typename CompressedMatrixDomain<Field_>::Unit *src, uint64_t dw, uint64_t ss, uint64_t K) const
{
	uint64_t baserow = 1;
	uint64_t oldbaserow = 0;
	int rcount = 1;
	for (uint64_t i = 0; i < K; ++i) {
		// copy new row
		for (uint64_t j = 0; j < dw; ++j) {
			dest[baserow*dw + j] = src[i*ss + j];
		}
		for (uint64_t k = 1; k < baserow; ++k) {
			for (uint64_t l = 0; l < dw; ++l) {
				dest[(baserow+k)*dw+l] = f.add(dest[k*dw+l], dest[baserow*dw + l]);
			}
		}
		++rcount;
		
		if ((Field_::addin_freq() != 0) and ((rcount >= Field_::addin_freq()) or (i = K-1))) {
			for (uint64_t j = oldbaserow; j < 2*baserow - oldbaserow; ++j) {
				for (uint64_t k = 0; k < dw; ++k) {
					dest[j*dw + k] = f.normalize(dest[j*dw + k]);
					//std::cout << j << std::endl;
				}
			}
			oldbaserow = 2*baserow;
			rcount = 1;
		}
		
		baserow *= 2;
	}
}

#endif