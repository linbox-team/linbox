#ifndef __compressed_matrix_inl__
#define __compressed_matrix_inl__

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::get_words (uint64_t loff, uint64_t sc, uint64_t nc)
{
	return (CompressedMatrix<Field_>::get_loff(loff, sc) + nc + CompressedMatrix<Field_>::Word::entries_per_word() - 1) / CompressedMatrix<Field_>::Word::entries_per_word();
}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::get_loff  (uint64_t loff, uint64_t sc)
{
	return (loff + sc) % CompressedMatrix<Field_>::Word::entries_per_word();
}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::get_roff  (uint64_t loff, uint64_t sc, uint64_t nc)
{
	return (loff + sc + nc) % CompressedMatrix<Field_>::Word::entries_per_word();
}

template <typename Field_>
CompressedMatrix<Field_>::CompressedMatrix  (typename CompressedMatrix<Field_>::Field &f_)
	: rows   (0)
	, cols   (0)
	, words  (0)
	, stride (0)
	, loff   (0)
	, roff   (0)
	, f      (f_)
	, alloc  (false)
	, rep    (nullptr)
{}

template <typename Field_>
CompressedMatrix<Field_>::CompressedMatrix  (typename CompressedMatrix<Field_>::Field &f_, uint64_t nr, uint64_t nc)
	: rows   (nr)
	, cols   (nc)
	, words  (CompressedMatrix<Field_>::get_words(0,0,nc))
	, stride (CompressedMatrix<Field_>::get_words(0,0,nc))
	, loff   (CompressedMatrix<Field_>::get_loff(0,0))
	, roff   (CompressedMatrix<Field_>::get_roff(0,0,nc))
	, f      (f_)
	, alloc  (true)
	, rep    (new typename CompressedMatrix<Field_>::Unit[nr * CompressedMatrix<Field_>::get_words(0,0,nc)])	
{
	if (rep == nullptr) {
		alloc = false;
		std::cout << "failed to allocate new CompressedMatrix" << std::endl;
	}
}

template <typename Field_>
CompressedMatrix<Field_>::CompressedMatrix  (const CompressedMatrix<Field_> &other)
	: rows   (other.rows)
	, cols   (other.cols)
	, words  (other.words)
	, stride (other.stride)
	, loff   (other.loff)
	, roff   (other.roff)
	, f      (other.f)
	, alloc  (false)
	, rep    (other.rep)
{}

template <typename Field_>
CompressedMatrix<Field_>::~CompressedMatrix ()
{
	if (alloc == true) {
		delete[] rep;
	}
}

template <typename Field_>
inline
CompressedMatrix<Field_>& 
CompressedMatrix<Field_>::operator = (const CompressedMatrix<Field_> &other)
{
	rows   = other.rows;
	cols   = other.cols;
	words  = other.words;
	stride = other.stride;
	loff   = other.loff;
	roff   = other.roff;
	f      = other.f;
	alloc  = false;
	rep    = other.rep;
	
	return *this;
}

template <typename Field_>
CompressedMatrix<Field_>::CompressedMatrix  (const CompressedMatrix<Field_> &other, uint64_t sr, uint64_t sc, uint64_t nr, uint64_t nc)
	: rows   (nr)
	, cols   (nc)
	, words  (CompressedMatrix<Field_>::get_words(other.loff, sc, nc))
	, stride (other.stride)
	, loff   (CompressedMatrix<Field_>::get_loff(other.loff, nc))
	, roff   (CompressedMatrix<Field_>::get_roff(other.loff, nc, sc))
	, f      (other.f)
	, alloc  (false)
	, rep    (other.rep + sr*((other.loff + sc) / CompressedMatrix<Field_>::Word::entries_per_word()))
{}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::rowdim () const
{
	return rows;
}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::coldim () const
{
	return cols;
}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::getwords  () const
{
	return words;
}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::getstride () const
{
	return stride;
}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::getloff   () const
{
	return loff;
}

template <typename Field_>
inline
uint64_t 
CompressedMatrix<Field_>::getroff   () const
{
	return roff;
}

template <typename Field_>
inline
typename CompressedMatrix<Field_>::Field& 
CompressedMatrix<Field_>::field ()
{
	return field;
}

template <typename Field_>
inline
CompressedMatrix<Field_>& 
CompressedMatrix<Field_>::submatrix (const CompressedMatrix<Field_> &other, uint64_t sr, uint64_t sc, uint64_t nr, uint64_t nc)
{
	if (alloc == true) {
		delete[] rep;
	}
	
	CompressedMatrix temp(other, sr, sc, nr, nc);
	*this = temp;
	return *this;
}

template <typename Field_>
inline
void 
CompressedMatrix<Field_>::clear()
{
	if (words == 1) {
		for (uint64_t i = 0; i < rows; ++i) {
			rep[i*stride] = rep[i*stride].deselect(loff, cols);
		}
	}
	else {
		auto sm = typename CompressedMatrix<Field_>::Word (0u, loff);
		auto em = typename CompressedMatrix<Field_>::Word (roff, CompressedMatrix<Field_>::Word::entries_per_word());
		for (uint64_t i = 0; i < rows; ++i) {
			rep[i*stride] &= em;
			for (uint64_t j = 1; j < words-1; ++j) {
				rep[i*stride + j].clear();
			}
			rep[i*stride+words-1] &= em;
		}
	}
}

template <typename Field_>
inline
typename CompressedMatrix<Field_>::Element 
CompressedMatrix<Field_>::getEntry (uint64_t ri, uint64_t ci) const
{
	auto wo = (ci + loff) / CompressedMatrix<Field_>::Word::entries_per_word();
	auto bo = (ci + loff) % CompressedMatrix<Field_>::Word::entries_per_word();
	
	return f.getEntry(rep[ri*stride + wo], bo);
}

template <typename Field_>
inline
void 
CompressedMatrix<Field_>::setEntry (unsigned ri, uint64_t ci, const typename CompressedMatrix<Field_>::Element eij)
{
	auto wo = (ci + loff) / CompressedMatrix<Field_>::Word::entries_per_word();
	auto bo = (ci + loff) % CompressedMatrix<Field_>::Word::entries_per_word();
	
	f.setEntry(rep[ri*stride + wo], eij, bo);
}

template <typename Field_>
inline
const typename CompressedMatrix<Field_>::Unit * 
CompressedMatrix<Field_>::getRow (uint64_t i) const
{
	return rep + i*stride;
}

template <typename Field_>
inline
void 
CompressedMatrix<Field_>::random()
{
	for (uint64_t i = 0; i < rows; ++i) {
		for (uint64_t j = 0; j < cols; ++j) {
			setEntry(i, j, f.rand_entry());
		}
	}
}

template <typename Field_>
inline
std::ostream& 
CompressedMatrix<Field_>::write (std::ostream &os) const
{
	for (uint64_t i = 0; i < rows; ++i) {
		for (uint64_t j = 0; j < cols; ++j) {
			os << getEntry(i, j);
		}
		os << std::endl;
	}
	return os;
}

#endif