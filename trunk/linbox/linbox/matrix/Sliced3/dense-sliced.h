#ifndef __DENSE_SLICED_H
#define __DENSE_SLICED_H

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>

/**
  The matrix class Sliced is defined.  
  It adheres to the LinBox dense matrix interface.

  It depends on SlicedBase, also defined here, which packs GF(3) elements 
  into a pair of ints. The int type is a template parameter.
*/

#include "dense-matrix.h"
//#include <linbox/util/timer.h>
//#include "sliced-stepper.h"

using namespace std;
using namespace LinBox;

/* SLICED BASE CODE 
   SlicedBase implements functios on a vector of GF(3) elements.
   T is an integer type of some length n. 
   The vectors are of length n, packed into two T words in sliced fashion.
   (Each element partakes of one bit from b0 and one bit from b1).

   Most of the SlicedBase functios are in C++ operator form. 
 */
template<class T>
struct SlicedBase
{
public:
	T& bits0(){return b0;}
	T& bits1(){return b1;}

	SlicedBase & operator=(const T &rhs){
		b0 = rhs;
		b1 = rhs;
		return *this;
	}

	SlicedBase & operator=(const SlicedBase &rhs){
		b0 = rhs.b0;
		b1 = rhs.b1;
		return *this;
	}

	//  only nontrivial multiplication in F(3) is *=2
	//  which is negation, which is word swap
	//  between the bits
	SlicedBase & operator*=(const T &two){
		//  this is multiplication by two only!!!!
		//  we never even read the arg.
		b1 ^= b0;
		return *this;
	}

	SlicedBase operator*(const T &two) const{
		return SlicedBase(*this) *= two;	
	}

	//  used for masking both words
	SlicedBase & operator&=(const T &mask){
		b1 &= mask;
		b0 &= mask;
		return *this;
	}

	SlicedBase operator & (const T &rhs){ 
		return SlicedBase(*this) &= rhs;	
	}

	//  used for shifting both words
	SlicedBase & operator>>=(const T &shift){
		b1 >>= shift;
		b0 >>= shift;
		return *this;
	}

	SlicedBase operator >> (const T &rhs){ 
		return SlicedBase(*this) >>= rhs;	
	}

	//  used for shifting both words
	SlicedBase & operator<<=(const T &shift){
		b1 <<= shift;
		b0 <<= shift;
		return *this;
	}

	SlicedBase operator << (const T &rhs){ 
		return SlicedBase(*this) <<= rhs;	
	}

	//  used for combining value into both words
	SlicedBase & operator|=(const T &val){
		b1 |= val;
		b0 |= val;
		return *this;
	}

	//  used for combining two units
	SlicedBase & operator|=(const SlicedBase &rhs){
		b1 |= rhs.b1;
		b0 |= rhs.b0;
		return *this;
	}

	SlicedBase & operator+=(const SlicedBase &rhs){
		T a = b0 ^ rhs.b1;
		T b = b1 ^ rhs.b0;
		T s = a ^ b1;
		T t = b ^ rhs.b1;
		b1 = a & b;
		b0 = s | t;

		return *this;
	}

	const SlicedBase operator+(const SlicedBase &rhs){
		return SlicedBase(*this) += rhs;	
	}
	
	//  comparison ops
	bool operator==(const SlicedBase &rhs){
		return b0 == rhs.bo && b1 == rhs.b1;
	}

	//  used for combining two units
	bool operator!=(const SlicedBase &rhs){
		return b0 != rhs.b0 || b1 != rhs.b1;
	}

	SlicedBase & zero(){
		b0 = 0;
		b1 = 0;
		return *this;
	}

//protected:
	T b0;
	T b1;
};

/**
	The Sliced Matrix class 
	_Domain must be  a GF(3) rep, BaseT must be an unsigned int type.
	TODO more docs, discuss row/col slicing
*/
template <class _Domain>
class Sliced : public DenseMatrix<SlicedBase<typename _Domain::Word_T> >
{
public:
	typedef _Domain Domain;
	typedef typename Domain::Scalar Scalar;
	typedef typename Domain::Word_T SlicedWord;
	typedef SlicedBase<SlicedWord> SlicedUnit;

	typedef DenseMatrix<SlicedBase<typename _Domain::Word_T> > Base_T;
	using Base_T::rawBegin;
	using Base_T::rawEnd;
	using Base_T::rowBegin;
	using Base_T::rowEnd;
	using Base_T::_stride;
	using Base_T::_rep;
	using Base_T::_alloc;
	typedef typename Base_T::RawIterator RawIterator;
	typedef Base_T Matrix;

	enum op_t { ADD, SMUL, AXPY, ZERO, COPY };   

	Sliced() : Matrix(), _domain(Domain()), _m(0), _n(0), _colPacked(false),
		_SIZE(8*sizeof(SlicedWord)), _sub(false) {}

	Sliced(Domain d) : Matrix(), _domain(d), _m(0), _n(0), _colPacked(false),
		_SIZE(8*sizeof(SlicedWord)), _sub(false) {}

	//  Constructor taking dims and bool defaulting ROWPACKED
	Sliced (Domain d, size_t m, size_t n, bool colP = false) :
		Matrix(), _domain(d), _m(m), _n(n), _colPacked(colP), _i(0), 
		_j(0), _loff(0), _roff(0), _sub(0) {
			Matrix::init(colP ? (m + _SIZE - 1)/_SIZE : m, 
					colP ? n : (n + _SIZE - 1)/_SIZE);
	}

	//  need to override parent size, because we are a different size
	void size(size_t m = 0, size_t n = 0, bool cp = false){
		_m = m; _n = n; _colPacked = cp;
		_i = _j = _roff = _loff = _sub = 0;

		if(_colPacked)
			Matrix::size((m + _SIZE - 1)/_SIZE, n);
		else
			Matrix::size(m, (n + _SIZE - 1)/_SIZE);
	}

	//  need to override parent inits, because we are a different size
	void init(size_t m = 0, size_t n = 0) {
		_m = m; _n = n;
		_i = _j = _roff = _loff = _sub = 0;
 		if(_colPacked)
			Matrix::size((m + _SIZE - 1)/_SIZE, n);
		else
			Matrix::init(m, (n + _SIZE - 1)/_SIZE);
	}

	void init(size_t m, size_t n, const Scalar &filler){
		init(m,n);

		//  TODO stepper functionality, more efficient.
		for(size_t i=0; i<m; ++i)
			for(size_t j=0; j<n; ++j)
				setEntry(i, j, filler);
	}

	Sliced & submatrix(Sliced &other, size_t i, size_t j, size_t m, size_t n){
		_m = m;
		_n = n;
		_colPacked = other._colPacked;
		_i = i;
		_sub = true;

		size_t headStart = other._loff ? _SIZE - other._loff : 0;
		_j = j + headStart;
		//  determine in which of parents' columns do we begin
		size_t firstCol =  _j/_SIZE;

		_loff = (_SIZE-(_j%_SIZE))%_SIZE;  // could prob. omit final mod?
		_roff = (_j+n)%_SIZE;
		
		// _roff really an offset? or just goes to "parent" matrix's normal end
		//if( (_j+n) == other.coldim() )
		//	_roff = 0;
			//  AND IF our super matrix is < SIZE, special case we don't account for
			/*  ^ this can cause the writing of a superfluous zero word in s_write_bin
				if last word is next to first word in a submat */

		//std::cerr << _j << " " << n << other.coldim();
		//cerr << "LOFF: " << _loff << " ROFF: " << _roff << endl;
		//std::cerr << "ijmn?" << i << " " << j << " " << m << " " << n << std::endl;
		//cerr << "firstCol: " << firstCol << " _j: " << _j << " cols(): " << _cols << " other.cols(): " << other.cols() << endl;
		
		//TODO colpacked version
		this->Matrix::submatrix(other, i, firstCol, m, (n + (_j%_SIZE) + (_SIZE - 1))/_SIZE);

		return *this;
	}

	//  blindly assumes "other" matrix is rowpacked...
	Sliced (Sliced &other, size_t i, size_t j, size_t m, size_t n) : Matrix(), _domain(other._domain) { 
		//  we currently can't handle submatrices that reside entirely within
		//  a single sliced unit, width-wise...
		//  meaning < _SIZE columns, while both borders are in the same word
		size_t begin = j+other._loff;
		size_t end = begin+n;
		size_t endTest = other._n < _SIZE ? other._n : _SIZE;
		//  check if we're small enough, then check if we're on a border (which is OK)
		if(n < _SIZE && begin%_SIZE && end%endTest){
			size_t test = other._loff ? other._loff : _SIZE;
			if((j+n) <= test){
				cerr << "Unsupported submatrix request.  Submatrix entirely within a word" << 
					" while not aligned with a left or right border" <<	endl;
				//exit(-1);
			}
		}

		//  passed the security checkpoint.  on to business
		submatrix(other, i, j, m, n);
	}

//  incl functions to handle specialized cases (where submats are not word-aligned)
#include "dense-sliced.inl"

	//  typical addin, barge right through with the raw iterator
	Sliced& addin(Sliced &other){
		if(_loff || _roff)
			return s_addin(other);

		RawIterator a = rawBegin();
		RawIterator c = rawEnd();
		RawIterator b = other.rawBegin();

		//for(; &(*a) != &(*c); ++a, ++b)
		for(; a != c; ++a, ++b)
			(*a) += (*b);
		
		return *this;
	}

	Sliced& smulin(Scalar &x){
		if(x == 1) 
			return *this;

		if(_loff || _roff)
			return s_smulin(x);

		if(x == 0)
			return zero();

		//  x == 2
		//int j = 0;
		//for(RawIterator a = rawBegin(); j < rows()*cols() ; j+=_SIZE, ++a)
		for(RawIterator a = rawBegin(); a != rawEnd(); ++a)
			(*a)*=2;

		return *this;
	}

	/*
	Sliced& smul(Scalar &x){
		return Sliced(*this).copy(*this).smulin(x);
	}
	*/

	void setEntry(size_t i, size_t j, const Scalar &a_ij){
		SlicedWord e = static_cast<SlicedWord>(a_ij);
		//  determine location
		//  TODO:  THIS SEEMS TO EXCPET _rep isn't adjusted for submatrix
		//  whereas in dense-matrix.h, submatrix adjusts _rep.  need to study 
		//  the pros/cons of either approach but CANNOT MIX THEM
		//  I think 
		//size_t word = (_i+i)*_stride + ((j+(_j%_SIZE))/_SIZE);
		size_t word = i*_stride + ((j+(_j%_SIZE))/_SIZE);
		//int w = i*Matrix::coldim() + ((j+(_j%_SIZE))/_SIZE);
		//int w = ((j+(_j%_SIZE))/_SIZE);
		//RawIterator word = rowBegin(i) + w;
		//word += w;

		size_t index = (_j+j) % _SIZE;

		SlicedWord b1 = (e & 2) >> 1;
		SlicedWord b0 = (e & 1) | b1;
		SlicedWord one = 1;

		//(*word) &= ~(one << index);
		//(*word).b0 |= b0 << index;
		//(*word).b1 |= b1 << index;

		_rep[word] &= ~(one << index);

		_rep[word].b0 |= b0 << index;
		_rep[word].b1 |= b1 << index;
	}

	Scalar& getEntry(Scalar &x, size_t i, size_t j){
		//int w = ((j+(_j%_SIZE))/_SIZE);
		//  TODO:  see note in setEntry
		size_t word = (_i+i)*_stride + ((j+(_j%_SIZE))/_SIZE);
		//RawIterator word = rowBegin(i) + w;
		//word += w;

		size_t index = (_j+j) % _SIZE;

		//cerr << "(w" << w << ",i" << index << ")";
		//std::cerr << (*word).b0 << "x" << (*word).b1 << std::endl;

		//int answer = (int)((((*word).b1 >> index) & 1) + (((*word).b0 >> index) & 1));
		size_t answer = (size_t)(((_rep[word].b1 >> index) & 1) + 
				((_rep[word].b0 >> index) & 1));
		//_domain.init(x, answer);
		x = answer;
		return x;
	}

	//  begin, end, scalar, other begin
	Sliced & axpyin(RawIterator &b, RawIterator &e, Scalar &s, RawIterator &ob){
		RawIterator x = b;
		RawIterator y = ob; 

		switch(static_cast<int>(s)){
			case 0:
				return *this;
			case 1: 
				if(_loff || _roff)
					return s_axpyin(b, s, ob);
				for(; x!=e; ++x,++y)
					(*x) += (*y);
				return *this;
			case 2:
				if(_loff || _roff)
					return s_axpyin(b, s, ob);
				for(; x!=e; ++x,++y){ 
					(*x) += (*y)*2;
				}
				return *this;
		}
		return *this;
	}

	//  MUL:
	//  become the product of two sliced matrices
	//  (only seems to work if row packed so far)
	//  (does not check for compatible sizes)
	//  (does NOT work yet for two submatrices)
	template <class Gettable>
	Sliced & mul(Gettable& A, Sliced& B){
		zero();

		Scalar a_ij;
		//  c&b - begin and end
		RawIterator c_b, c_e, b_b; //, b_e;

		size_t count = 0;
		//  axpyin to C individual entries of A with rows of B
		for(; count < rowdim(); ++count){
			c_b = rowBegin(count);
			c_e = rowEnd(count);

			//Timer t;
			//t.clear(); t.start();
			for(size_t len = 0; len < A.coldim(); ++len){
				//  element of A goes down rows of A.  
				//  (could improve on speed of method to get this value [step?])
				a_ij = A.getEntry(a_ij, count, len);

				//  this is axpy'd to THAT row of B.
				b_b = B.rowBegin(len);
				axpyin(c_b, c_e, a_ij, b_b);
			}
			//t.stop();  std::cerr << t << std::endl;
		}
		return *this;
	}

	std::ostream& write(std::ostream &os = std::cerr, size_t offset = 0){
		Scalar t;
		for(size_t i = 0; i<_m; i++){
			for(size_t q=0; q<offset; ++q)
				os << "  ";
			for(size_t j = 0; j<_n; j++){
				os << (size_t)getEntry(t, i, j);
				/*  DEBUG _SIZE-BIT BOUNDARIES 
				if(offset && j && ! ((_j+j+1)%_SIZE)){
				if(j && ! ((_j+j+1)%_SIZE))
					os << "|";
				else
					os << " ";
				}
				else
				*/
					os << " ";
			}
			os << endl;
		}
		os << endl;
		return os;
	}

	//  returns an (approximately) random _SIZE-bit word
	SlicedWord randomLL(){
		SlicedWord r = 0;
		for(size_t i=0; i<(_SIZE/8); i++)
			r |= ((SlicedWord)(rand()%256) << (i << 3));
		return r;
	}

	//  completely randomizes sliced block entries
	Sliced& random(size_t seed = 0){
		if(seed) srand(seed);
		//srand(seed);

		RawIterator a = rawBegin();
		//  TODO:  problem here is that 
		//  this will cause about 1/2 0's, 1/4 1's & 2's
		for(; a != rawEnd(); ++a){
			(*a).b0 = randomLL();
			(*a).b1 = randomLL() & (*a).b0;

			//std::cerr << (*a).b0 << "x" << (*a).b1 << std::endl;
		}
		return *this;
	}

	//  completely zeros out a sliced block
	Sliced& zero(){
		if(_sub && (_loff || _roff)){
			Scalar t;
			//  TODO (can do [[[slightly]]] better if we specialize for 0
			return s_smulin(_domain.init(t,0));
		}

		//  there has to be a way to use a 
		//  low level mem* function to zero out a memory block
		RawIterator a = rawBegin();
		RawIterator b = rawEnd();
		for(;a != b; ++a)
			(*a).zero();

		return *this;
	}

	bool isEqual(Sliced &other){
        if (rows() != other.rows() || cols() != other.cols()) {
			//std::cout << "shape mismatch" << std::endl; 
			return false;
		}
	//	if(_loff || _roff || other._loff || other._roff){
			//  very slow... could be improved if both mats are aligned
			Scalar x, y;
			for(size_t i = 0; i<rows(); ++i)
				for(size_t j = 0; j<cols(); ++j){
					//std::cerr << getEntry(x, i, j) << "vs" << other.getEntry(y, i, j) << std::endl;
					if(getEntry(x, i, j) != other.getEntry(y, i, j)) {
						//std::cout << "entry mismatch " << i << " " << j << std::endl; 
						return false;
					}
				}
	/*  fails for very small (1 by 1) matrix.
		}
		else{
			//  otherwise just check _SIZE-bits at a time
			//  TODO could be faster w/ memcmp, but not applicable for breaks "aligned" submats
			for(RawIterator a = rawBegin(), b=other.rawBegin(); a != rawEnd(); ++a, ++b){
				//int count = 0;
				if((*a) != (*b)){
					//cerr << (*a).b0 << "ver" << (*b).b0 << endl;
					//cerr << (*a).b1 << "ver" << (*b).b1 << endl;
					return false;
				}
			}
		}
	*/
		return true;
	}

	bool isZero(){
		if(_loff || _roff){
			Scalar x;
			for(size_t i = 0; i<rows(); ++i)
				for(size_t j = 0; j<cols(); ++j)
					if(getEntry(x, i, j))
						return false;
		}
		else{
			//  otherwise just check _SIZE-bits at a time
			//  could prob be faster
			for(RawIterator a = rawBegin(); a != rawEnd(); ++a){
				if((*a).b0 || (*a).b1)
					return false;
			}
		}
		return true;
	}

	Sliced& deepcopy(Sliced &other){
		if(_sub) return s_copy(other);

		*this = other;
		return *this;
	}

	std::ostream& writeRep(std::ostream &os = std::cerr){
		for(RawIterator i=rawBegin(); i!=rawEnd(); ++i){
			os << "0th bits: " << (SlicedWord)(*i).b0 << " 1st bits: " << (SlicedWord)(*i).b1 << endl;	
		}
		os << endl;
		return os;
	}

	size_t memSize(){
		size_t bytesPerRow = ((_n + _SIZE-1)/_SIZE) * 2 * sizeof(SlicedWord);
		return bytesPerRow  * _m;
	}

	std::ostream& writeBinary(std::ostream &os){
		if(_sub){ //std::cerr << "\n\nNear death!\n\n" << std::endl; 
			return s_wb(os); } // TODO fix

		size_t bytes = memSize();
		//std::cerr << "WRITE " << bytes << " BYTES." << std::endl;
		os.write((char *)&(*_rep), bytes);
		
		return os;
	}

	//  assumes we're NOT a submatrix (read into contiguous block)
	std::istream& readBinary(std::istream& is){
		if(_sub){ //std::cerr << "\n\nDeath!\n\n" << std::endl;
			return s_rb(is); }
		size_t bytes = memSize();
		is.read((char *)&(*_rep), bytes);
		// std::cerr << "READ " << bytes << " BYTES." << std::endl;
		return is;
	}

	bool writeBinaryFile(const char *file){
		ofstream bin;
		bin.open(file, ios::out | ios::binary);
		if(!bin){ std::cerr << "failure opening " << file << std::endl; return false; }
		writeBinary(bin);
		bin.close();
		return true;
	}

	bool readBinaryFile(const char *file){
		ifstream in;
		in.open(file, ios::in | ios::binary);
		if(!in){ std::cerr << "failure opening " << file << std::endl; return false; }
		readBinary(in);
		in.close();
		return true;
	}

    void mmapFile(int fd){
		if(_alloc){
			std::cerr << "alloc'd matrix trying to mmap." << std::endl;
			close(fd);
			return;
		}
        size_t bytes = memSize();

        //_rep = (Sliced::SlicedUnit *)mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		//TODO MAP_HUGETLB
        _rep = (Sliced::SlicedUnit *)mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_SHARED , fd, 0);
        if((void *)_rep == MAP_FAILED) {
            close(fd);
            perror("Error mapping the file.");
        }
    }

    void mmapBinaryFile(const char *file){
        int fd = open(file, O_RDWR);
        if(fd == -1) {
            perror("Error opening file to mmap.");
        }
		mmapFile(fd);
    }

    void munmapBinaryFile(){
        size_t bytes = memSize();
        munmap((void *)_rep, bytes);
    }

	_Domain& domain(){ return _domain; }

	size_t rows(){ return _m; }
	size_t cols(){ return _n; }
	size_t rowdim(){ return _m; }
	size_t coldim(){ return _n; }

	size_t l(){ return _loff; }
	size_t r(){ return _roff; }

	//  pointer/offset info for debugging
	void pinfo(){
		std::cerr << "Matrix @ "; 
		RawIterator i = rawBegin();	
		i.pinfo();
		std::cerr << "\t" << _m << " x " << _n << std::endl;
		std::cerr << "\tLO: " << _loff << " RO: " << _roff << std::endl;
		std::cerr << std::endl;
	}
		
	//  size info for debugging
	void sinfo(){
		std::cerr << "Matrix " << _m << " by " << _n << std::endl;
		std::cerr << "... " << _m * _n / 4 << " bytes." << std::endl;
	}

private:
	_Domain _domain;
	size_t _m, _n;
	bool _colPacked;
	//  SUBMATRIX
	size_t _i, _j;
	size_t _loff, _roff; //left & right offsets
	size_t _SIZE;
	bool _sub;
}; // Sliced

#endif // __DENSE_SLICED_H
