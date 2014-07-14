#ifndef __LINBOX_sliced3_SMI_H
#define __LINBOX_sliced3_SMI_H

namespace LinBox{

template<typename E>
struct SubMatIterator { // periodic gaps iterator
	E * _me, * _end;
	size_t _s1; //stride from end of a row to beginning of next.
	size_t _s2; //stride from one row the next (same position in each).

	SubMatIterator () : _me(0), _end(0), _s1(0), _s2(0) {
	}

	SubMatIterator (E* loc, size_t c = 0, size_t stride = 0) 
		: _me(loc), _end(loc+c), _s1(stride - c), _s2(stride) {
//			std::cerr << "created with " << c << std::endl;
	}

	SubMatIterator(const SubMatIterator& p) 
		: _me(p._me), _end(p._end), _s1(p._s1), _s2(p._s2) {
	}

	SubMatIterator& operator=(const SubMatIterator& p) {
		_me = p._me; _end = p._end; _s1 = p._s1; _s2 = p._s2;
		return *this;
	}

	/*  debugging */
	void pinfo(){
		std::cerr << "Iterator info; me: " << _me << " end: " << _end << std::endl;
	}

	E& operator*() const { return *_me; 
	}

	// untested
	SubMatIterator operator+(int n) {
		SubMatIterator p(*this); return p += n;
	}

	SubMatIterator& operator+=(size_t n) {
		size_t q = n/(_s2 - _s1); // q is number of whole rows to jump
		std::ptrdiff_t r = n - q*(_s2 - _s1);
		size_t m = q*_s2;
		_me += m + r; _end += m;
		if (_end <= _me) { _me += _s1; _end += _s2; }
		return *this;
	} 

	SubMatIterator& operator++() { 
		if (++_me == _end) { _me += _s1;  _end += _s2; } 
		return *this; 
		/*
		bool a = (++_me == _end);
		_me += (_s1 && a);
		_end += (_s2 && a)
		*/
	}

	// untested
	SubMatIterator operator++(int) {
		SubMatIterator p(*this); ++*this; return p;
	} 

	// untested
	SubMatIterator& operator--() {
		//if (--_me == _end + (_s2-_s1) - 1) { _me -= _s1; _end -= _s2; }		
		//Above broken. I think this fixes it:  Bryan
		if (--_me == (_end - (_s2-_s1+1))) { _me -= _s1; _end -= _s2; }		
		return *this;
	}

	// untested
	SubMatIterator operator--(int) {
		SubMatIterator p = *this; --*this; return p;
	} 

	bool operator==(const SubMatIterator&  b) const { return _me == b._me; 
	}

	bool operator!=(const SubMatIterator&  b) const { return _me != b._me; 
	}
}; // SubMatIterator

template<typename E>
struct ConstSubMatIterator { // periodic gaps iterator
	E * _me, * _end;
	size_t _s1, _s2;

	ConstSubMatIterator () : _me(0), _end(0), _s1(0), _s2(0) {
	}
	ConstSubMatIterator (E * loc, size_t c = 0, size_t stride = 0) : _me(loc), _end(loc+c), _s1(stride - c), _s2(stride) {
	}

	ConstSubMatIterator(const ConstSubMatIterator& p) : _me(p._me), _end(p._end), _s1(p._s1), _s2(p._s2) {
	}

	const E& operator*() const { 
		return *_me; 
	}

	ConstSubMatIterator& operator++() { 
		++_me; if (_me == _end) { _me += _s1;  _end += _s2; } return *this; 
	}

	bool operator==(const ConstSubMatIterator&  b) const { 
		return _me == b._me; 
	}

	bool operator!=(const ConstSubMatIterator&  b) const { 
		return _me != b._me; 
	}

}; // ConstSubMatIterator

/*
std::ostream& operator<< (std::ostream &out, const SubMatIterator& smi){
	smi.pinfo();
}
*/

}; // LinBox

#endif //__LINBOX_sliced3_SMI_H
