#ifndef __SLICED_STEPPER_H
#define __SLICED_STEPPER_H

#include "dense-sliced.h"
#include "sliced-domain.h"
#include <linbox/field/modular.h>

namespace LinBox{

// To fill a sliced vector
struct stepper {

	typedef SlicedDomain<LinBox::Modular<uint8_t> > Domain;
	typedef Sliced<Domain > Matrix;
	typedef Matrix::Scalar Scalar;
	typedef Matrix::RawIterator RawIterator;
	typedef Matrix::SlicedUnit SlicedUnit;

	Matrix &_A;
	RawIterator _r;
	size_t _i; 
	const static size_t SIZE = sizeof(Domain::Word_T)*8;
	Scalar _store[SIZE];

	stepper(Matrix & A) : _A(A), _r(_A.rawBegin()), _i(SIZE-1) { }

	virtual inline void flush(){
			//_r.pinfo();
			//std::cerr << std::endl;
			SlicedUnit &t = (*_r);
			t.zero();
			for(size_t i = _i + 1; i < SIZE; ++i){
				//std::cerr << (int)_store[i] << " ";
				/*
				t <<= 1;
				if(_store[i] > 1) t |= 1;
				else t.b0 |= _store[i];
				*/
				t <<= 1;
				t.b1 |= ((_store[i] & 2) >> 1); 
				t.b0 |= ((_store[i] & 1) | t.b1);
				//std::cerr << (int)t.b0 <<  " & " << (int)t.b1 << std::endl;
			} 
			//std::cerr << "----------------" << std::endl;
			_i = SIZE-1;
			++_r;
	}

	inline void step(Scalar e) { 
		//std::cerr << (int)e << " ";
		_store[_i--] = e;
		if(_i > SIZE) flush();
	}

	inline void row(){
		if(_i != SIZE-1) flush();
	}

};

} // LinBox

#endif // __SLICED_STEPPER_H
