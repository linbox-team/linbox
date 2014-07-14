#ifndef SLICED_INL
#define SLICED_INL
// This is designed to be included in the middle of the Sliced class declaration.

/*
	void splitOp(SlicedUnit &a, op_t t, bool first, int offset, int arg = 0, SlicedUnit &b = 0);
	Sliced& s_addin(Sliced &other);
	Sliced& s_smulin(Scalar &x);
	Sliced & s_axpyin(RawIterator &b, RawIterator &e, Scalar &s, RawIterator &ob);
*/

	//  addin, axpyin, 
	void compatible(){
		
	}

	//  if we're a submatrix, we must not touch bits outside
	//  our boundary if we begin or end in the middle of a sliced word
	//  this helper function will do the operation
	//  a is the input/output sliced word.  t is the operation
	//  "last" is false if [0,offset] part of word, true if [offset,_SIZE-1]
	//  arg is in case of smul or axpy, the scalar multiple
	//  b is the other word in case of a binary operation
	void splitOp(SlicedUnit &a, op_t t, bool last, int offset, Scalar &arg, SlicedUnit &b){
		//std::cerr << "DSINL" << std::endl;

		SlicedWord computeMask;
		SlicedWord one = (SlicedWord)1;

		if(last)
			computeMask = (one << offset) - 1; // 2^offset - 1
		else
			computeMask = ~((one << (_SIZE-offset)) - 1);

		/*
		cerr << "arg: " << arg << endl;
		cerr << "COMPUTEMASK: " << computeMask << endl;
		cerr << "0: " << a.b0 << "  1: " << a.b1 << endl;
		*/

		SlicedUnit save = a & (~computeMask);

		switch(t){
			case AXPY:
				if(arg == 2){
					a += b*arg;
					break;
				}
				//  else arg will be 1, fall to add
			case ADD:
				a += b;
				break;
			case SMUL:
				if(arg){
					a *= arg;
					break;
				}
				//  else arg will be 0, fall to zero
			case ZERO:
				a.zero();
				break;
			case COPY:
				a = b;
				break;
		}

		a &= computeMask; // zero out bits that aren't "ours"
		a |= save;  //  retrieve bits that we "never touched"
	}

	//  specialized addin, care must be taken at beg. and end of rows
	Sliced& s_addin(Sliced &other){
		RawIterator r = rawBegin();
		RawIterator o = other.rawBegin();
		Scalar t;

		size_t i, j;
		//  loop all rows
		for(i=0; i<rows(); ++i){
			j = 0;
			//  do first elt
			if(_loff){
				splitOp(*r, ADD, false, _loff, t, *o);
				++r; ++o;
				j+=_loff;
			}
			//  do "middle" elts
			for(; j<cols()-_roff; j+=_SIZE, ++r, ++o)
				(*r) += (*o);
			//  do "last" elt
			if(_roff){
				splitOp(*r, ADD, true, _roff, t, *o);
				++r; ++o;
			}
		}
		return *this;
	}

	//  specialized smulin, care must be taken at beg. and end of rows
	Sliced& s_smulin(Scalar &x){
		RawIterator r = rawBegin();
		size_t i, j;
		//  loop all rows
		for(i=0; i<rows(); ++i){
			j = 0;
			//  do first elt
			if(_loff){
				splitOp(*r, SMUL, false, _loff, x, *r);
				++r;
				j+=_loff;
			}
			//  do "middle" elts
			if(x == 0){
				for(; j<cols()-_roff; j+=_SIZE, ++r){
					//r.pinfo();
					(*r).zero();
				}
			}
			else{
				for(; j<cols()-_roff; j+=_SIZE, ++r){
					(*r) *= x;
				}
			}
			//  do "last" elt
			if(_roff){
				//cerr << "ROFF" <<  _roff << endl;
				//r.pinfo();
				splitOp(*r, SMUL, true, _roff, x, *r);
				++r;
			}
		}
		return *this;
	} 

	//  arguments: begin iter, scalar, other mat's begin iter
	Sliced & s_axpyin(RawIterator &r, Scalar &s, RawIterator &o){
		size_t j = 0;
		//  NO ROWS, WE ASSUME TO WORK ON A SINGLE ROW
			j = 0;
			//  do first elt
			if(_loff){
				splitOp(*r, AXPY, false, _loff, s, *o);
				++r; ++o;
				j+=_loff;
			}
			//  do "middle" elts
			if(s == 1){
				for(; j<cols()-_roff; j+=_SIZE, ++r, ++o)
					(*r) += (*o);
			}
			else if(s == 2){
				for(; j<cols()-_roff; j+=_SIZE, ++r, ++o)
					(*r) += (*o)*2;
			}
			//  do "last" elt
			if(_roff)
				splitOp(*r, ADD, true, _roff, s, *o);
		//  END OF (NO) ROWS
		return *this;
	}

	//  specialized copy, care must be taken at beg. and end of rows
	Sliced & s_copy(Sliced &other){
		RawIterator r = rawBegin();
		RawIterator o = other.rawBegin();

		Scalar t;

		/* for an undeveloped,
			more complicated but likely faster
			memcpy for middle of rows
		int midSection = ((cols()-(_loff+_roff))+_SIZE-1)/_SIZE;
		int bytes = midSection * 2 * sizeof(SlicedWord);
		int data_offset = 0;
		*/

		size_t i, j;
		//  loop all rows
		for(i=0; i<rows(); ++i){
			j = 0;
			//  do first elt
			if(_loff){
				splitOp(*r, COPY, false, _loff, t, *o);
				++r; ++o;
				j+=_loff;
			}
			//  do "middle" elts // TODO could be memcpy -> faster
			for(; j<cols()-_roff; j+=_SIZE, ++r, ++o)
				(*r) = (*o);
			//  do "last" elt
			if(_roff){
				splitOp(*r, COPY, true, _roff, t, *o);
				++r; ++o;
			}
		}
		return *this;
	}

	//  IF THERE IS A LEFT OFFSET (RIGHT OFFSET POSSIBLE, TOO)
	//  specialized copy, care must be taken at beg. and end of rows
	std::ostream & s_write_bin(std::ostream &os){
		// how many bits do we have room for?
		int A = _SIZE - _loff; 
		int spill = _roff ? _roff - A : _SIZE - A;

		SlicedWord one = (SlicedWord)1;
		//  lower bits mask
		SlicedWord mask_lowbits = (one << A) - one;
		//  for last word, move _roff bits if _roff < A
		SlicedWord last_mask;
		last_mask = spill < 0 ? (one << _roff) - one : mask_lowbits;

		//  storing temp values
		SlicedUnit lower, upper;
		int size = 2 * sizeof(SlicedWord);

		//  for rare case of only one sliced-word per row (< _SIZE elts)
		bool onlyOne = false;  //  robusto

		RawIterator r, e;
		//  loop all rows
		for(size_t i=0; i<rows(); ++i){

			r = rowBegin(i);  // first in the row
			e = ((rowEnd(i))); // e is after the last in the row
			--e;  // e is now the last in the row
			if(r == e) onlyOne=true; //  if there's only one sliced word in this row
			else --e;  //  e will now be the 2nd-to-last in the row

			//  iterate until 2nd to last word:
			for( ; r!=e; ){
				//  get "end" of "current" word and make it our beginning
				lower = (*r);
				lower >>=A;
				//  move to next word
				++r;
				//  get "beginning" of "next" word and make it our end
				upper = (*r);
				(upper &= mask_lowbits) <<= _loff;
				lower |= upper;	 //  lower now has values we want
				os.write((char *)&lower, size);  // write to file
			}
			//  process the end of the row
			lower = (*r);
			lower >>= A;  
			//  no further data
			if(onlyOne){
				os.write((char *)&lower, size);
			}
			else{ // one more word to split up
				++r;
				upper = (*r);
				(upper &= last_mask) <<= _loff;
				lower |= upper;
				os.write((char *)&lower, size);
				//  process final word, if necessary
				if(spill > 0){
					lower = (*r);
					lower >>= A;
					//  mask to "our" bits
					if(_roff)
						lower &= ((one << spill) - one);
					os.write((char *)&lower, size); // write 2nd-to-final
				}
			}
		}
		return os;
	}

	//  IF THERE IS A RIGHT-OFFSET ONLY
	//  specialized copy, care must only be taken at end of rows
	std::ostream & s_write_bin_r(std::ostream &os){
		//  TODO  - I can't just declare these
		RawIterator r = rowBegin(0); RawIterator e = rowEnd(0);

		SlicedWord one = (SlicedWord)1;
		//  lower bits mask
		SlicedWord right_mask = (one << _roff) - one;

		//  storing last word's temp value
		SlicedUnit final;

		//  size of a sliced unit
		int size = 2 * sizeof(SlicedWord);
		//  accounts for all but the last sliced unit (the -1)
		int bytesPerRow = (((_n + _SIZE-1)/_SIZE) - 1) * size;

		//  loop all rows
		for(size_t i=0; i<rows(); ++i){
			//  set up pointers
			r = rowBegin(i);
			e = (--(rowEnd(i)));  // the last sliced unit

			//  block write all but the last sliced unit:
			os.write((char *)&(*r), bytesPerRow); 

			//  process the last unit:
			final = (*e);
			final &= right_mask;  //  cut down to size
			os.write((char *)&final, size);
		}
		return os;
	}

   //  if we're a submatrix
   std::ostream & s_wb(std::ostream &os){
      if(_loff) return s_write_bin(os);
      else if(_roff) return s_write_bin_r(os);

      //  else we're totally aligned, 
      //  safe to write row at a time
      int bytesPerRow = ((_n + _SIZE-1)/_SIZE) * 2 * sizeof(SlicedWord);

      for(size_t i = 0; i < _m; ++i){
         os.write((char *)(&*rowBegin(i)), bytesPerRow);
      }

      return os;
   }

	//  ------------------------------------------------------------
	//  TODO: there is a major major hack here
	//  TODO: and it's serious.
	//  right now we're assuming we WANT zeroes in the last word
	//  there is NO MAINTAINING of the right half of the last word
	//  the way the fn is written.

	//  IF THERE IS A RIGHT-OFFSET ONLY
	//  specialized copy, care must only be taken at end of rows
	std::istream & s_read_bin_r(std::istream &is){
		//std::cerr << "Well, here we are..." << std::endl;
		//  TODO  - I can't just declare these
		RawIterator r = rowBegin(0); RawIterator e = rowEnd(0);

		//SlicedWord one = (SlicedWord)1;
		//  lower bits mask
		//SlicedWord right_mask = (one << _roff) - one;
		//  storing last word's temp value
		//SlicedUnit final;

		//  size of a sliced unit
		int size = 2 * sizeof(SlicedWord);
		//  accounts for all but the last sliced unit (the -1)
		int bytesPerRow = (((_n + _SIZE-1)/_SIZE) - 1) * size;

		//  loop all rows
		for(size_t i=0; i<rows(); ++i){
			//  set up pointers
			r = rowBegin(i);
			e = (--(rowEnd(i)));  // the last sliced unit

			//  block read all but the last sliced unit:
			is.read((char *)&(*r), bytesPerRow); 

			//  read the last one (should have zeros automatically
			is.read((char *)&(*e), size);

			// should have zeroes automatically... 
			//  following lines (form above) not needed

			//  process the last unit:
			//final = (*e);
			//final &= right_mask;  //  cut down to size
			//os.write((char *)&final, size);
		}
		return is;
	}

	// if we're a submatrix
	std::istream& s_rb(std::istream &is){
		if(_loff){ std::cerr << "MEGA DEATH" << std::endl; exit(-1); }
		if(_roff) return s_read_bin_r(is);

		//  else we're totally aligned, 
		//  safe to write row at a time
		int bytesPerRow = ((_n + _SIZE-1)/_SIZE) * 2 * sizeof(SlicedWord);

		for(size_t i = 0; i < _m; ++i){
			is.read((char *)(&*rowBegin(i)), bytesPerRow);
		}

		return is;
	}

/* //  MUL:
	Sliced & s_mul(Sliced& A, Sliced& B){
		return *this;
	}
*/

#endif // __SLICED_INL
