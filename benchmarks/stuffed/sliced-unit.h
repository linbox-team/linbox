template <typename T, unsigned int N, unsigned int D>
class SlicedUnit
{
	public: // typedefs
		using Word = CompressedWord<T,N,1>;
		using Base = typename Word::Base;
	
	public: // static helpers
		static unsigned int get_depth();
	
	public: // data
		std::array<Word,D> w;
	
	public: // constructors
		SlicedUnit () = default;
		SlicedUnit (const SlicedUnit &other) = default;
		~SlicedUnit () = default;
		
	
	public: // operators
		SlicedUnit& operator =   (const SlicedUnit &rhs) = default;
	
		SlicedUnit  operator &   (const Word rhs) const;
		SlicedUnit& operator &=  (const Word rhs);
		SlicedUnit  operator &   (const Base rhs) const;
		SlicedUnit& operator &=  (const Base rhs);
		
		SlicedUnit  operator |   (const SlicedUnit rhs) const;
		SlicedUnit& operator |=  (const SlicedUnit rhs);
		
		SlicedUnit  operator <<  (const unsigned int rhs) const;
		SlicedUnit& operator <<= (const unsigned int rhs);
		SlicedUnit  operator >>  (const unsigned int rhs) const;
		SlicedUnit& operator >>= (const unsigned int rhs);
		
		SlicedUnit  operator ~   () const;
		
		bool operator == (const SlicedUnit rhs) const;
		bool operator != (const SlicedUnit rhs) const;
	
	public: // operations
		void clear  ();
		bool isZero () const;
						
		SlicedUnit& negin  ();				
						
		SlicedUnit  nand   (const Word rhs) const;
		SlicedUnit& nandin (const Word rhs);
		SlicedUnit  nand   (const Base rhs) const;
		SlicedUnit& nandin (const Base rhs);
	
		SlicedUnit  lshift   (const unsigned int rhs) const; // element-wise
		SlicedUnit  rshift   (const unsigned int rhs) const;
		SlicedUnit& lshiftin (const unsigned int rhs);
		SlicedUnit& rshiftin (const unsigned int rhs);
						
		SlicedUnit  select     (unsigned int idx, unsigned int l) const;
		SlicedUnit& selectin   (unsigned int idx, unsigned int l);
		SlicedUnit  deselect   (unsigned int idx, unsigned int l) const;
		SlicedUnit& deselectin (unsigned int idx, unsigned int l);
		
		void swapEntry (unsigned int i, unsigned int j);
		void swapWord  (unsigned int i, unsigned int j);
	
	public: // debug
	
};