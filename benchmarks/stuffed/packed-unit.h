template <typename T, unsigned int N, unsigned int R, unsigned int L>
class PackedUnit
{
	public: // typedefs
		using Word = CompressedWord<T,N,R+L>;
		using Base = typename Word::Base;
	
	public: // static helpers
		static unsigned int get_r();
		static unsigned int get_l();
		static Word r_mask();
		static Word l_mask();
		static Word word_entry_mask();
		static Base base_entry_mask();
	
	public: // data
		std::array<Word,1> w;
	
	public: // constructors
		PackedUnit () = default;
		PackedUnit (const PackedUnit &other) = default;
		~PackedUnit () = default;
		
	
	public: // operators
		PackedUnit& operator =   (const PackedUnit &rhs) = default;
	
		PackedUnit  operator &   (const Word rhs) const;
		PackedUnit& operator &=  (const Word rhs);
		PackedUnit  operator &   (const Base rhs) const;
		PackedUnit& operator &=  (const Base rhs);
		
		PackedUnit  operator |   (const PackedUnit rhs) const;
		PackedUnit& operator |=  (const PackedUnit rhs);
		
		PackedUnit  operator <<  (const unsigned int rhs) const;
		PackedUnit& operator <<= (const unsigned int rhs);
		PackedUnit  operator >>  (const unsigned int rhs) const;
		PackedUnit& operator >>= (const unsigned int rhs);
		
		PackedUnit  operator ~   () const;
		
		bool operator == (const PackedUnit rhs) const;
		bool operator != (const PackedUnit rhs) const;
	
	public: // operations
		void clear  ();
		bool isZero () const;
						
		PackedUnit& negin  ();				
						
		PackedUnit  nand   (const Word rhs) const;
		PackedUnit& nandin (const Word rhs);
		PackedUnit  nand   (const Base rhs) const;
		PackedUnit& nandin (const Base rhs);
	
		PackedUnit  lshift   (const unsigned int rhs) const; // element-wise
		PackedUnit  rshift   (const unsigned int rhs) const;
		PackedUnit& lshiftin (const unsigned int rhs);
		PackedUnit& rshiftin (const unsigned int rhs);
						
		PackedUnit  select     (unsigned int idx, unsigned int l) const;
		PackedUnit& selectin   (unsigned int idx, unsigned int l);
		PackedUnit  deselect   (unsigned int idx, unsigned int l) const;
		PackedUnit& deselectin (unsigned int idx, unsigned int l);
		
		void swapEntry (unsigned int i, unsigned int j);
		void swapWord  (unsigned int i, unsigned int j);
	
	public: // debug
	
};