#include "field-reader-analyzer.h"

#include "linbox-config.h"

// include just about every field in LinBox :-)
#include "linbox/field/gf2.h"
#include "linbox/field/modular.h"
#include "linbox/field/unparametric.h"
#include "linbox/field/gmp-rational.h"

#ifdef HAVE_NTL
#include "linbox/field/ntl.h"
#endif

#ifdef HAVE_GIVARO
#include "linbox/field/givaro.h"
#endif

#ifdef HAVE_LIDIA
#include "linbox/field/lidia.h"
#endif


namespace LinBox {



	// Default FieldReaderAnalyzer constructor.  Creats an empty Reader
	// and sets all fields to their "false" values
	FieldReaderAnalyzer::FieldReaderAnalyzer() {
		_isField = false;
		_fieldType = FieldReaderAnalyzer::NotField;
	}

	// Main FieldReaderAnalyzer constructor.  This is essentially 
	// equivalent to calling the "reset" method, so just call that
	FieldReaderAnalyzer::FieldReaderAnalyzer(Reader &R) {
		reset(R);
	}

	// reset - Takes in a Reader and "resets" the analyzer to use the
	// new field.  This method "clones" the Reader and resets all the 
	// flags.  It returns true or false based on whether the Reader
	// actually contains a field
	bool FieldReaderAnalyzer::reset(Reader &R) {
		long e;
		string s;
		_R = R;

		if(!_R.expectTagName("field") ||!_R.expectAttributeNum("cardinality", _card) || !_R.expectChildTag()) {
			_isField = false;
			_fieldType = FieldReaderAnalyzer::NotField;
			_card = integer(-1);
			return false;
		}
		_isField = true;
		if(!_R.checkAttributeString("implDetail", _implDetail))
			_implDetail = s;
		
		// now check the Field type
		_R.traverseChild();
		if(_R.checkTagName("finite")) {
			if(!_R.expectChildTag()) {
				// uh-oh, mal-formed field
				_isField = false;
				_fieldType = FieldReaderAnalyzer::NotField;
				_card = integer(-1);
				return false;
			}
			_R.traverseChild();
			if(!_R.expectTagName("characteristic")) {
					// uh-oh, mal-formed field
				_isField = false;
				_fieldType = FieldReaderAnalyzer::NotField;
				_card = integer(-1);
				return false;
			}
			_R.upToParent();
			// check for extended field
			if(_R.getNextChild()) { 
				if(!_R.expectChildTag()) {
					// uh-oh, mal-formed field
					_isField = false;
					_fieldType = FieldReaderAnalyzer::NotField;
					_card = integer(-1);
						return false;
				}
				_R.traverseChild();
				if(!_R.expectTagName("extension") || !_R.expectChildTag()) {
					// uh-oh, mal-formed field
					_isField = false;
					_fieldType = FieldReaderAnalyzer::NotField;
					_card = integer(-1);
					return false;
				}
				_R.traverseChild();
				if(!_R.expectTagNum(e)) {
					// uh-oh, mal-formed field
					_isField = false;
					_fieldType = FieldReaderAnalyzer::NotField;
					_card = integer(-1);
					return false;
				}
				if(e > 1) {
					_fieldType = FieldReaderAnalyzer::FiniteExtended;
				}
				else
					_fieldType = FieldReaderAnalyzer::Finite;
				_R.Up(2);
				_R.Left(1);
			}
			else
				_fieldType = FieldReaderAnalyzer::Finite;
		}
		else if(_R.checkTagName("rational")) {
			_fieldType = FieldReaderAnalyzer::Rational;
		}
		else if(_R.checkTagName("real")) {
			_fieldType = FieldReaderAnalyzer::Real;
		}
		else if(_R.checkTagName("integer")) {
			_fieldType = FieldReaderAnalyzer::Integer;
		}
		else 
			_fieldType = FieldReaderAnalyzer::Unknown;
		
	
		_R.Up(1);
		return true;
	}
					
	// isField - Returns true based on whether or not the contained Reader
	// is actually a field.  This is established when the Reader is first
	// given to the Analyzer, so it just returns a boolean value
	bool FieldReaderAnalyzer::isField() const {
		return _isField;
	}
			
	// whatType - Returns the generic type code of the object.  This is
	// again set when a reader is supplied with the reset method, so
	// this just returns a peice of state
	int FieldReaderAnalyzer::whatType() const {
		return _fieldType;
	}

	// implDetail - Returns the value of the optional _implDetail
	// attribute of the object.  If there is no _implDetail, the empty 
	// string (the value of the STL string class after the default
	// constructor was called) is returned
	//
	string FieldReaderAnalyzer::implDetail() const {
		return _implDetail;
	}

	// cardinality - Returns the cardinality of the field.  Returns
	// 0 if the field is of infinite cardinality and -1 if the 
	// object isn't a field
	//
	integer FieldReaderAnalyzer::cardinality() const {
		return _card;
	}

	// makeField - The first of two.  This method is templatized on
	// the type of the functor and the type of the UserData data
	// structure to be passed into the functor as an auxilary 
	// parameter.
	// This method will create an instance of the field described by
	// the Reader and calls the functor on the pointer to this object
	// Note, in any case but a FiniteExtended Field, LinBox has SOME
	// kind of default class it can fall back to (even if it is 
	// something like Unparametric<double> for real)
	// Note, if the Reader doesn't describe a Field, the exception
	// ReaderNotAField is thrown
	//
	// This method works as follows:  The type of the field is first
	// case analyzed using a case statement.  If the implDetail was
	// available, attempt to match this to a field in LinBox.  If
	// this is impossible (which can happen if the original field was
	// derived from some optional plug-in package like Givaro, but this
	// package isn't available on the system used to build this Field,
	// the default for each type is used (at current, there is no
	// FiniteExtended Field in LinBox without a plug-in package), so 
	// the default in this case is to throw "NoFiniteExtendedInLinBox"
	//
	template<class FunctorType, class DSType>
	void* FieldReaderAnalyzer::makeField(const FunctorType &FT, DSType &AuxUserData) const
	{       

		switch(whatType()) {

		// For this type, there are the follwing Native
		// LinBox types:
		//   - GF2
		//   - Modular<type> where type is one of
		//          - uint8
		//          - uint16
		//          - uint32
		//          - integer (This is the default for both
		//                     Modular & for all Finite fields in LB) 
		//
		// There are also the following types that come from plug-in
		// packages:
		//   - From NTL:
		//      - NTL::ZZ_p (defaults to modular<integer>)
		//      - NTL::zz_p (this one defaults to modular<uint32>)
		//      
		//   - From Givaro:
		//      - GivaroZpz (defaults to modular<integer>)
		//          The above is templatized on either Std16, Std32,
		//          or Log16.  Detect each in the _implDetail attrib
		//          Also, as LinBox has no native Table based finite,
		//          there is no default, so throw the exception
		//          NoTableFiniteInLinBox 
			
		case FieldReaderAnalyzer::Finite: {

			if(_implDetail == "gf2") {
				GF2* p = new GF2(_R);
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "modular-uint8") {
				Modular<uint8>* p = new Modular<uint8>(_R);
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "modular-uint16") {
				Modular<uint16>* p = new Modular<uint16>(_R);
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "modular-uint32") {
				Modular<uint32>* p = new Modular<uint32>(_R);
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "modular-integer") {
				Modular<integer> *p = new Modular<integer>(_R);
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "modular") {
				Modular<integer> *p = new Modular<integer>(_R);
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "ntl-zzp") {
#ifdef HAVE_NTL // if we have NTL
				UnparametricField<NTL::zz_p> *p = new UnparametricField<NTL::zz_p>(_R);
#else // the default
				Modular<uint32> *p = new Modular<uint32>(_R);
#endif
				return FT(p, AuxUserData);

			}
			else if(_implDetail == "ntl-ZZp") {
#ifdef HAVE_NTL // if we have NTL
				UnparametricField<NTL::ZZ_p> *p = new UnparametricField<NTL::ZZ_p>(_R);
#else
				Modular<integer> *p = new Modular<integer>(_R);
#endif
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "givaro-zpz-std16") {
#ifdef HAVE_GIVARO // if we have Givaro
				GivaroZpz<Std16> *p = new GivaroZpz<Std16>(_R);
#else
				Modular<uint16> *p = new Modular<uint16>(_R);
#endif
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "givaro-zpz-std32") {
#ifdef HAVE_GIVARO
				GivaroZpz<Std32> *p = new GivaroZpz<Std32>(_R);
#else
				Modular<uint32> *p = new Modular<uint32>(_R);
#endif
				return FT(p, AuxUserData);
			}
			else if(_implDetail == "givaro-zpz-log16") {
#ifdef HAVE_GIVARO
				GivaroZpz<Log16> *p = new GivaroZpz<Log16>(_R);
#else
				// crap, we're screwed, just throw the
				// exception and bail
				throw NoTableFiniteInLinBox;
#endif
				return FT(p, AuxUserData);
			}
			else { // we take our chances w/ Modular<integer>
				Modular<integer> *p = new Modular<integer>(_R);
				return FT(p, AuxUserData);
			}
			break;
		}

		// LinBox doesn't have any native support for FiniteExtension
		// fields, so if there is nothing, the  NoFiniteExtenedInLinBox
		// exception is thrown.  The follwing plug-in packages
		// are supported:
		//  - From NTL - NTL::ZZ_pE
		//  - From Givaro - GivaroGfq
		//  - From Lidia - LidiaGfq
		// By default if there is nothing to fall back on, so if
		// there isn't a field to match _implDetail, try to use any 
		// field.  If that doesn't work, just throw the error
			
		case FieldReaderAnalyzer::FiniteExtended : {

			if(_implDetail == "ntl-ZZpE") {
#ifdef HAVE_NTL
				UnparametricField<NTL::ZZ_pE>* p = new UnparametricField<NTL::ZZ_pE>(_R);
				return FT(p, AuxUserData);
#endif // if we don't have the package, fall through to the 
				// "we'll use anything!" stage
			}
			else if(_implDetail == "givaro-gfq") {
#ifdef HAVE_GIVARO
				GivaroGfq *p = new GivaroGfq(_R);
				return FT(p, AuxUserData);
#endif
			}
			else if(_implDetail == "lidia-gfq") {
#ifdef HAVE_LIDIA
				LidiaGfq *p = new LidiaGfq(_R);
				return FT(p, AuxUserData);
#endif
			}
			// Now we are in desperation mode, just use the first
			// plug-in package that we've got
			//
			else {
#ifdef HAVE_GIVARO
				GivaroGfq *p = new GivaroGfq(_R);
				return FT(p, AuxUserData);
#elif HAVE_LIDIA
				LidiaGfq *p = new LidiaGfq(_R);
				return FT(p, AuxUserData);
#elif HAVE_NTL // Our last resort
				UnparametricField<NTL::ZZ_pE> *p = new UnparametricField<NTL::ZZ_pE>(_R);
				return FT(p, AuxUserData);
#else // uh-oh, we're cooked
				throw NoFiniteExtendedInLinBox;
#endif
			}
			break;
		}
		// In this case there is only one rational class in
		// LinBox, GMPRationalField, so just use it

		case FieldReaderAnalyzer::Rational : {
			GMPRationalField *p = new GMPRationalField(_R);
			return FT(p, AuxUserData);
			break;
		}

		// For this type, there are the following Native LinBox
		// types: UnparametricField<double> (not a good one)
		// From plug-in packages, there is NTL, with NTL::RR
		// If NTL is not present, use UnparametricField
		case FieldReaderAnalyzer::Real : {

			// There's no need to check implDetail in this case,
			// because there's only one field which properly
			// produces a real field, NTL::RR.  Let's hope we
			// have it
#ifdef HAVE_NTL
			UnparametricField<NTL::RR>* p = new UnparametricField<NTL::RR>(_R);
#else
			UnparametricField<double>* p = new UnparametricField<double>(_R);
#endif
			return FT(p, AuxUserData);
			break;
		}

		// Integer case.  At present there is no "integer" field, it
		// is assumed that calculations will be carried out in the 
		// integer domain by reducing the field mod a group of primes
		// then building the result back up using the CRT
		// In any event, this functionality has yet to be added to
		// LinBox, so just use UnparametricField<integer>
		case FieldReaderAnalyzer::Integer : {

			UnparametricField<integer> *p = new UnparametricField<integer>(_R);
			return FT(p, AuxUserData);
			break;
		}

		// We don't know what it is, but this is the best we can approximate
		// go with it, and try.  If the compiler can't sort it out, there's no chance
		// so we need to add a new type.
		//
		case FieldReaderAnalyzer::Unknown : {
			UnparametricField<integer> *p = new UnparametricField<integer>(_R);
			return FT(p, AuxUserData);
			break;
		}	

		case FieldReaderAnalyzer::NotField:
		default: 
			throw ReaderNotAField();
			break;
		}

		// we never get here, but for good measure
		return NULL;
	}


	// makeField - The second of two.  This one is exactly like the
	// first above, except that it essentially doesn't expect a user
	// supplied Data Structure (so as a result it doesn't call the
	// functor with such a Data Type, doesn't take a user data structure,
	// and isn't templatized for a user data structure
	//
	template<class FunctorType>
	void* FieldReaderAnalyzer::makeField(const FunctorType &FT) const
	{       

		switch(whatType()) {

		// For this type, there are the follwing Native
		// LinBox types:
		//   - GF2
		//   - Modular<type> where type is one of
		//          - uint8
		//          - uint16
		//          - uint32
		//          - integer (This is the default for both
		//                     Modular & for all Finite fields in LB) 
		//
		// There are also the following types that come from plug-in
		// packages:
		//   - From NTL:
		//      - NTL::ZZ_p (defaults to modular<integer>)
		//      - NTL::zz_p (this one defaults to modular<uint32>)
		//      
		//   - From Givaro:
		//      - GivaroZpz (defaults to modular<integer>)
		//          The above is templatized on either Std16, Std32,
		//          or Log16.  Detect each in the _implDetail attrib
		//          Also, as LinBox has no native Table based finite,
		//          there is no default, so throw the exception
		//          NoTableFiniteInLinBox 
			
		case FieldReaderAnalyzer::Finite: {

			if(_implDetail == "gf2") {
				GF2* p = new GF2(_R);
				return FT(p);
			}
			else if(_implDetail == "modular-uint8") {
				Modular<uint8>* p = new Modular<uint8>(_R);
				return FT(p);
			}
			else if(_implDetail == "modular-uint16") {
				Modular<uint16>* p = new Modular<uint16>(_R);
				return FT(p);
			}
			else if(_implDetail == "modular-uint32") {
				Modular<uint32>* p = new Modular<uint32>(_R);
				return FT(p);
			}
			else if(_implDetail == "modular-integer") {
				Modular<integer> *p = new Modular<integer>(_R);
				return FT(p);
			}
			else if(_implDetail == "modular") {
				Modular<integer> *p = new Modular<integer>(_R);
				return FT(p);
			}
			else if(_implDetail == "ntl-zzp") {
#ifdef HAVE_NTL // if we have NTL
				UnparametricField<NTL::zz_p> *p = new UnparametricField<NTL::zz_p>(_R);
#else // the default
				Modular<uint32> *p = new Modular<uint32>(_R);
#endif
				return FT(p);

			}
			else if(_implDetail == "ntl-ZZp") {
#ifdef HAVE_NTL // if we have NTL
				UnparametricField<NTL::ZZ_p> *p = new UnparametricField<NTL::ZZ_p>(_R);
#else
				Modular<integer> *p = new Modular<integer>(_R);
#endif
				return FT(p);
			}
			else if(_implDetail == "givaro-zpz-std16") {
#ifdef HAVE_GIVARO // if we have Givaro
				GivaroZpz<Std16> *p = new GivaroZpz<Std16>(_R);
#else
				Modular<uint16> *p = new Modular<uint16>(_R);
#endif
				return FT(p);
			}
			else if(_implDetail == "givaro-zpz-std32") {
#ifdef HAVE_GIVARO
				GivaroZpz<Std32> *p = new GivaroZpz<Std32>(_R);
#else
				Modular<uint32> *p = new Modular<uint32>(_R);
#endif
				return FT(p);
			}
			else if(_implDetail == "givaro-zpz-log16") {
#ifdef HAVE_GIVARO
				GivaroZpz<Log16> *p = new GivaroZpz<Log16>(_R);
#else
				// crap, we're screwed, just throw the
				// exception and bail
				throw NoTableFiniteInLinBox;
#endif
				return FT(p);
			}
			else { // we take our chances w/ Modular<integer>
				Modular<integer> *p = new Modular<integer>(_R);
				return FT(p);
			}
			break;
		}

		// LinBox doesn't have any native support for FiniteExtension
		// fields, so if there is nothing, the  NoFiniteExtenedInLinBox
		// exception is thrown.  The follwing plug-in packages
		// are supported:
		//  - From NTL - NTL::ZZ_pE
		//  - From Givaro - GivaroGfq
		//  - From Lidia - LidiaGfq
		// By default if there is nothing to fall back on, so if
		// there isn't a field to match _implDetail, try to use any 
		// field.  If that doesn't work, just throw the error
			
		case FieldReaderAnalyzer::FiniteExtended : {

			if(_implDetail == "ntl-ZZpE") {
#ifdef HAVE_NTL
				UnparametricField<NTL::ZZ_pE>* p = new UnparametricField<NTL::ZZ_pE>(_R);
				return FT(p);
#endif // if we don't have the package, fall through to the 
				// "we'll use anything!" stage
			}
			else if(_implDetail == "givaro-gfq") {
#ifdef HAVE_GIVARO
				GivaroGfq *p = new GivaroGfq(_R);
				return FT(p);
#endif
			}
			else if(_implDetail == "lidia-gfq") {
#ifdef HAVE_LIDIA
				LidiaGfq *p = new LidiaGfq(_R);
				return FT(p);
#endif
			}
			// Now we are in desperation mode, just use the first
			// plug-in package that we've got
			//
			else {
#ifdef HAVE_GIVARO
				GivaroGfq *p = new GivaroGfq(_R);
				return FT(p);
#elif HAVE_LIDIA
				LidiaGfq *p = new LidiaGfq(_R);
				return FT(p);
#elif HAVE_NTL // Our last resort
				UnparametricField<NTL::ZZ_pE> *p = new UnparametricField<NTL::ZZ_pE>(_R);
				return FT(p);
#else // uh-oh, we're cooked
				throw NoFiniteExtendedInLinBox;
#endif
			}
			break;
		}

		// In this case there is only one rational class in
		// LinBox, GMPRationalField, so just use it

		case FieldReaderAnalyzer::Rational : {
			GMPRationalField *p = new GMPRationalField(_R);
			return FT(p);
			break;
		}

		// For this type, there are the following Native LinBox
		// types: UnparametricField<double> (not a good one)
		// From plug-in packages, there is NTL, with NTL::RR
		// If NTL is not present, use UnparametricField
		case FieldReaderAnalyzer::Real : {

			// There's no need to check implDetail in this case,
			// because there's only one field which properly
			// produces a real field, NTL::RR.  Let's hope we
			// have it
#ifdef HAVE_NTL
			UnparametricField<NTL::RR>* p = new UnparametricField<NTL::RR>(_R);
#else
			UnparametricField<double>* p = new UnparametricField<double>(_R);
#endif
			return FT(p);
			break;
		}

		// Integer case.  At present there is no "integer" field, it
		// is assumed that calculations will be carried out in the 
		// integer domain by reducing the field mod a group of primes
		// then building the result back up using the CRT
		// In any event, this functionality has yet to be added to
		// LinBox, so just use UnparametricField<integer>
		case FieldReaderAnalyzer::Integer : {

			UnparametricField<integer> *p = new UnparametricField<integer>(_R);
			return FT(p);
			break;
		}

		// we don't know what type we have, this is an anonymous type.  This is the best I can
		// do to recover without dying.  If this doesn't work, a new type should definately be
		// added. 
		// 	
		case FieldReaderAnalyzer::Unknown: {
			UnparametricField<integer> *p = new UnparametricField<integer>(_R);
			return FT(p);
			break;
		}


		case FieldReaderAnalyzer::NotField:
		default:
			throw ReaderNotAField();
			break;
		}

		// we never get here, but for good measure
		return NULL;
	}

} // namespace LinBox {

#endif // #ifndef __FIELD_READER_ANALYZER
