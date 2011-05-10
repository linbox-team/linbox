/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_field_reader_analyzer_H
#define __LINBOX_field_reader_analyzer_H

#include <string>

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/integer.h"


namespace LinBox
{

	// empty class definition for error throw below
	class ReaderNotAField {};
	class NoFiniteExtendedInLinBox {};
	class NoTableFiniteInLinBox {};

	class FieldReaderAnalyzer {

	public:
		FieldReaderAnalyzer();
		FieldReaderAnalyzer(Reader &);
		~FieldReaderAnalyzer() {}

		// This method resets the analyzer with a new reader
		bool reset(Reader &);

		// These are the "generic field types" that are supported in
		// LinBox.  Each one of them has a default (sort of, though
		// we can't really do FiniteExtended without NTL, Givaro or
		// Lidia, w/o NTL the best we can do for Real is
		// Unparametric<double>, though w/ integer we could do
		// Unparametric<integer>.  Note that of course NotField
		// doesn't represent a field but the fact that the
		// user screwed up
		//
		static const int Finite = 0;
		static const int FiniteExtended = 1;
		static const int Rational = 2;
		static const int Real = 3;
		static const int Integer = 4;
		static const int Unknown = 5;

		static const int NotField = -1;

		// This method returns true if the Reader last used to
		// reset the analyzer was actually a field
		bool isField() const;

		// This method returns one of the values above, depending
		// on which field was last used to set the analyzer
		int whatType() const;

		// This method returns the std::string of the implDetail attribute
		// of the field, if the object analyzed was a field.  If it
		// wasn't a field, return the empty std::string (ie, the std::string
		// created by "std::string s;"
		//
		std::string implDetail() const;

		// returns the field's cardinality if the field is actually
		// a field, returns -1 otherwise
		integer cardinality() const;


		// Okay, this is the complicated method.
		// This method tries to create an instance of the field
		// given in the Reader on the heap.  As arguments it takes
		// a functor templatized on the field type and a data
		// structure type, and an instance of this data structure
		// type to be passed as the second argument of the functor.
		// This DSType is meant to allow the user to pass some kind
		// of data structure in to their functor call
		//
		// Note, if the Reader given to be analyzed does not describe
		// a field (or the reader attempts to initalize the field
		// and fails for whatever reason), an exception is thrown
		// Another special note for the author of functors to
		// work in this thing:  it is expected that the operator()
		// of this functor (or the function itself, if it's a function)
		// be templatized on the field type and take a Field* as it's
		// first argument, and expect to take a reference to a
		// DSType as it's second argument, and to return a void*
		// It is also assumed that the responsibility for cleaning
		// up the field created is passed on to the functor (which
		// can in turn pass this responsibility on to whoever it likes)
		//
		// Sorry, I didn't mean to make this so complicated, but
		// I can think of 3 very different uses for this thing,
		// and short of implementing this in Scheme, this is the best
		// I can think of to make the FieldReaderAnalyzer work in
		// all of them
		//
		template<class FunctorType, class DSType>
		void* makeField(FunctorType &, DSType &) const;

		// This one is the same as the one above, but is designed
		// for functors that don't need another data argument to
		// be passed in.
		//
		template<class FunctorType>
		void* makeField(FunctorType &) const;


	private:
		mutable Reader _R;
		bool _isField;
		int _fieldType;
		std::string _implDetail;
		integer _card;
	};

}



// include just about every field in LinBox :-)
#include "linbox/linbox-config.h"

#include "linbox/field/gf2.h"
#include "linbox/field/modular.h"
#include "linbox/field/unparametric.h"
#include "linbox/field/gmp-rational.h"

#ifdef __LINBOX_HAVE_NTL
#include "linbox/field/ntl.h"
#endif

#if 0
#ifdef __LINBOX_HAVE_GIVARO
#include "linbox/field/givaro.h"
#endif

#ifdef __LINBOX_HAVE_LIDIA
#include "linbox/field/lidia.h"
#endif
#endif


namespace LinBox
{



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
		std::string s;
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
	// std::string (the value of the STL std::string class after the default
	// constructor was called) is returned
	//
	std::string FieldReaderAnalyzer::implDetail() const {
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
	void* FieldReaderAnalyzer::makeField(FunctorType &FT, DSType &AuxUserData) const
	{

		switch(whatType()) {

			// For this type, there are the follwing Native
			// LinBox types:
			//   - GF2
			//   - Modular<type> where type is one of
			//          - uint8_t
			//          - uint16_t
			//          - uint32_t
			//          - integer (This is the default for both
			//                     Modular & for all Finite fields in LB)
			//
			// There are also the following types that come from plug-in
			// packages:
			//   - From NTL:
			//      - NTL::ZZ_p (defaults to modular<integer>)
			//      - NTL::zz_p (this one defaults to modular<uint32_t>)
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
								  Modular<uint8_t>* p = new Modular<uint8_t>(_R);
								  return FT(p, AuxUserData);
							  }
							  else if(_implDetail == "modular-uint16") {
								  Modular<uint16_t>* p = new Modular<uint16_t>(_R);
								  return FT(p, AuxUserData);
							  }
							  else if(_implDetail == "modular-uint32_t") {
								  Modular<uint32_t>* p = new Modular<uint32_t>(_R);
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
#ifdef __LINBOX_HAVE_NTL // if we have NTL
								  UnparametricField<NTL::zz_p> *p = new UnparametricField<NTL::zz_p>(_R);
#else // the default
								  Modular<uint32_t> *p = new Modular<uint32_t>(_R);
#endif
								  return FT(p, AuxUserData);

							  }
							  else if(_implDetail == "ntl-ZZp") {
#ifdef __LINBOX_HAVE_NTL // if we have NTL
								  UnparametricField<NTL::ZZ_p> *p = new UnparametricField<NTL::ZZ_p>(_R);
#else
								  Modular<integer> *p = new Modular<integer>(_R);
#endif
								  return FT(p, AuxUserData);
							  }
#if 0
							  else if(_implDetail == "givaro-zpz-std16") {
#ifdef __LINBOX_HAVE_GIVARO // if we have Givaro
								  GivaroZpz<Std16> *p = new GivaroZpz<Std16>(_R);
#else
								  Modular<uint16_t> *p = new Modular<uint16_t>(_R);
#endif
								  return FT(p, AuxUserData);
							  }
							  else if(_implDetail == "givaro-zpz-std32") {
#ifdef __LINBOX_HAVE_GIVARO
								  GivaroZpz<Std32> *p = new GivaroZpz<Std32>(_R);
#else
								  Modular<uint32_t> *p = new Modular<uint32_t>(_R);
#endif
								  return FT(p, AuxUserData);
							  }
							  else if(_implDetail == "givaro-zpz-log16") {
#ifdef __LINBOX_HAVE_GIVARO
								  GivaroZpz<Log16> *p = new GivaroZpz<Log16>(_R);
#else
								  // crap, we're screwed, just throw the
								  // exception and bail
								  throw NoTableFiniteInLinBox;
#endif
								  return FT(p, AuxUserData);
							  }
#endif
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
#ifdef __LINBOX_HAVE_NTL
									   UnparametricField<NTL::ZZ_pE>* p = new UnparametricField<NTL::ZZ_pE>(_R);
									   return FT(p, AuxUserData);
#endif // if we don't have the package, fall through to the
									   // "we'll use anything!" stage
								   }
#if 0
								   else if(_implDetail == "givaro-gfq") {
#ifdef __LINBOX_HAVE_GIVARO
									   GivaroGfq *p = new GivaroGfq(_R);
									   return FT(p, AuxUserData);
#endif
								   }
								   else if(_implDetail == "lidia-gfq") {
#ifdef __LINBOX_HAVE_LIDIA
									   LidiaGfq *p = new LidiaGfq(_R);
									   return FT(p, AuxUserData);
#endif
								   }
#endif
								   // Now we are in desperation mode, just use the first
								   // plug-in package that we've got
								   //

								   else {
#if 0
#ifdef __LINBOX_HAVE_GIVARO
									   GivaroGfq *p = new GivaroGfq(_R);
									   return FT(p, AuxUserData);
#elif __LINBOX_HAVE_LIDIA
									   LidiaGfq *p = new LidiaGfq(_R);
									   return FT(p, AuxUserData);
#endif
#endif

#ifdef __LINBOX_HAVE_NTL // Our last resort
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
#ifdef __LINBOX_HAVE_NTL
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
	void* FieldReaderAnalyzer::makeField(FunctorType &FT) const
	{

		switch(whatType()) {

			// For this type, there are the follwing Native
			// LinBox types:
			//   - GF2
			//   - Modular<type> where type is one of
			//          - uint8_t
			//          - uint16_t
			//          - uint32_t
			//          - integer (This is the default for both
			//                     Modular & for all Finite fields in LB)
			//
			// There are also the following types that come from plug-in
			// packages:
			//   - From NTL:
			//      - NTL::ZZ_p (defaults to modular<integer>)
			//      - NTL::zz_p (this one defaults to modular<uint32_t>)
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
								  Modular<uint8_t>* p = new Modular<uint8_t>(_R);
								  return FT(p);
							  }
							  else if(_implDetail == "modular-uint16") {
								  Modular<uint16_t>* p = new Modular<uint16_t>(_R);
								  return FT(p);
							  }
							  else if(_implDetail == "modular-uint32_t") {
								  Modular<uint32_t>* p = new Modular<uint32_t>(_R);
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
#ifdef __LINBOX_HAVE_NTL // if we have NTL
								  UnparametricField<NTL::zz_p> *p = new UnparametricField<NTL::zz_p>(_R);
#else // the default
								  Modular<uint32_t> *p = new Modular<uint32_t>(_R);
#endif
								  return FT(p);

							  }
							  else if(_implDetail == "ntl-ZZp") {
#ifdef __LINBOX_HAVE_NTL // if we have NTL
								  UnparametricField<NTL::ZZ_p> *p = new UnparametricField<NTL::ZZ_p>(_R);
#else
								  Modular<integer> *p = new Modular<integer>(_R);
#endif
								  return FT(p);
							  }
#if 0
							  else if(_implDetail == "givaro-zpz-std16") {
#ifdef __LINBOX_HAVE_GIVARO // if we have Givaro
								  GivaroZpz<Std16> *p = new GivaroZpz<Std16>(_R);
#else
								  Modular<uint16_t> *p = new Modular<uint16_t>(_R);
#endif
								  return FT(p);
							  }
							  else if(_implDetail == "givaro-zpz-std32") {
#ifdef __LINBOX_HAVE_GIVARO
								  GivaroZpz<Std32> *p = new GivaroZpz<Std32>(_R);
#else
								  Modular<uint32_t> *p = new Modular<uint32_t>(_R);
#endif
								  return FT(p);
							  }
							  else if(_implDetail == "givaro-zpz-log16") {
#ifdef __LINBOX_HAVE_GIVARO
								  GivaroZpz<Log16> *p = new GivaroZpz<Log16>(_R);
#else

								  // crap, we're screwed, just throw the
								  // exception and bail
								  throw NoTableFiniteInLinBox;
#endif
								  return FT(p);
							  }
#endif
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
#ifdef __LINBOX_HAVE_NTL
									   UnparametricField<NTL::ZZ_pE>* p = new UnparametricField<NTL::ZZ_pE>(_R);
									   return FT(p);
#endif // if we don't have the package, fall through to the
									   // "we'll use anything!" stage
								   }
#if 0
								   else if(_implDetail == "givaro-gfq") {
#ifdef __LINBOX_HAVE_GIVARO
									   GivaroGfq *p = new GivaroGfq(_R);
									   return FT(p);
#endif
								   }
								   else if(_implDetail == "lidia-gfq") {
#ifdef __LINBOX_HAVE_LIDIA
									   LidiaGfq *p = new LidiaGfq(_R);
									   return FT(p);
#endif
								   }
#endif
								   // Now we are in desperation mode, just use the first
								   // plug-in package that we've got
								   //
								   else {
#if 0
#ifdef __LINBOX_HAVE_GIVARO
									   GivaroGfq *p = new GivaroGfq(_R);
									   return FT(p);
#elif __LINBOX_HAVE_LIDIA
									   LidiaGfq *p = new LidiaGfq(_R);
									   return FT(p);
#endif
#endif
#ifdef __LINBOX_HAVE_NTL // Our last resort
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
#ifdef __LINBOX_HAVE_NTL
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

} // namespace LinBox

#endif //__LINBOX_field_reader_analyzer_H

