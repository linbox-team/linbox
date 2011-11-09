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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_reader_blackbox_factory_H
#define __LINBOX_reader_blackbox_factory_H

#include "linbox/linbox-config.h"

#include "linbox-reader.h"
#include "linbox/blackbox/archetype.h"

#include <string>

namespace LinBox
{

	class UnSupportedMatrixType {};


	template<class Vector>
	class ReaderBlackBoxFactory {
	public:
		typedef BlackboxArchetype<Vector> BlackBox;

		ReaderBlackBoxFactory();
		ReaderBlackBoxFactory(Reader &R);

		bool reset(Reader &R);


		static const int Dense = 0;
		static const int Sparse = 1;
		static const int Special = 2;
		static const int NotBlackBox = 3;

		bool isBlackBox() const;
		bool hasField() const;
		int whatType() const;

		std::string implDetail() const;

		BlackBox* makeBlackBox();

		template<class Field>
		BlackBox* makeBlackBox(Field *F);

		template<class Field>
		void *operator()(Field *F) {
			// we don't need the field, just the type
			delete F;
			return (void*) makeBlackBox<Field>(F);
		}


	private:
		bool _isBlackbox, _hasField;

		int _type;
		std::string _implDetail;
		Reader _R;
	};

}


#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/nag-sparse.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/dif.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/transpose.h"

#ifdef __LINBOX_HAVE_NTL
#include "linbox/blackbox/ntl-toeplitz.h"
#endif

#include "field-reader-analyzer.h"

namespace LinBox
{


	template<class Vector>
	ReaderBlackBoxFactory<Vector>::ReaderBlackBoxFactory() {

		_isBlackbox = false;
		_hasField = false;
		_type = NotBlackBox;
	}

	template<class Vector>
	ReaderBlackBoxFactory<Vector>::ReaderBlackBoxFactory(Reader &R) {

		reset(R);
	}

	template<class Vector>
	bool ReaderBlackBoxFactory<Vector>::reset(Reader &R) {
		std::string s;


		_R = R;
		if(!_R.expectTagName("MatrixOver") | !_R.expectChildTag()) {
			_isBlackbox = false;
			_hasField = false;
			_type = ReaderBlackBoxFactory::NotBlackBox;
			_implDetail = s;
			return false;
		}
		if(!_R.checkAttributeString("implDetail", _implDetail))
			_implDetail = s;

		_R.traverseChild();
		if(_R.checkTagName("field")) {
			_R.upToParent();
			_hasField = true;
			_R.getNextChild();
			if(!_R.expectChildTag()) {
				_isBlackbox = false;
				_hasField = false;
				_type = ReaderBlackBoxFactory::NotBlackBox;
				_implDetail = s;
				return false;
			}
			_R.traverseChild();
		}
		else {
			_hasField = false;
		}

		if(_R.checkTagName("matrix"))
			_type = ReaderBlackBoxFactory::Dense;
		else if(_R.checkTagName("sparseMatrix"))
			_type = ReaderBlackBoxFactory::Sparse;
		else
			_type = ReaderBlackBoxFactory::Special;

		_R.Up(1).Left(1);
		_isBlackbox = true;
		return true;
	}

	template<class Vector>
	bool ReaderBlackBoxFactory<Vector>::isBlackBox() const {
		return _isBlackbox;
	}

	template<class Vector>
	bool ReaderBlackBoxFactory<Vector>::hasField() const {
		return _hasField;
	}

	template<class Vector>
	int ReaderBlackBoxFactory<Vector>::whatType() const {
		return _type;
	}

	template<class Vector>
	std::string ReaderBlackBoxFactory<Vector>::implDetail() const {
		return _implDetail;
	}


	// This is the first method.  It returns a valid pointer if
	// the operation was successful, and NULL otherwise.
	// This is the most general method.  It first takes in a Reader and
	// checks that the Reader describes a Blackbox (with the MatrixOver
	// tag.  If it has a field, a FieldReaderAnalyzer is used to attempt
	// to read the Field and generate the correct type.  If there is no
	// field, this function builds the Blackbox matrix

	template<class Vector>
	BlackboxArchetype<Vector>* ReaderBlackBoxFactory<Vector>::makeBlackBox()
	{
		FieldReaderAnalyzer FA;

		// If there is an error, return NULL
		if(!_isBlackbox) {
			return NULL;
		}

		if(_hasField) {

			FA.reset(_R.Down(1));
			_R.Up(1);

			// calls the overloaded operator() method of this
			// class for the proper use of the field analyzer
			return (BlackBox*) FA.makeField(*this);
		}

		else {

			// go on the major types
			_R.Down(1);
			if(_R.checkTagName("permutation")) {
				Permutation<Vector>* p = new Permutation<Vector>(_R.Up(1));
				return p;
			}

			else if(_R.checkTagName("compose")) {
				Compose<Vector>* p =
				new Compose<Vector>(_R.Up(1));
				return p;
			}
			/*
			   else if(_R.checkTagName("transpose")) {
			   Transpose<Vector>* p = new Transpose<Vector>(_R.Up(1));
			   return p;
			   }
			   */
			else {
				_R.Up(1);
				throw UnSupportedMatrixType();
			}
		}

		return NULL; // we of course don't get here, but to be safe
	}

	// don't play w/ the pointer parameter, it has been deleted and points
	// to garbage
	//
	template<class Vector>
	template<class Field>
	BlackboxArchetype<Vector>*  ReaderBlackBoxFactory<Vector>::makeBlackBox(Field *F)
	{
		typedef typename Field::Element Element;

		if(!_isBlackbox) {
			return NULL;
		}

		switch(_type) {


			// for this type, we have the following blackboxes:
			// SparseMatrix (templatized on various types)
			// and TriplesBB.  The default is
			// SparseMatrix templatized on dense vectors

		case ReaderBlackBoxFactory<Vector>::Sparse : {


								     if(_implDetail == "sparse-sequence") {
									     SparseMatrix<Field, Vector, vector<pair<size_t, Element> > >* p = new SparseMatrix<Field, Vector, vector<pair<size_t, Element> > >(_R);

									     return p;
								     }
								     else if(_implDetail == "sparse-associative") {
									     SparseMatrix<Field, Vector, map<size_t, Element> >* p = new SparseMatrix<Field, Vector, map<size_t, Element> >(_R);
									     return p;
								     }
								     else if(_implDetail == "sparse-parallel") {
									     SparseMatrix<Field, Vector, pair<vector<size_t>, vector<Element> > >* p = new SparseMatrix<Field, Vector, pair<vector<size_t>, vector<Element> > >(_R);
									     return p;
								     }
								     else if(_implDetail == "triplesbb") {
									     TriplesBB<Field, Vector>* p = new TriplesBB<Field, Vector>(_R);
									     return p;
								     }
								     // This is our fall-back case
								     else {
									     SparseMatrix<Field, Vector, vector<pair<size_t, Element> > >* p = new SparseMatrix<Field, Vector, vector<pair<size_t, Element> > >(_R);
									     return p;
								     }
							     }
							     break;

		case ReaderBlackBoxFactory<Vector>::Dense : {
								    DenseMatrix<Field, Vector>* p = new DenseMatrix<Field, Vector>(_R);
								    return p;
							    }
							    break;

		case ReaderBlackBoxFactory<Vector>::Special : {
								      if(_implDetail == "scalar") {
									      ScalarMatrix<Field, Vector>* p = new ScalarMatrix<Field, Vector>(_R);
									      return p;
								      }
								      else if(_implDetail == "zero-one") {
									      ZeroOne<Field, Vector>* p = new ZeroOne<Field, Vector>(_R);
									      return p;
								      }
								      else if(_implDetail == "diagonal-dense") {
									      Diagonal<Field, Vector>* p = new Diagonal<Field, Vector>(_R);
									      return p;
								      }
								      else if(_implDetail == "diagonal-sequence") {
									      Diagonal<Field, Vector>* p = new Diagonal<Field, Vector>(_R);
									      return p;
								      }
								      else if(_implDetail == "diagonal-associative") {
									      Diagonal<Field, Vector>* p = new Diagonal<Field, Vector>(_R);
									      return p;
								      }
								      else if(_implDetail == "ntl-toeplitz") {
									      Toeplitz<Field, Vector>* p = new Toeplitz<Field, Vector>(_R);
									      return p;
								      }
								      /*
									 else if(_implDetail == "sum") {
									 Sum<Field, Vector>* p = new Sum<Field, Vector>(_R);
									 return p;
									 }
									 else if(_implDetail == "dif") {
									 Dif<Field, Vector>* p = new Dif<Field, Vector>(_R);
									 return p;
									 }
									 else if(_implDetail == "submatrix") {
									 Submatrix<Field, Vector>* p = new Submatrix<Field, Vector>(_R);
									 return p;
									 }
									 else if(_implDetail == "inverse") {
									 Inverse<Field, Vector>* p = new Inverse<Field, Vector>(_R);
									 return p;
									 }
									 else if(_implDetail == "moore-penrose") {
									 MoorePenrose<Field, Vector>* p = new MoorePenrose<Field, Vector>(_R);
									 return p;
									 }
									 */
								      else {
									      throw UnSupportedMatrixType();
								      }
							      }
							      break;

							      // we shouldn't ever get here
		case ReaderBlackBoxFactory::NotBlackBox :
		default:

							      return NULL;
							      break;
		}

		return NULL; // or here
	}

}


#endif //__LINBOX_reader_blackbox_factory_H

