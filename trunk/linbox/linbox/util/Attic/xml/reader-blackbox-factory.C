#include "reader-blackbox-factory.h"

#include "linbox-config.h"
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

#ifdef HAVE_NTL
#include "linbox/blackbox/ntl-toeplitz.h"
#endif

#include "field-reader-analyzer.h"

namespace LinBox {


	template<class Vector>
	ReaderBlackBoxFactory<Vector>::ReaderBlackBoxFactory() {

		_isBlackbox = false;
		_hasField = false;
		_type = NotBlackBox;
	}

	template<class Vector>
	ReaderBlacBoxFactory<Vector>::ReaderBlackboxFactory(Reader &R) {

		reset(R);
	}

	template<class Vector>
	bool ReaderBlackBoxFactory<Vector>::reset(Reader &R) {
		string s;


		_R = R;
		if(!_R.expectTagName("MatrixOver") | !_R.expectChildTag()) {
			_isBlackBox = false;
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
				_isBlackBox = false;
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
		_isBlackBox = true;
		return true;
	}

	template<class Vector>
	bool ReaderBlackBoxFactory<Vector>::isBlackBox() const {
		return _isBlackBox;
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
	string ReaderBlackBoxFactory<Vector>::implDetail() const {
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
	BlackboxArchetype<Vector>* ReaderBlackBoxFactory<Vector>::makeBlackBox(Reader &R) const 
	{
		FieldAnalyzer FA;

		if(_haveField) {

			FA.reset(_R.Down(1));
			_R.Up(1);

			// calls the overloaded operator() method of this
			// class for the proper use of the field analyzer
			return (BlackBox*) FA.makeField(*this, _R);
		}

		else {

			// go on the major types
			_R.Down(1);
			if(_R.checkTagName("permutation")) {
				Permutation<Vector>* p = new Permutation<Vector>(_R.Up(1));
				return p;
			}
			else if(_R.checkTagName("compose")) {
				Compose<Vector>* p = new Compose<Vector>(_R.Up(1));
				return p;
			}
			else if(_R.checkTagName("transpose")) {
				Transpose<Vector>* p = new Transpose<Vector>(_R.Up(1));
				return p;
			}
			else {
				throw UnsupportedMatrixType;
			}
		}

		return NULL; // we of course don't get here, but to be safe
	}



}
