#ifndef __READER_BLACKBOX_FACTORY
#define __READER_BLACKBOX_FACTORY

#include "linbox-reader.h"

using LinBox::Reader;

namespace LinBox {

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

		string implDetail() const;


		BlackBox* makeBlackBox(Reader &) const;

		template<class Field>
		void* makeBlackBox(Reader &) const;
		
		template<class Field>
		void *operator()(Field *F, Reader &R) const {
			// we don't need the field, just the type
			delete F;
			return (void*) makeBlackBox<Field>(R);
		}


	private:
		bool _isBlackbox, _hasField;

		int _type;
		string _implDetail;
		Reader _R;
	};


				

		
		
