
#ifndef _GIVARO_POLY_MOD_POLY_H
#define _GIVARO_POLY_MOD_POLY_H

#include <givaro/extension.h>
#include <givaro/givindeter.h>
#include "linbox/field/Givaro/givaro-field.h"

namespace LinBox {

template <class BaseField>
class GivaroPolyModPoly {
public:
	typedef Givaro::Extension<GivaroField<BaseField> > Parent_t;
	typedef typename Parent_t::PolElement Element;

	typedef GivaroField<BaseField> Domain_t;
	typedef typename GivaroField<BaseField>::Element Type_t;

	GivaroPolyModPoly(BaseField& F,Element p,Givaro::Indeter Y="Y") :
		F_(&F),
		FactorDom_(*F_,Y),
		ExtensionField_(FactorDom_,GivaroField<BaseField>(p)) {}

	Parent_t* getExtension() {
		return &ExtensionField_;
	}

private:

	GivaroField<BaseField> *F_;

	Givaro::Poly1FactorDom<GivaroField<BaseField>,Givaro::Dense> FactorDom_;

	Givaro::Extension<GivaroField<BaseField> > ExtensionField_;

};

}

#endif // _GIVARO_POLY_MOD_POLY_H
