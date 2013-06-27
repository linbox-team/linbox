/*
 * examples/fields/ex-fields.C
 *
 * Copyright (C) 2001, 2002, 2010 G. Villard
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/** \file examples/fields/ex-fields.C
 * \author Gilles Villard
 * \brief Using a template function with two distinct fields.
 */

// ---------------------------------------------

#include "linbox-config.h"

#include <iostream>

// ---------------------------------------------
#include "linbox/field/modular.h"
#include "linbox/field/ntl.h"

using namespace LinBox;
using namespace std;

/**  The template function "fct" reads two elements "a" and "b" of the
 *  field "K" from the standard input and writes a/b on the standard output */

template <class Field>
void divide_ex(const Field&  K)
{

	/* "K" is a field domain (a C++ object) of type "Field" (here the template
	 *  parameter). The type of the elements of "K" is obtained through
	 * "Field" by typename Field::Element */

	typedef typename Field::Element K_elt;

	K_elt a,b,r;

	K.init(a); K.init(b); K.init(r);
	cout << "division example: enter two numbers: ";
	K.read(cin,a);  K.read(cin,b);
	K.div(r,a,b);
	K.write( K.write(cout<< "the quotient in ") << " is ",
		 r) << endl;

}

// ---------------------------------------------
/// no command line args
int main()
{

	/* Using the parameterized domain capabilities, several domains
	 * representing integers modulo may be used simultaneously. */

	Modular<uint32_t> D(3), K(7);

	divide_ex(D);  divide_ex(K);

	// NTL arbitrary precision real field
	// (Could be parameterized by the precision)

	// UnparametricField<NTL::RR> K2;
	NTL_RR K2 ;
	NTL::RR::SetPrecision(500);
	NTL::RR::SetOutputPrecision(50);

	// NTL modulo p field

	//UnparametricField<NTL::zz_p> K2;
	//NTL::zz_p::init(553);

	divide_ex(K2);

	return 0;
};


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

