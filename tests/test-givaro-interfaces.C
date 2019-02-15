/* tests/test-givaro-zpz.C
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file  tests/test-givaro-zpz.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>


#include "givaro/modular.h"
#include "givaro/modular-extended.h"
#include "givaro/montgomery.h"
#include "givaro/gfq.h"

#include "test-common.h"
#include "test-field.h"

using namespace Givaro;
using namespace LinBox;
#include <typeinfo>       // operator typeid

template<typename Inttype, typename FField>
bool test_type_archetype(const Inttype& q) {
	bool pass(true);
    FField * K = new FField(q);
	FieldArchetype FA(K);
    std::string str("Testing archetype with envelope of ");
    str += typeid(*K).name();
    str += " field";    
	if (!testField<FieldArchetype> (FA, str.c_str()) )
		pass = false;
	delete K;
    return pass;
}

template<typename FField>
bool test_archetype(const integer& q) { 
    return test_type_archetype<integer,FField>(q);
}

template<typename Inttype>
bool test_mongt(const integer& q) {
    return test_type_archetype<Inttype,Montgomery<Inttype>>( (Inttype)q );
}



int main (int argc, char **argv)
{
    static integer q = 10733;
	static size_t n = 10000;
	static int iterations = 10;
	static int e = 3;
	static int trials = 10000;
	static int categories = 1000;
	static int hist_level = 10;

    
    static Argument args[] = {
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
        { 'e', "-e E", "Use GF(q^e) for the extension field [1].",  TYPE_INT,     &e },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
    };
    
    parseArguments (argc, argv, args);
    

    bool pass(true);
    
    pass &= test_archetype<Modular< uint16_t>>(101);
    pass &= test_archetype<Modular< uint32_t>>(101);
	pass &= test_archetype<Modular< Givaro::Log16>>(101);
	pass &= test_archetype<GFqDom<int64_t>>(101);
    pass &= test_archetype<Modular< int16_t>> (101);
	pass &= test_archetype<Modular< int32_t>>(101);
    pass &= test_archetype<ModularBalanced<float>>(101);
    pass &= test_archetype<ModularBalanced<double>>(101);
    pass &= test_archetype<ModularBalanced<int32_t>>(101);
    pass &= test_archetype<ModularBalanced<int64_t>>(101);
    pass &= test_archetype<ModularExtended<float>>(101);
    pass &= test_archetype<ModularExtended<double>>(101);
    pass &= test_mongt<int32_t>(101);
    pass &= test_mongt<RecInt::ruint<7>>(101);

    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
