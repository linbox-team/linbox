/* tests/test-givaro-interfaces.C
 * Copyright (C) 2019 The LinBox group
 *
 * Written by Jean-Guillaume Dumas
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file  tests/test-givaro-interfaces.C
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
    static integer q = 101;
    
    static Argument args[] = {
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
        END_OF_ARGUMENTS
    };
    
    parseArguments (argc, argv, args);
    

    bool pass(true);
    
    pass &= test_archetype<Modular< uint16_t>>(q);
    pass &= test_archetype<Modular< uint32_t>>(q);
    pass &= test_archetype<Modular< Givaro::Log16>>(q);
    pass &= test_archetype<GFqDom<int64_t>>(q);
    pass &= test_archetype<Modular< int16_t>> (q);
    pass &= test_archetype<Modular< int32_t>>(q);
    pass &= test_archetype<ModularBalanced<float>>(q);
    pass &= test_archetype<ModularBalanced<double>>(q);
    pass &= test_archetype<ModularBalanced<int32_t>>(q);
    pass &= test_archetype<ModularBalanced<int64_t>>(q);
    pass &= test_archetype<ModularExtended<float>>(q);
    pass &= test_archetype<ModularExtended<double>>(q);
    pass &= test_mongt<int32_t>(q);
    pass &= test_mongt<RecInt::ruint<7>>(q);

    return pass ? 0 : -1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
