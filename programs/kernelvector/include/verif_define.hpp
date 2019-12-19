/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include <recint/ruint.h>
#include <givaro/montgomery-ruint.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/utils/flimits.h>

#include "common_common.hpp"

#include <set>
#include <list>
#include <vector>

//------------------//
//----- Timers -----//

#include <fflas-ffpack/utils/timer.h>

using TTimer = FFLAS::OMPTimer;

//--------------------//
//----- Typedefs -----//

#ifdef __USE_128bits
using FieldElement_t = RecInt::ruint128;
using FieldModulo_t = RecInt::ruint128;
#else
using FieldElement_t = RecInt::ruint64;
using FieldModulo_t = RecInt::ruint64;
#endif

using Field_t   = Givaro::Montgomery<FieldElement_t>;
using Proj_t    = std::set<Index_t>;
using Data_t    = Field_t::Element_ptr;
using DMatrix_t = std::vector<FieldElement_t>;
using SMatrix_t = FFLAS::Sparse<Field_t, FFLAS::SparseMatrix_t::HYB_ZO>;


