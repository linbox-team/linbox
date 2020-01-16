/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include "common_sdmb.hpp"
#include "common_sdmp.hpp"

// #define __USE_BINARY_SDM

#if defined(__USE_BINARY_SDM)
    #define SDM_EXTENSION "sdmb"
    #define SDM_READ_HEADER(...) readSDMBHeader(__VA_ARGS__)
    #define SDM_READ_MATRIX(...) readSDMBMatrix(__VA_ARGS__)
    #define SDM_WRITE_HEADER(...) writeSDMBHeader(__VA_ARGS__)
    #define SDM_WRITE_MATRIX(...) writeSDMBMatrix(__VA_ARGS__)
#else
    #define SDM_EXTENSION "sdmp"
    #define SDM_READ_HEADER(...) readSDMPHeader(__VA_ARGS__)
    #define SDM_READ_MATRIX(...) readSDMPMatrix(__VA_ARGS__)
    #define SDM_READ_VECTOR(...) readSDMPVector(__VA_ARGS__)
    #define SDM_WRITE_HEADER(...) writeSDMPHeader(__VA_ARGS__)
    #define SDM_WRITE_MATRIX(...) writeSDMPMatrix(__VA_ARGS__)
#endif

