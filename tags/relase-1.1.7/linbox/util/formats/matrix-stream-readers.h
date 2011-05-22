/* Copyright (C) LinBox
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

#ifndef __LINBOX_matrix_stream_readers_H
#define __LINBOX_matrix_stream_readers_H

/*! \file matrix-stream-readers.h
 * Here is where all the formats (each of which is a subclass of
 * MatrixStreamReader) are defined, in two places:
 *
 * First, the macro __MATRIX_STREAM_READERDEFS is put in the init() function of
 * MatrixStream and it creates all the format readers.  For each format reader,
 * a line of this macro should read: addReader( new MyReaderType() );
 *
 * Second, so those statements actually compile, the file containing each format
 * reader should be included with a line of the form: #include "my-reader.h"
 */

#define __MATRIX_STREAM_READERDEFS \
	addReader( new SMSReader<Field>() ); \
	addReader( new SparseRowReader<Field>() ); \
	addReader( new MatrixMarketReader<Field>() ); \
	addReader( new MapleReader<Field>() ); \
	addReader( new DenseReader<Field>() );

#include "sms.h"
#include "sparse-row.h"
#include "generic-dense.h"
#include "matrix-market.h"
#include "maple.h"

#endif //__LINBOX_matrix_stream_readers_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
