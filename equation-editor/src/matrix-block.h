/* -*- mode: c; style: linux -*- */

/* matrix-block.h
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifndef __MATRIX_BLOCK_H
#define __MATRIX_BLOCK_H

#include <gnome.h>

#include "block.h"

BEGIN_GNOME_DECLS

#define MATRIX_BLOCK(obj)          GTK_CHECK_CAST (obj, matrix_block_get_type (), MatrixBlock)
#define MATRIX_BLOCK_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, matrix_block_get_type (), MatrixBlockClass)
#define IS_MATRIX_BLOCK(obj)       GTK_CHECK_TYPE (obj, matrix_block_get_type ())

typedef struct _MatrixBlock MatrixBlock;
typedef struct _MatrixBlockClass MatrixBlockClass;
typedef struct _MatrixBlockPrivate MatrixBlockPrivate;

struct _MatrixBlock 
{
	Block parent;

	MatrixBlockPrivate *p;
};

struct _MatrixBlockClass 
{
	BlockClass block_class;
};

guint matrix_block_get_type         (void);

GtkObject *matrix_block_new         (void);

END_GNOME_DECLS

#endif /* __MATRIX_BLOCK_H */
