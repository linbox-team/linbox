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

/**
 * matrix_block_new:
 * @rows: Initial number of rows
 * @cols: Initial number of columns
 * 
 * Construct a new MatrixBlock with the given number of rows and columns
 *
 * Return value: The newly constructed matrix block
 **/

GtkObject *matrix_block_new         (guint rows, guint cols);

/**
 * matrix_block_insert_row:
 * @block: object
 * @position: Row number before which to insert new row, -1 to append
 * 
 * Insert a row into the given position; if position is negative or greater
 * than the number of rows, append to the end of the matrix
 **/

void matrix_block_insert_row        (MatrixBlock *block,
				     gint position);

/**
 * matrix_block_insert_col:
 * @block: object
 * @position: Column number before which to insert new column, -1 to append
 * 
 * Insert a column into the matrix at the given position. If position is
 * negative or greater than the number of columns, append the column to the
 * end.
 **/

void matrix_block_insert_col        (MatrixBlock *block,
				     gint position);

/**
 * matrix_block_remove_row:
 * @block: object
 * @position: Row number to remove, -1 to remove the last row
 * 
 * Remove the given row from the matrix, unrefing all objects therein. If
 * position is negative, remove the last row
 **/

void matrix_block_remove_row        (MatrixBlock *block,
				     gint position);

/**
 * matrix_block_remove_col:
 * @block: object
 * @position: Column number to remove, or -1 to remove the last column
 * 
 * Remove the given column from the matrix, unrefing all objects therein. If
 * position is negative, remove the last column
 **/

void matrix_block_remove_col        (MatrixBlock *block,
				     gint position);

/**
 * matrix_block_set_math_object:
 * @block: object
 * @row: Row at which to place object (0 <= row < no. rows)
 * @col: Column at which to place object (0 <= col < no. cols)
 * @math_object: The new math object to place at that coordinate
 * 
 * Set the object at (row, col) to the given math object. Unref the existing
 * math object, if it exists, and ref the new object
 **/

void matrix_block_set_math_object   (MatrixBlock *block,
				     guint row, guint col,
				     MathObject *math_object);

/**
 * matrix_block_get_math_object:
 * @block: object
 * @row: The row from which to fetch the object (0 < row < no. rows)
 * @col: The column from which to fetch the object (0 < col < no. cols)
 * 
 * Fetch the math object at (row, col)
 * 
 * Return value: The math object at (row, col)
 **/

MathObject *matrix_block_get_math_object (MatrixBlock *block,
					  guint row, guint col);

END_GNOME_DECLS

#endif /* __MATRIX_BLOCK_H */
