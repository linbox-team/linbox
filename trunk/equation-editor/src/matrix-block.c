/* -*- mode: c; style: linux -*- */

/* matrix-block.c
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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "matrix-block.h"
#include "matrix-block-layout.h"

enum {
	ARG_0,
	ARG_ROWS,
	ARG_COLS
};

struct _MatrixBlockPrivate 
{
	guint rows;
	guint cols;
	MathObject ***objects;
};

static BlockClass *parent_class;

static MatrixBlockLayout *layout;

static void matrix_block_init        (MatrixBlock *matrix_block);
static void matrix_block_class_init  (MatrixBlockClass *class);

static void matrix_block_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void matrix_block_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void matrix_block_finalize    (GtkObject *object);

static const Layout *matrix_block_get_layout (MathObject *math_object);

static void matrix_block_foreach     (Block *block,
				      BlockIteratorCB callback,
				      gpointer data);

static void destroy_object_array     (MatrixBlock *matrix_block);
static void setup_object_array       (MatrixBlock *matrix_block);

guint
matrix_block_get_type (void)
{
	static guint matrix_block_type = 0;

	if (!matrix_block_type) {
		GtkTypeInfo matrix_block_info = {
			"MatrixBlock",
			sizeof (MatrixBlock),
			sizeof (MatrixBlockClass),
			(GtkClassInitFunc) matrix_block_class_init,
			(GtkObjectInitFunc) matrix_block_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		matrix_block_type = 
			gtk_type_unique (block_get_type (), 
					 &matrix_block_info);
	}

	return matrix_block_type;
}

static void
matrix_block_init (MatrixBlock *matrix_block)
{
	matrix_block->p = g_new0 (MatrixBlockPrivate, 1);
}

static void
matrix_block_class_init (MatrixBlockClass *class) 
{
	GtkObjectClass *object_class;
	MathObjectClass *math_object_class;
	BlockClass *block_class;

	gtk_object_add_arg_type ("MatrixBlock::rows",
				 GTK_TYPE_INT,
				 GTK_ARG_READWRITE,
				 ARG_ROWS);

	gtk_object_add_arg_type ("MatrixBlock::cols",
				 GTK_TYPE_INT,
				 GTK_ARG_READWRITE,
				 ARG_COLS);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = matrix_block_finalize;
	object_class->set_arg = matrix_block_set_arg;
	object_class->get_arg = matrix_block_get_arg;

	parent_class = BLOCK_CLASS
		(gtk_type_class (block_get_type ()));

	math_object_class = MATH_OBJECT_CLASS (class);
	math_object_class->get_layout = matrix_block_get_layout;

	block_class = BLOCK_CLASS (class);
	block_class->foreach = matrix_block_foreach;

	if (layout == NULL)
		layout = ROW_BLOCK_LAYOUT (matrix_block_layout_new ());
}

static void
matrix_block_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MatrixBlock *matrix_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (object));

	matrix_block = MATRIX_BLOCK (object);

	switch (arg_id) {
	case ARG_ROWS:
		if (matrix_block->p->objects != NULL)
			destroy_object_array (matrix_block);

		matrix_block->p->rows = GTK_VALUE_INT (*arg);

		if (matrix_block->p->cols > 0)
			setup_object_array (matrix_block);

		break;

	case ARG_COLS:
		if (matrix_block->p->objects != NULL)
			destroy_object_array (matrix_block);

		matrix_block->p->cols = GTK_VALUE_INT (*arg);

		if (matrix_block->p->rows > 0)
			setup_object_array (matrix_block);

		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
matrix_block_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MatrixBlock *matrix_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (object));

	matrix_block = MATRIX_BLOCK (object);

	switch (arg_id) {
	case ARG_ROWS:
		GTK_VALUE_INT (*arg) = matrix_block->p->rows;
		break;

	case ARG_COLS:
		GTK_VALUE_INT (*arg) = matrix_block->p->cols;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
matrix_block_finalize (GtkObject *object) 
{
	MatrixBlock *matrix_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (object));

	matrix_block = MATRIX_BLOCK (object);

	destroy_object_array (matrix_block);
	g_free (matrix_block->p);
}

GtkObject *
matrix_block_new (guint rows, guint cols) 
{
	return gtk_object_new (matrix_block_get_type (),
			       "rows", rows,
			       "cols", cols,
			       NULL);
}

/**
 * matrix_block_insert_row:
 * @block: 
 * @position: 
 * 
 * Insert a row into the given position; if position is negative or greater
 * than the number of rows, append to the end of the matrix
 **/

void
matrix_block_insert_row (MatrixBlock *block, gint position)
{
	gint i;

	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (block));

	if (position > block->p->rows) position = -1;

	block->p->rows++;
	block->p->objects = g_renew (int *, block->p->objects, block->p->rows);

	if (position >= 0) {
		for (i = block->p->rows; i > position; i--)
			block->p->objects[i] = block->p->objects[i - 1];
		block->p->objects[position] = 
			g_new0 (MathObject *, block->p->cols);
	} else {
		block->p->objects[block->p->rows - 1] =
			g_new0 (MathObject *, block->p->cols);
	}
}

/**
 * matrix_block_insert_col:
 * @block: 
 * @position: 
 * 
 * Insert a column into the matrix at the given position. If position is
 * negative or greater than the number of columns, append the column to the
 * end.
 **/

void
matrix_block_insert_col (MatrixBlock *block, gint position)
{
	gint i, j;

	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (block));

	if (position > block->p->cols) position = -1;

	block->p->cols++;

	for (i = 0; i < block->p->rows; i++) {
		block->p->objects[i] = 
			g_renew (MathObject *, block->p->objects[i],
				 block->p->cols);

		if (position >= 0) {
			for (j = block->p->cols; j > position; j--)
				block->p->objects[i][j] =
					block->p->objects[i][j - 1];
			block->p->objects[i][position] = NULL;
		} else {
			block->p->objects[i][block->p->cols - 1] = NULL;
		}
	}
}

/**
 * matrix_block_remove_row:
 * @block: 
 * @position: 
 * 
 * Remove the given row from the matrix, unrefing all objects therein. If
 * position is negative, remove the last row
 **/

void
matrix_block_remove_row (MatrixBlock *block, gint position)
{
	gint i;

	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (block));
	g_return_if_fail (position > block->p->rows);

	block->p->rows--;

	if (position >= 0)
		for (i = position; i < block->p->rows; i++)
			block->p->objects[i] = block->p->objects[i + 1];

	block->p->objects = g_renew (int *, block->p->objects, block->p->rows);
}

/**
 * matrix_block_remove_col:
 * @block: 
 * @position: 
 * 
 * Remove the given column from the matrix, unrefing all objects therein. If
 * position is negative, remove the last column
 **/

void
matrix_block_remove_col (MatrixBlock *block, gint position)
{
	gint i, j;

	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (block));
	g_return_if_fail (position > block->p->cols);

	block->p->cols--;

	for (i = 0; i < block->p->rows; i++) {
		if (position >= 0) {
			for (j = poition; j < block->p->cols; j++)
				block->p->objects[i][j] =
					block->p->objects[i][j + 1];
		}

		block->p->objects[i] = 
			g_renew (MathObject *, block->p->objects[i],
				 block->p->cols);
	}
}

/**
 * matrix_block_set_math_object:
 * @block: 
 * @row: Row at which to place object (0 <= row < no. rows)
 * @col: Column at which to place object (0 <= col < no. cols)
 * @math_object: 
 * 
 * Set the object at (row, col) to the given math object. Unref the existing
 * math object, if it exists, and ref the new object
 **/

void
matrix_block_set_math_object (MatrixBlock *block,
			      guint row, guint col,
			      MathObject *math_object)
{
	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (block));
	g_return_if_fail (row > block->p->rows);
	g_return_if_fail (col > block->p->cols);

	if (block->p->objects[row][col] != NULL)
		gtk_object_unref (block->p->objects[row][col]);

	block->p->objects[row][col] = math_object;

	if (math_object != NULL)
		gtk_object_ref (math_object);
}

/**
 * matrix_block_get_math_object:
 * @block: 
 * @row: The row from which to fetch the object (0 < row < no. rows)
 * @col: The column from which to fetch the object (0 < col < no. cols)
 * 
 * Fetch the math object at (row, col)
 * 
 * Return value: 
 **/

MathObject *
matrix_block_get_math_object (MatrixBlock *block, guint row, guint col)
{
	g_return_val_if_fail (block != NULL, NULL);
	g_return_val_if_fail (IS_MATRIX_BLOCK (block), NULL);
	g_return_val_if_fail (row > block->p->rows, NULL);
	g_return_val_if_fail (col > block->p->cols, NULL);

	return block->objects[row][col];
}

static const Layout *
matrix_block_get_layout (MathObject *math_object) 
{
	return LAYOUT (layout);
}

static void
matrix_block_foreach (Block *block, BlockIteratorCB callback, gpointer data)
{
	MatrixBlock *matrix_block;
	int i, j;
	gboolean breaking = FALSE;

	matrix_block = MATRIX_BLOCK (block);

	for (i = 0; i < matrix_block->p->rows; i++) {
		for (j = 0; j < matrix_block->p->cols; j++) {
			if (callback (block, matrix_block->p->objects[i][j],
				      data)) 
			{
				breaking = TRUE;
				break;
			}
		}

		if (breaking) break;
	}
}

static void
destroy_object_array (MatrixBlock *matrix_block) 
{
	int i, j;

	if (matrix_block->p->objects == NULL) return;

	for (i = 0; i < matrix_block->p->rows; i++) {
		for (j = 0; j < matrix_block->p->cols; j++) {
			if (matrix_block->p->objects[i][j] != NULL)
				gtk_object_unref (GTK_OBJECT
						  (matrix_block->p->
						   objects[i][j]));
		}

		g_free (matrix_block->p->objects[i]);
	}

	g_free (matrix_block->p->objects);
}

static void
setup_object_array (MatrixBlock *matrix_block) 
{
	int i, j;

	if (matrix_block->p->objects != NULL)
		destroy_object_array (matrix_block);

	matrix_block->p->objects = 
		g_new0 (MathObject **, matrix_block->p->rows);

	for (i = 0; i < matrix_block->p->rows; i++)
		matrix_block->p->objects = 
			g_new0 (MathObject *, matrix_block->p->cols);
}
