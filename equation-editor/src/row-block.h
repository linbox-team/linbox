/* -*- mode: c; style: linux -*- */

/* row-block.h
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

#ifndef __ROW_BLOCK_H
#define __ROW_BLOCK_H

#include <gnome.h>

#include "block.h"

BEGIN_GNOME_DECLS

#define ROW_BLOCK(obj)          GTK_CHECK_CAST (obj, row_block_get_type (), RowBlock)
#define ROW_BLOCK_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, row_block_get_type (), RowBlockClass)
#define IS_ROW_BLOCK(obj)       GTK_CHECK_TYPE (obj, row_block_get_type ())

typedef struct _RowBlock RowBlock;
typedef struct _RowBlockClass RowBlockClass;
typedef struct _RowBlockPrivate RowBlockPrivate;

struct _RowBlock 
{
	Block parent;

	RowBlockPrivate *p;
};

struct _RowBlockClass 
{
	BlockClass block_class;
};

guint row_block_get_type            (void);

/**
 * row_block_new:
 * @void: 
 * 
 * Factory method
 * 
 * Return value: New empty row block object
 **/

GtkObject *row_block_new            (void);

/**
 * row_block_insert:
 * @row_block: 
 * @math_object: 
 * @before: The object before which to insert, NULL to insert at the end
 * 
 * Insert the given math object before the object `before'; append to the end
 * if the `before' is not located in the row block
 **/

void row_block_insert               (RowBlock *row_block,
				     MathObject *math_object,
				     MathObject *before);

/**
 * row_block_insert_at:
 * @row_block: 
 * @math_object: 
 * @position: 
 * 
 * Insert a math object at the given position
 **/

void row_block_insert_at            (RowBlock *row_block,
				     MathObject *math_object,
				     gint position);

/**
 * row_block_delete_at:
 * @row_block: 
 * @position: The position whose object to delete
 * 
 * Deletes the math object at the given position (the math object ptr is
 * included to verify that the math object being deleted is correct.
 **/

void row_block_delete_at            (RowBlock *row_block,
				     gint position);


/**
 * row_block_get_object_at:
 * @row_block: 
 * @position: Position at which to get object (0 <= position < length)
 * 
 * Get the object at the specified position
 * 
 * Return value: 
 **/

MathObject *row_block_get_object_at (RowBlock *row_block,
				     gint position);

/**
 * row_block_get_position_of:
 * @row_block: 
 * @object: 
 * 
 * Get the position of the given object in the row block
 * 
 * Return value: The position of the given object in the row block, or -1 if
 * the object does not exist
 **/

gint row_block_get_position_of      (RowBlock *row_block,
				     MathObject *object);

/**
 * row_block_get_length:
 * @row_block: 
 * 
 * Get the number of objects in the row block
 * 
 * Return value: 
 **/

guint row_block_get_length          (RowBlock *row_block);

END_GNOME_DECLS

#endif /* __ROW_BLOCK_H */
