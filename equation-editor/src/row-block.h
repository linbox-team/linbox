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

struct _RowBlock 
{
	Block parent;
};

struct _RowBlockClass 
{
	BlockClass block_class;
};

guint row_block_get_type         (void);

GtkObject *row_block_new         (void);

END_GNOME_DECLS

#endif /* __ROW_BLOCK_H */
