/* -*- mode: c; style: linux -*- */

/* cursor.h
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

#ifndef __CURSOR_H
#define __CURSOR_H

#include <gnome.h>

#include "math-expression.h"

BEGIN_GNOME_DECLS

#define CURSOR(obj)          GTK_CHECK_CAST (obj, cursor_get_type (), Cursor)
#define CURSOR_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, cursor_get_type (), CursorClass)
#define IS_CURSOR(obj)       GTK_CHECK_TYPE (obj, cursor_get_type ())

typedef struct _Cursor Cursor;
typedef struct _CursorClass CursorClass;
typedef struct _CursorPrivate CursorPrivate;

struct _Cursor 
{
	GtkObject parent;

	CursorPrivate *p;
};

struct _CursorClass 
{
	GtkObjectClass gtk_object_class;
};

guint       cursor_get_type                      (void);

GtkObject  *cursor_new                           (MathExpression *expr);

MathObject *cursor_get_current_object            (Cursor *cursor);
gint        cursor_get_insertion_point           (Cursor *cursor);
MathObject *cursor_get_object_at_insertion_point (Cursor *cursor);

void        cursor_move_left                     (Cursor *cursor);
void        cursor_move_right                    (Cursor *cursor);

void        cursor_move_to_beginning             (Cursor *cursor);
void        cursor_move_to_end                   (Cursor *cursor);

END_GNOME_DECLS

#endif /* __CURSOR_H */
