/* -*- mode: c; style: linux -*- */

/* controller.h
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *            Matt Spilich, Anthony Asher, Rob Wede
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

#ifndef __CONTROLLER_H
#define __CONTROLLER_H

#include <gnome.h>
#include "math-object.h"

BEGIN_GNOME_DECLS

#define CONTROLLER(obj)          GTK_CHECK_CAST (obj, controller_get_type (), Controller)
#define CONTROLLER_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, controller_get_type (), ControllerClass)
#define IS_CONTROLLER(obj)       GTK_CHECK_TYPE (obj, controller_get_type ())

typedef struct _Controller Controller;
typedef struct _ControllerClass ControllerClass;
typedef struct _ControllerPrivate ControllerPrivate;

struct _Controller 
{
	GtkObject parent;

	ControllerPrivate *p;
};

struct _ControllerClass 
{
	GtkObjectClass gtk_object_class;
};

guint controller_get_type         (void);

GtkObject *controller_new         (void);

void controller_insert	(Controller *controller, 
			GdkEventKey *event);

void controller_initialize (Controller *controller, MathObject *toplevel);

END_GNOME_DECLS

#endif /* __CONTROLLER_H */
