/* erdos.c
 * Copyright (C) 1999  Chris Lahey <clahey@umich.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "erdos.h"
#include <glib.h>

ErdosObjectClass *
erdos_object_get_class (void)
{
  static ErdosObjectClass *object_type = NULL;

  if (!object_type)
    {
      object_type = g_malloc(sizeof(ErdosObjectClass));
      erdos_object_class_init(object_type);
    }

  return object_type;
}

void
erdos_object_class_init (ErdosObjectClass *klass)
{
  klass->finalize = NULL;
  klass->clone = NULL;
}

void
erdos_object_init (ErdosObject *object)
{
  ((ErdosObjectClass *)(((ErdosObject *)object)->klass)) = erdos_object_get_class();
  object->ref_count = 0;
}

void erdos_object_ref (ErdosObject *object)
{
  if ( object )
    {
      object->ref_count ++;
    }
}

void erdos_object_unref (ErdosObject *object)
{
  if ( object )
    {
      object->ref_count --;
      if ( object->ref_count <= 0 )
	{
	  if ( object->klass->finalize )
	    {
	      object->klass->finalize(object);
	    }
	  g_free( object );
	}
    }
}

ErdosObject *erdos_object_cow (ErdosObject *object)
{
  if ( object->ref_count > 1 && object->klass->clone )
    {
      ErdosObject *newobj = object->klass->clone(object);
      object->ref_count --;
      newobj->ref_count = 1;
      return newobj;
    }
  else
    return object;
}
