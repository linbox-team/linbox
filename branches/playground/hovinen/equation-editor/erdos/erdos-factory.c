/* erdos-factory.c
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


#include "erdos-factory.h"
#include "erdos-mn.h"
#include "erdos-mo.h"
#include "erdos-mrow.h"
#include "erdos-mfrac.h"
#include "erdos-mstyle.h"

#include <stdlib.h>
#include <string.h>

/*
 * Get a value for a node either carried as an attibute or as
 * the content of a child.
 */
static char *
xmlGetContent (xmlNodePtr node)
{
	char *ret;
	xmlNodePtr child;

	ret = (char *) xmlNodeGetContent (node);
	if (ret != NULL)
		return ret;
	child = node->childs;

	while (child != NULL) {
		if (!strcmp (child->name, "text")) {
		        /*
			 * !!! Inefficient, but ...
			 */
			ret = xmlGetContent(child);
			if (ret != NULL)
			    return (ret);
		}
		child = child->next;
	}

	return NULL;
}

static char *
strip_white_space(char *string)
{
  char *read_head = string;
  char *write_head = string;
  int last_white = 1 /* TRUE */;
  if ( ! string )
    return string;
  for ( ; *read_head; read_head ++ )
    {
      switch( *read_head )
	{
	case 0x20:
	case 0x09:
	case 0x0a:
	case 0x0d:
	  if ( ! last_white )
	    {
	      *(write_head++) = 0x20;
	    }
	  last_white = 1 /* TRUE */;
	  break;
	default:
	  *(write_head++) = *read_head;
	  last_white = 0 /* FALSE */;
	  break;
	}
    }
  if ( write_head != string && *(write_head - 1) == 0x20 )
    write_head--;
  *write_head = 0;
  return string;
}

ErdosExpression *process_xml(xmlNode *tree)
{
  ErdosExpression *return_val = NULL;
  if ( ! strcmp( tree->name, "mn" ) )
    {
      char *tempstring = strdup( xmlGetContent( tree ) );
      strip_white_space( tempstring );
      return_val = erdos_mn_new( tempstring );
      free( tempstring );
    }
  else if ( ! strcmp( tree->name, "mo" ) )
    {
      char *tempstring = strdup( xmlGetContent( tree ) );
      strip_white_space( tempstring );
      return_val = erdos_mo_new( tempstring );
      free( tempstring );
    }
  else if ( ! strcmp( tree->name, "mstyle" ) )
    {
      xmlNode *children = tree->childs;
      return_val = erdos_mstyle_new(process_xml(children));
    }
  else if ( ! strcmp( tree->name, "mrow" ) )
    {
      ErdosMrow *row = (ErdosMrow *) erdos_mrow_new();
      xmlNode *children;
      for ( children = tree->childs; children; children = children->next )
	{
	  erdos_mrow_append( row, process_xml(children) );
	}
      return_val = (ErdosExpression *)row;
    }
  else if ( ! strcmp( tree->name, "mfrac" ) )
    {
      ErdosExpression *numerator, *denominator;
      xmlNode *children_n;
      xmlNode *children_d;
      children_n = tree->childs;
      if ( children_n )
	children_d = children_n->next;
      if ( children_n && children_d)
	{
	  numerator = process_xml(children_n);
	  denominator = process_xml(children_d);
	  return_val = (ErdosExpression *) erdos_mfrac_new(numerator, denominator);
	}
      else
	  return_val = NULL;
    }
  if ( return_val && tree->properties )
    {
      xmlAttr *attribute = tree->properties;
      while ( attribute )
	{
	  char *name = attribute->name;
	  char *val = xmlGetContent(attribute->val);

	  erdos_expression_set_attribute(return_val, name, val);
	  
	  attribute = attribute->next;
	}
    }
  return return_val;
}
