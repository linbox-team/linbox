
#include "erdos-utils.h"
#include <glib.h>
#include <string.h>

ErdosLineThickness *erdos_line_thickness_new(char *string)
{
  ErdosLineThickness thickness;
  thickness.type = ERDOS_LINE_THICKNESS_MEDIUM;
  thickness.value = 0.0;
  strcpy( thickness.unit, "" );
  if ( string == NULL )
    return NULL;
  else if ( ! strcmp( string, "thick" ) )
    {
      thickness.type = ERDOS_LINE_THICKNESS_THICK;
    }
  else if ( ! strcmp( string, "medium" ) )
    {
      thickness.type = ERDOS_LINE_THICKNESS_MEDIUM;
    }
  else if ( ! strcmp( string, "thin" ) )
    {
      thickness.type = ERDOS_LINE_THICKNESS_THIN;
    }
  else
    {
      int percentage;
      int error = erdos_parse_number_optional_unit(string, &thickness.value, thickness.unit, &percentage);
      if ( error )
	return NULL;
      if ( percentage )
	thickness.type = ERDOS_LINE_THICKNESS_PERCENTAGE;
      else if ( *thickness.unit )
	thickness.type = ERDOS_LINE_THICKNESS_UNIT;
      else
	thickness.type = ERDOS_LINE_THICKNESS_MULTIPLIER;
    }
  return erdos_line_thickness_clone( &thickness );
}

double erdos_line_thickness_calculate_actual(ErdosLineThickness *thickness, ErdosRenderer *renderer, ErdosStyle *style)
{
  if ( thickness == NULL )
    {
      return erdos_renderer_value_lookup(renderer, style, 0x88);
    }
  else
    {
      switch ( thickness->type )
	{
	case ERDOS_LINE_THICKNESS_THICK:
	  return erdos_renderer_value_lookup(renderer, style, 0x88) * 1.5;
	case ERDOS_LINE_THICKNESS_THIN:
	  return erdos_renderer_value_lookup(renderer, style, 0x88) * .5;
	case ERDOS_LINE_THICKNESS_MULTIPLIER:
	  return erdos_renderer_value_lookup(renderer, style, 0x88) * thickness->value;
	case ERDOS_LINE_THICKNESS_PERCENTAGE:
	  return erdos_renderer_value_lookup(renderer, style, 0x88) * thickness->value / 100.0;
	case ERDOS_LINE_THICKNESS_UNIT:
	  return erdos_renderer_unit_convert(renderer, style, thickness->value, thickness->unit);
	case ERDOS_LINE_THICKNESS_MEDIUM: /* roll over. */
	default:
	  return erdos_renderer_value_lookup(renderer, style, 0x88);
	}
    }
}

ErdosLineThickness *erdos_line_thickness_clone(ErdosLineThickness *thickness)
{
  if ( thickness == NULL )
    return NULL;
  else
    {
      ErdosLineThickness *newthickness = g_new(ErdosLineThickness, 1);
      *newthickness = *thickness;
      return newthickness;
    }
}

static void
strip_white(char **string)
{
  char *read_head = *string;
  if ( ! string )
    return;
  while( 1 )
    {
      switch( *read_head )
	{
	case 0x20:
	case 0x09:
	case 0x0a:
	case 0x0d:
	  break;
	default:
	  *string = read_head;
	  return;
	  break;
	}
    }
}

static int
read_double(char **string, double *value)
{
  int return_val = -1;
  *value = 0;
  for ( ; **string; (*string) ++ )
    {
      if ( **string >= '0' && **string <= '9' )
	{
	  *value *= 10.0;
	  *value += **string - '0';
	  return_val = 0;
	}
      else
	break;
    }
  if ( **string == '.' )
    {
      double multiplier = 0.1;
      (*string) ++;
      for ( ; **string; (*string) ++ )
	{
	  if ( **string >= '0' && **string <= '9' )
	    {
	      *value += multiplier * (**string - '0');
	      multiplier /= 10.0;
	      return_val = 0;
	    }
	  else
	    break;
	}
      if ( return_val != 0 )
	(*string) --;
    }
  return return_val;
}

int
erdos_parse_number_optional_unit(char *input, double *value, char unit[3], int *percentage)
{
  int negative = 0;
  int error = 0;
  *percentage = 0;
  *unit = 0;
  strip_white(&input);
  if ( *input == '-' )
    {
      negative = 1;
      input ++;
    }
  error = read_double(&input, value);
  strip_white(&input);
  if ( *input )
    {
      if ( ! strncmp( input, "%", 1 ) )
	{
	  *percentage = 1;
	  input ++;
	  strip_white(&input);
	  if ( *input )
	    return -1;
	  else
	    return 0;
	}
      switch( ((*input) << 8) + *(input + 1) )
	{
	case ('p' << 8) + 't':
	case ('p' << 8) + 'x':
	case ('e' << 8) + 'm':
	case ('e' << 8) + 'x':
	case ('i' << 8) + 'n':
	case ('c' << 8) + 'm':
	case ('m' << 8) + 'm':
	case ('p' << 8) + 'c':
	  strncpy( unit, input, 2 );
	  strip_white(&input);
	  if ( *input )
	    return -1;
	  else
	    return 0;
	default:
	  return -1;
	}
    }
  return 0;
}
