
#ifndef __ERDOS_UTILS_H__
#define __ERDOS_UTILS_H__

#include "erdos-types.h"
#include "erdos-renderer.h"
#include "erdos-style.h"


ErdosLineThickness *erdos_line_thickness_new(char *string);
ErdosLineThickness *erdos_line_thickness_clone(ErdosLineThickness *thickness);
double erdos_line_thickness_calculate_actual(ErdosLineThickness *thickness, ErdosRenderer *renderer, ErdosStyle *style);

/* 0 if successful */
/* percentage is a true false value of whether to divide a value without units by 100. */
int erdos_parse_number_optional_unit(char *input, double *value, char unit[3], int *percentage);

#endif /* __ERDOS_UTILS_H__ */
