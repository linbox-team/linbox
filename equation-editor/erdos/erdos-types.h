
#ifndef __ERDOS_TYPES_H__
#define __ERDOS_TYPES_H__

typedef enum _ErdosForm ErdosForm;
typedef enum _ErdosLineThicknessType ErdosLineThicknessType;
typedef struct _ErdosLineThickness ErdosLineThickness;


enum _ErdosForm
{
  ERDOS_PREFIX,
  ERDOS_INFIX,
  ERDOS_SUFFIX
};


enum _ErdosLineThicknessType
{
  ERDOS_LINE_THICKNESS_MULTIPLIER,
  ERDOS_LINE_THICKNESS_PERCENTAGE,
  ERDOS_LINE_THICKNESS_UNIT,
  ERDOS_LINE_THICKNESS_THIN,
  ERDOS_LINE_THICKNESS_MEDIUM,
  ERDOS_LINE_THICKNESS_THICK
};

struct _ErdosLineThickness
{
  ErdosLineThicknessType type;
  double value;
  char unit[3];
};

#endif /* __ERDOS_TYPES_H__ */
