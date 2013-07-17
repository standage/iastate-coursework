#ifndef HW1_GEOMETRY
#define HW1_GEOMETRY

#include <math.h>
#include "Vector.h"

// 3D point
typedef struct
{
  float x;
  float y;
  float z;
} Point;

// Determine the Euclidean distance between the 3D points
float euclidean_distance(Point *p, Point *q);

// Determine the angle between the points/vectors; solution based on answer found on thread at
// http://stackoverflow.com/questions/1211212/how-to-calculate-an-angle-from-three-points
float get_angle_from_points(Point *p, Point *q, Point *r);
float get_angle_from_vectors(Vector *v1, Vector *v2);

// Determine the angle between two planes--the first plane defined by (p1, p2, p3)
// and the second plane defined by (p2, p3, p4)
float get_dihedral_angle(Point *p1, Point *p2, Point *p3, Point *p4);

#endif
