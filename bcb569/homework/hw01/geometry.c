#include "Geometry.h"

float euclidean_distance(Point *p, Point *q)
{
  float sum_squared_distances = 0.0;

  sum_squared_distances += pow(q->x - p->x, 2);
  sum_squared_distances += pow(q->y - p->y, 2);
  sum_squared_distances += pow(q->z - p->z, 2);

  return sqrt(sum_squared_distances);
}

float get_angle_from_points(Point *p, Point *q, Point *r)
{
  Vector v1 = { p->x - q->x, p->y - q->y, p->z - q->z };
  Vector v2 = { r->x - q->x, r->y - q->y, r->z - q->z };
  return get_angle_from_vectors(&v1, &v2);
}

float get_angle_from_vectors(Vector *v1, Vector *v2)
{
  float dp = vector_dot_product(v1, v2);
  float v1n = vector_norm(v1);
  float v2n = vector_norm(v2);
  
  float angle_radians = acos(dp/(v1n * v2n));
  float angle_degrees = angle_radians * 180 / 3.14159265;
  return angle_degrees;
}

float get_dihedral_angle(Point *p1, Point *p2, Point *p3, Point *p4)
{
  Vector v1 = { p1->x - p2->x, p1->y - p2->y, p1->z - p2->z };
  Vector v2 = { p3->x - p2->x, p3->y - p2->y, p3->z - p2->z };
  Vector *n1 = vector_cross_product(&v1, &v2);
  
  Vector v3 = { p2->x - p3->x, p2->y - p3->y, p2->z - p3->z };
  Vector v4 = { p4->x - p3->x, p4->y - p3->y, p4->z - p3->z };
  Vector *n2 = vector_cross_product(&v3, &v4);
  
  float angle = get_angle_from_vectors(n1, n2);
  vector_delete(n1);
  vector_delete(n2);
  
  return angle;
}
