#include "Vector.h"

Vector *vector_copy(Vector *v)
{
  Vector *u = vector_new();
  u->x = v->x;
  u->y = v->y;
  u->z = v->z;
  return u;
}

Vector *vector_cross_product(Vector *v1, Vector *v2)
{
  Vector *u = vector_new();
  u->x = (v1->y*v2->z) - (v1->z*v2->y);
  u->y = (v1->z*v2->x) - (v1->x*v2->z);
  u->z = (v1->x*v2->y) - (v1->y*v2->x);
  return u;
}

void vector_delete(Vector *v)
{
  free(v);
  v = NULL;
}

float vector_dot_product(Vector *v1, Vector *v2)
{
  return (v1->x*v2->x) + (v1->y*v2->y) + (v1->z*v2->z);
}

Vector *vector_new()
{
  Vector *v = (Vector *)malloc( sizeof(Vector) );
  return v;
}

float vector_norm(Vector *v)
{
  return sqrt(vector_dot_product(v, v));
}

void vector_normalize(Vector *v)
{
  vector_scale(v, 1/vector_norm(v));
}

Vector *vector_rotate(Vector *v, Vector *k, float theta)
{
  float cos_theta = cos(theta);
  float sin_theta = sin(theta);
  
  Vector *k_cross_v = vector_cross_product(k, v);
  vector_scale(k_cross_v, sin_theta);
  
  Vector *k_k_dot_v = vector_copy(k);
  vector_scale(k_k_dot_v, vector_dot_product(k, v)*(1 - cos_theta));
  
  Vector *v_rot = vector_new();
  v_rot->x = (v->x * cos(theta)) + k_cross_v->x + k_k_dot_v->x;
  v_rot->y = (v->y * cos(theta)) + k_cross_v->y + k_k_dot_v->y;
  v_rot->z = (v->z * cos(theta)) + k_cross_v->z + k_k_dot_v->z;
  
  vector_delete(k_cross_v);
  vector_delete(k_k_dot_v);
  
  return v_rot;
}

void vector_scale(Vector *v, float s)
{
  v->x = (v->x * s);
  v->y = (v->y * s);
  v->z = (v->z * s);
}
