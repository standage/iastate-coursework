#ifndef HW1_VECTOR
#define HW1_VECTOR

#include <math.h>
#include <stdlib.h>

// 3D vector
typedef struct
{
  float x;
  float y;
  float z;
} Vector;

// Make a copy of the given vector
Vector *vector_copy(Vector *v);

// Compute the normal vector of the plane formed by the two given vectors
Vector *vector_cross_product(Vector *v1, Vector *v2);

// Destructor: free the memory previously used by this vector
void vector_delete(Vector *v);

// Compute the cosine of the angle between the two given vectors
float vector_dot_product(Vector *v1, Vector *v2);

// Constructor: allocate memory for a vector
Vector *vector_new();

// Compute the magnitude of the vector
float vector_norm(Vector *v);

// Normalize this vector
void vector_normalize(Vector *v);

// Rotate this vector theta degrees around the axis k (Rodrigues' rotation formula)
Vector *vector_rotate(Vector *v, Vector *k, float theta);

// Scale this vector by a factor of s
void vector_scale(Vector *v, float s);

#endif
