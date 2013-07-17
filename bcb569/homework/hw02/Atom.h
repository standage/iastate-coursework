#ifndef HW2_ATOM
#define HW2_ATOM

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "Array.h"

#define POINTS_PER_ATOM 500

// van des Waals radii as published by Bondi
#define PROBE_RADIUS    1.40
#define VDWR_CARBON     1.70
#define VDWR_HYDROGEN   1.20
#define VDWR_NITROGEN   1.55
#define VDWR_OXYGEN     1.52
#define VDWR_SULFUR     1.83

// Simple data structure for 3D points
typedef struct
{
  double x;
  double y;
  double z;
} Point;

// Data structure for aggregating relevant data for each atom in the protein
typedef struct
{
  Point p;
  double radius; // the effective radius, including that of the probe
  char element;
  int residue_sequence_number;
  double sasa;
  Array *surface;
  Array *accessible;
} Atom;

// Distance between the two points
double euclidean_distance(Point *p, Point *q);

// Free memory previously used by this atom object
void atom_delete(Atom *a);

// Randomly sample n points uniformly from this atom's van der Waals surface
void atom_sample_surface(Atom *a, int n);

// Allocate memory for a new atom object
Atom *atom_new(Point *p, char element, int rsn);

// Parse relevant information from an ATOM record in PDB format 
Atom *parse_atom_from_pdb(const char *line);

#endif
