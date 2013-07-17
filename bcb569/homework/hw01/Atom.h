#ifndef HW1_ATOM
#define HW1_ATOM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Geometry.h"

// Field names for this data structure are based on PDB documentation for ATOM records
// PDB format documentation: http://www.wwpdb.org/documentation/format33/sect9.html
typedef struct
{
  int   serial_number;
  char  atom_name[5];
  char  alternate_location;
  char  residue_name[4];
  char  chain_id;
  int   residue_sequence_number;
  char  insertion_code;
  Point p;
  float occupancy;
  float temp_factor;
  char  element[3];
} Atom;

// Destructor: free memory previously used by this atom
void atom_delete(Atom *a);

// Constructor: allocate memory for this atom, set its values based on the provided PDB `ATOM' entry
Atom *atom_new_from_pdb(const char *pdb_string);

// Print `ATOM' record to the given file in PDB format
void atom_print(Atom *a, FILE *outstream);

#endif
