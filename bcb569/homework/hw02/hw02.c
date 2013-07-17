#include <stdlib.h>
#include <time.h>
#include "Array.h"
#include "Atom.h"

#define PI 3.14159265
#define NUM_ARGS 2

int main(int argc, char **argv)
{
  // Check command-line arguments
  if(argc != NUM_ARGS)
  {
    fprintf( stderr, "Error: expected %d arguments, got %d\n", NUM_ARGS - 1,
             argc - 1 );
    return 1;
  }
  
  // Seed the random number generator
  srand( time(NULL) );
  
  // Open input file
  const char *infile = argv[1];
  FILE *instream = fopen(infile, "r");
  if(instream == NULL)
  {
    fprintf(stderr, "file fail\n");
    return 1;
  }
  
  // Load atoms from PDB file, sample points from the surface of each atom
  Array *atoms = array_new( sizeof(Atom *) );
  char buffer[1024] = "";
  while(fgets(buffer, 1024, instream))
  {
    if(strlen(buffer) >= 4 && strncmp(buffer, "ATOM", 4) == 0)
    {
      Atom *a = parse_atom_from_pdb(buffer);
      array_add(atoms, a);
    }
  }
  fclose(instream);
  int num_atoms  = array_size(atoms);
  int num_points = num_atoms * POINTS_PER_ATOM;
  fprintf( stderr, "Sampled %d points (%d atoms, %d points per atom)\n",
           num_points, num_atoms, POINTS_PER_ATOM );

  // Classify each point as accessible or buried, compute approximate surface
  // area from accessible points
  Atom *last_atom = (Atom *)array_get(atoms, num_atoms - 1);
  int num_residues = last_atom->residue_sequence_number;
  double *sasa_per_residue = (double *)malloc( sizeof(double)*num_residues );
  fputs("Finding accessible points\n", stderr);
  int i;
  for(i = 0; i < num_atoms; i++)
  {
    Atom *current_atom = (Atom *)array_get(atoms, i);
    int j;
    for(j = 0; j < array_size(current_atom->surface); j++)
    {
      Point *p = (Point *)array_get(current_atom->surface, j);
      int k;
      int buried = 0;
      for(k = 0; k < num_atoms; k++)
      {
        if(k != i)
        {
          Atom *other_atom = (Atom *)array_get(atoms, k);
          double distance = euclidean_distance(p, &other_atom->p);
          if(distance < other_atom->radius)
          {
            buried = 1;
            break;
          }
        }
      }
      if(!buried)
        array_add(current_atom->accessible, p);
    }
    if((i+1) % 50 == 0)
      fprintf(stderr, "    %3d atoms done\n", i+1);
  }
  
  // Calculate residue SASAs from atom SASAs, and the whole protein SASA from
  // residue SASAs
  double sasa_protein = 0.0;
  int j = 0;
  for(i = 0; i < num_residues; i++)
  {
    Atom *a = array_get(atoms, j);
    while(a->residue_sequence_number == i+1 && j < num_atoms)
    {
      double c = 4.0 * PI / POINTS_PER_ATOM;
      double num_accessible = (double)array_size(a->accessible);
      a->sasa = c * num_accessible * a->radius * a->radius;
      sasa_per_residue[i] += a->sasa;
      sasa_protein += a->sasa;
      if(++j < num_atoms)
        a = array_get(atoms, j);
    }
  }
  
  // Print results
  j = 0;
  printf( "Solvent-accessible surface area for '%s': %.4lf\n", infile,
          sasa_protein );
  for(i = 0; i < num_residues; i++)
  {
    printf("  Residue %2d: %.4lf\n", i+1, sasa_per_residue[i]);
    Atom *a = array_get(atoms, j);
    while(a->residue_sequence_number == i+1 && j < num_atoms)
    {
      printf("    Atom %3d: %.4lf\n", j+1, a->sasa);
      if(++j < num_atoms)
        a = array_get(atoms, j);
    }
  }

  // Free memory
  free(sasa_per_residue);
  for(i = 0; i < array_size(atoms); i++)
  {
    Atom *a = (Atom *)array_get(atoms, i);
    atom_delete(a);
  }
  array_delete(atoms);
  return 0;
}
