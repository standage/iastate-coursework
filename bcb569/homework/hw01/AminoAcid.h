#ifndef HW1_AMINO_ACID
#define HW1_AMINO_ACID

#include "Array.h"
#include "Atom.h"

// Data structure representing a single peptide
typedef struct
{
  Atom *N;
  Atom *CA;
  Atom *C;
  Array *atoms;
} AminoAcid;

// Associate the given atom with the given amino acid
void amino_acid_add_atom(AminoAcid *aa, Atom *a);

// Destructor: free the memory previously used by this amino acid
void amino_acid_delete(AminoAcid *aa);

// Find an atom based on its name
Atom *amino_acid_find_atom(AminoAcid *aa, const char *atom_name);
Atom *amino_acid_find_atom_pattern(AminoAcid *aa, const char *name_pattern);

// Constructor: allocate memory for this amino acid
AminoAcid *amino_acid_new();

// Number of atoms for this amino acid
int amino_acid_num_atoms(AminoAcid *aa);

// Print side chain torsion angles
void amino_acid_print_chis(AminoAcid *aa, FILE *outstream);

#endif
