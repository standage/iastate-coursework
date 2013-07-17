#include "AminoAcid.h"

void amino_acid_add_atom(AminoAcid *aa, Atom *a)
{
  array_add(aa->atoms, a);
  if(strcmp(a->atom_name, "N") == 0)
    aa->N = a;
  else if(strcmp(a->atom_name, "CA") == 0)
    aa->CA = a;
  else if(strcmp(a->atom_name, "C") == 0)
    aa->C = a;
}

void amino_acid_delete(AminoAcid *aa)
{
  int i;
  for(i = 0; i < array_size(aa->atoms); i++)
  {
    Atom *a = (Atom *)array_get(aa->atoms, i);
    atom_delete(a);
  }
  array_delete(aa->atoms);
  aa->N  = NULL;
  aa->CA = NULL;
  aa->C  = NULL;
  free(aa);
}

Atom *amino_acid_find_atom(AminoAcid *aa, const char *atom_name)
{
  int i;
  for(i = 0; i < array_size(aa->atoms); i++)
  {
    Atom *a = (Atom *)array_get(aa->atoms, i);
    if(strcmp(a->atom_name, atom_name) == 0)
      return a;
  }
  return NULL;
}

Atom *amino_acid_find_atom_pattern(AminoAcid *aa, const char *name_pattern)
{
  int i;
  for(i = 0; i < array_size(aa->atoms); i++)
  {
    Atom *a = (Atom *)array_get(aa->atoms, i);
    int name_length = strlen(a->atom_name);
    int j;
    for(j = 0; j < name_length; j++)
    {
      if(strcmp(a->atom_name + j, name_pattern) == 0)
        return a;
    }
  }
  return NULL;
}

AminoAcid *amino_acid_new()
{
  AminoAcid *aa = (AminoAcid *)malloc( sizeof(AminoAcid) );
  aa->N  = NULL;
  aa->CA = NULL;
  aa->C  = NULL;
  aa->atoms = array_new();
  return aa;
}

int amino_acid_num_atoms(AminoAcid *aa)
{
  return array_size(aa->atoms);
}

void amino_acid_print_chis(AminoAcid *aa, FILE *outstream)
{
  Atom *N  = amino_acid_find_atom(aa, "N");
  Atom *CA = amino_acid_find_atom(aa, "CA");
  Atom *CB = amino_acid_find_atom(aa, "CB");
  Atom *starG = amino_acid_find_atom_pattern(aa, "G");
  if(starG != NULL)
  {
    float chi1 = get_dihedral_angle(&N->p, &CA->p, &CB->p, &starG->p);
    printf("        %-6s %.3f\n", "X_1:", chi1);
  }
  
  Atom *starG1 = amino_acid_find_atom_pattern(aa, "G1");
  if(starG1 != NULL)
  {
    float chi11 = get_dihedral_angle(&N->p, &CA->p, &CB->p, &starG1->p);
    printf("        %-6s %.3f\n", "X_1,1:", chi11);
  }
  
  Atom *starG2 = amino_acid_find_atom_pattern(aa, "G2");
  if(starG2 != NULL)
  {
    float chi12 = get_dihedral_angle(&N->p, &CA->p, &CB->p, &starG2->p);
    printf("        %-6s %.3f\n", "X_1,2:", chi12);
  }
  
  Atom *starD = amino_acid_find_atom_pattern(aa, "D");
  if(starG != NULL && starD != NULL)
  {
    float chi2 = get_dihedral_angle(&CA->p, &CB->p, &starG->p, &starD->p);
    printf("        %-6s %.3f\n", "X_2:", chi2);
  }
  
  Atom *starD1 = amino_acid_find_atom_pattern(aa, "D1");
  if(starG != NULL && starD1 != NULL)
  {
    float chi21 = get_dihedral_angle(&CA->p, &CB->p, &starG->p, &starD1->p);
    printf("        %-6s %.3f\n", "X_2,1:", chi21);
  }
  
  Atom *starD2 = amino_acid_find_atom_pattern(aa, "D2");
  if(starG != NULL && starD2 != NULL)
  {
    float chi22 = get_dihedral_angle(&CA->p, &CB->p, &starG->p, &starD2->p);
    printf("        %-6s %.3f\n", "X_2,2:", chi22);
  }
  
  Atom *starE = amino_acid_find_atom_pattern(aa, "E");
  if(starG != NULL && starD != NULL && starE != NULL)
  {
    float chi3 = get_dihedral_angle(&CB->p, &starG->p, &starD->p, &starE->p);
    printf("        %-6s %.3f\n", "X_3:", chi3);
  }
  
  Atom *starE1 = amino_acid_find_atom_pattern(aa, "E1");
  if(starG != NULL && starD != NULL && starE1 != NULL)
  {
    float chi31 = get_dihedral_angle(&CB->p, &starG->p, &starD->p, &starE1->p);
    printf("        %-6s %.3f\n", "X_3,1:", chi31);
  }
  
  Atom *starE2 = amino_acid_find_atom_pattern(aa, "E2");
  if(starG != NULL && starD != NULL && starE2 != NULL)
  {
    float chi32 = get_dihedral_angle(&CB->p, &starG->p, &starD->p, &starE2->p);
    printf("        %-6s %.3f\n", "X_3,2:", chi32);
  }
  
  Atom *starZ = amino_acid_find_atom_pattern(aa, "Z");
  if(starG != NULL && starD != NULL && starE != NULL && starZ != NULL)
  {
    float chi4 = get_dihedral_angle(&starG->p, &starD->p, &starE->p, &starZ->p);
    printf("        %-6s %.3f\n", "X_4:", chi4);
  }
  
  Atom *NH1 = amino_acid_find_atom(aa, "NH1");
  if(starD != NULL && starE != NULL && starZ != NULL && NH1 != NULL)
  {
    float chi51 = get_dihedral_angle(&starD->p, &starE->p, &starZ->p, &NH1->p);
    printf("        %-6s %.3f\n", "X_5,1:", chi51);
  }
  
  Atom *NH2 = amino_acid_find_atom(aa, "NH2");
  if(starD != NULL && starE != NULL && starZ != NULL && NH2 != NULL)
  {
    float chi52 = get_dihedral_angle(&starD->p, &starE->p, &starZ->p, &NH2->p);
    printf("        %-6s %.3f\n", "X_5,2:", chi52);
  }
}
