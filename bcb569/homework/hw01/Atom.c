#include "Atom.h"

void atom_delete(Atom *a)
{
  free(a);
  a = NULL;
}

Atom *atom_new_from_pdb(const char *pdb_string)
{
  size_t i, j, l;
  char buffer[64];
  Atom *a = (Atom *)malloc( sizeof(Atom) );

  // Parse serial number
  i = 6;
  while(pdb_string[i] == ' ')
    i++;
  l = 10 - i + 1;
  strncpy(buffer, pdb_string + i, l);
  buffer[l] = '\0';
  a->serial_number = atoi(buffer);

  // Parse atom name
  i = 12;
  while(pdb_string[i] == ' ')
    i++;
  j = i;
  while(pdb_string[j+1] != ' ' && j < 16)
    j++;
  l = j - i + 1;
  strncpy(a->atom_name, pdb_string + i, l);
  a->atom_name[l] = '\0';
  
  // Parse alternate location
  a->alternate_location = pdb_string[16];

  // Parse residue name
  strncpy(a->residue_name, pdb_string + 17, 3);
  a->residue_name[3] = '\0';
  
  // Parse chain id
  a->chain_id = pdb_string[21];

  // Parse residue sequence number
  i = 22;
  while(pdb_string[i] == ' ')
    i++;
  l = 25 - i + 1;
  strncpy(buffer, pdb_string + i, l);
  buffer[l] = '\0';
  a->residue_sequence_number = atoi(buffer);

  // Parse insertion code
  a->insertion_code = pdb_string[26];

  // Parse x,y,z coordinates
  strncpy(buffer, pdb_string + 30, 8);
  a->p.x = atof(buffer);
  strncpy(buffer, pdb_string + 38, 8);
  a->p.y = atof(buffer);
  strncpy(buffer, pdb_string + 46, 8);
  a->p.z = atof(buffer);

  // Parse occupancy
  strncpy(buffer, pdb_string + 54, 6);
  a->occupancy = atof(buffer);

  // Parse temp factor
  strncpy(buffer, pdb_string + 60, 6);
  a->temp_factor = atof(buffer);

  // Parse element
  if(pdb_string[76] == ' ')
  {
    a->element[0] = pdb_string[77];
    a->element[1] = '\0';
  }
  else
  {
    strncpy(a->element, pdb_string + 76, 2);
    a->element[2] = '\0';
  }
  
  return a;
}

void atom_print(Atom *a, FILE *outstream)
{
  fprintf( outstream, "ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
           a->serial_number, a->atom_name, a->alternate_location, a->residue_name,
           a->chain_id, a->residue_sequence_number, a->insertion_code, a->p.x, a->p.y, a->p.z,
           a->occupancy, a->temp_factor, a->element );
}
