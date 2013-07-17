/*
  Daniel S. Standage
  BCB 569: Structural genome informatics
  Fall 2011
  Homework 1
*/
#include <assert.h>
#include <stdio.h>

#include "AminoAcid.h"

#define EXPECTED_ARG_COUNT 3
#define BUFFER_SIZE 128
#define false 0
#define true 1
typedef int bool;


int main(int argc, char **argv)
{
  int i;

  // Grab and check PDB file from command-line arguments
  if(argc != EXPECTED_ARG_COUNT)
  {
    fprintf(stderr, "Error: expected %d argument(s), found %d\n", EXPECTED_ARG_COUNT - 1, argc - 1);
    return 1;
  }
  FILE *instream = fopen(argv[1], "r");
  if(instream == NULL)
  {
    fprintf(stderr, "Error: could not open PDB file '%s'\n", argv[1]);
    return 1;
  }

  // Parse atoms from PDB file
  Array *polypeptide = array_new();
  AminoAcid *aa = amino_acid_new();
  char buffer[BUFFER_SIZE];
  bool first_n = true;
  while(fgets(buffer, BUFFER_SIZE, instream))
  {
    if(strlen(buffer) > 4 && strncmp(buffer, "ATOM", 4) == 0)
    {
      Atom *a = atom_new_from_pdb(buffer);
      if(strcmp(a->atom_name, "N")  == 0)
      {
        if(first_n)
          first_n = false;
        else
        {
          array_add(polypeptide, aa);
          aa = amino_acid_new();
        }
      }
      amino_acid_add_atom(aa, a);
    }
  }
  array_add(polypeptide, aa);
  fclose(instream);
  size_t n = array_size(polypeptide);
  printf("\nFound n = %lu residues\n", n);

  // Calculate bond lengths, bond angles, alpha carbon distances
  printf("Calculating the following in a single pass:\n");
  printf("    3n-1 = %lu bond lengths\n", 3*n - 1);
  printf("    3n-2 = %lu bond angles\n", 3*n - 2);
  printf("     n-1 = %lu alpha carbon distances\n\n", n - 1);
  
  float total_bond_length_N_CA = 0.0f;
  float total_bond_length_CA_C = 0.0f;
  float total_bond_length_C_N  = 0.0f;
  float *bond_lengths_N_CA = (float *)malloc(n*sizeof(float));
  float *bond_lengths_CA_C = (float *)malloc(n*sizeof(float));
  float *bond_lengths_C_N  = (float *)malloc((n-1)*sizeof(float));
  size_t current_bond_length_N_CA = 0;
  size_t current_bond_length_CA_C = 0;
  size_t current_bond_length_C_N  = 0;
  
  float total_bond_angle_N_CA_C = 0.0f;
  float total_bond_angle_CA_C_N = 0.0f;
  float total_bond_angle_C_N_CA = 0.0f;
  float *bond_angles_N_CA_C = (float *)malloc(n*sizeof(float));
  float *bond_angles_CA_C_N = (float *)malloc((n-1)*sizeof(float));
  float *bond_angles_C_N_CA = (float *)malloc((n-1)*sizeof(float));
  size_t current_bond_angle_N_CA_C = 0;
  size_t current_bond_angle_CA_C_N = 0;
  size_t current_bond_angle_C_N_CA = 0;
  
  float total_carbon_distance = 0.0f;
  float *carbon_distances = (float *)malloc((n - 1)*sizeof(float));
  size_t current_carbon_distance = 0;
  for(i = 0; i < n - 1; i++)
  {
    AminoAcid *p_i = (AminoAcid *)array_get(polypeptide, i);
    AminoAcid *p_j = (AminoAcid *)array_get(polypeptide, i+1);
    
    bond_lengths_N_CA[current_bond_length_N_CA] = euclidean_distance(&p_i->N->p, &p_i->CA->p);
    total_bond_length_N_CA += bond_lengths_N_CA[current_bond_length_N_CA++];
    bond_lengths_CA_C[current_bond_length_CA_C] = euclidean_distance(&p_i->CA->p, &p_i->C->p);
    total_bond_length_CA_C += bond_lengths_CA_C[current_bond_length_CA_C++];
    bond_lengths_C_N[current_bond_length_C_N] = euclidean_distance(&p_i->C->p, &p_j->N->p);
    total_bond_length_C_N += bond_lengths_C_N[current_bond_length_C_N++];

    bond_angles_N_CA_C[current_bond_angle_N_CA_C] = get_angle_from_points(&p_i->N->p, &p_i->CA->p, &p_i->C->p);
    total_bond_angle_N_CA_C += bond_angles_N_CA_C[current_bond_angle_N_CA_C++];
    bond_angles_CA_C_N[current_bond_angle_CA_C_N] = get_angle_from_points(&p_i->CA->p, &p_i->C->p, &p_j->N->p);
    total_bond_angle_CA_C_N += bond_angles_CA_C_N[current_bond_angle_CA_C_N++];
    bond_angles_C_N_CA[current_bond_angle_C_N_CA] = get_angle_from_points(&p_i->C->p, &p_j->N->p, &p_j->CA->p);
    total_bond_angle_C_N_CA += bond_angles_C_N_CA[current_bond_angle_C_N_CA++];

    carbon_distances[current_carbon_distance] = euclidean_distance(&p_i->CA->p, &p_j->CA->p);
    total_carbon_distance += carbon_distances[current_carbon_distance++];

    if(i == n - 2)
    {
      bond_lengths_N_CA[current_bond_length_N_CA] = euclidean_distance(&p_j->N->p, &p_j->CA->p);
      total_bond_length_N_CA += bond_lengths_N_CA[current_bond_length_N_CA++];
      bond_lengths_CA_C[current_bond_length_CA_C] = euclidean_distance(&p_j->CA->p, &p_j->C->p);
      total_bond_length_CA_C += bond_lengths_CA_C[current_bond_length_CA_C++];

      bond_angles_N_CA_C[current_bond_angle_N_CA_C] = get_angle_from_points(&p_j->N->p, &p_j->CA->p, &p_j->C->p);
      total_bond_angle_N_CA_C += bond_angles_N_CA_C[current_bond_angle_N_CA_C++];
    }
  }
  // Sanity check
  assert(current_bond_length_N_CA + current_bond_length_CA_C + current_bond_length_C_N == (3*n - 1));
  assert(current_bond_angle_N_CA_C + current_bond_angle_CA_C_N + current_bond_angle_C_N_CA == (3*n - 2));
  assert(current_carbon_distance == (n - 1));

  // Print results for bond lengths, bond angles, alpha carbon distances
  puts("Results:                         mean          std. dev.");

  puts("    Bond lengths:");
  float mean_bond_length = total_bond_length_N_CA / n;
  float total_bond_length_squared_diff = 0.0f;
  for(i = 0; i < current_bond_length_N_CA; i++)
    total_bond_length_squared_diff += pow(bond_lengths_N_CA[i] - mean_bond_length, 2);
  float bond_length_sd = sqrt(total_bond_length_squared_diff / n);
  printf("        N_i - CA_i:              %.3fÅ        %.3fÅ\n", mean_bond_length, bond_length_sd);
  mean_bond_length = total_bond_length_CA_C / n;
  total_bond_length_squared_diff = 0.0f;
  for(i = 0; i < current_bond_length_CA_C; i++)
    total_bond_length_squared_diff += pow(bond_lengths_CA_C[i] - mean_bond_length, 2);
  bond_length_sd = sqrt(total_bond_length_squared_diff / n);
  printf("        CA_i - C_i:              %.3fÅ        %.3fÅ\n", mean_bond_length, bond_length_sd);
  mean_bond_length = total_bond_length_C_N / n;
  total_bond_length_squared_diff = 0.0f;
  for(i = 0; i < current_bond_length_C_N; i++)
    total_bond_length_squared_diff += pow(bond_lengths_C_N[i] - mean_bond_length, 2);
  bond_length_sd = sqrt(total_bond_length_squared_diff / n);
  printf("        C_i - N_i+1:             %.3fÅ        %.3fÅ\n", mean_bond_length, bond_length_sd);

  puts("    Bond angles:");
  float mean_bond_angle = total_bond_angle_N_CA_C / n;
  float total_bond_angle_squared_diff = 0.0f;
  for(i = 0; i < current_bond_angle_N_CA_C; i++)
    total_bond_angle_squared_diff += pow(bond_angles_N_CA_C[i] - mean_bond_angle, 2);
  float bond_angle_sd = sqrt(total_bond_angle_squared_diff / n);
  printf("        N_i - CA_i - C_i:        %.3fº      %.3fº\n", mean_bond_angle, bond_angle_sd);
  mean_bond_angle = total_bond_angle_CA_C_N / n;
  total_bond_angle_squared_diff = 0.0f;
  for(i = 0; i < current_bond_angle_CA_C_N; i++)
    total_bond_angle_squared_diff += pow(bond_angles_CA_C_N[i] - mean_bond_angle, 2);
  bond_angle_sd = sqrt(total_bond_angle_squared_diff / n);
  printf("        CA_i - C_i - N_i+1:      %.3fº      %.3fº\n", mean_bond_angle, bond_angle_sd);
  mean_bond_angle = total_bond_angle_C_N_CA / n;
  total_bond_angle_squared_diff = 0.0f;
  for(i = 0; i < current_bond_angle_C_N_CA; i++)
    total_bond_angle_squared_diff += pow(bond_angles_C_N_CA[i] - mean_bond_angle, 2);
  bond_angle_sd = sqrt(total_bond_angle_squared_diff / n);
  printf("        C_i - N_i+1 - CA_i+1:    %.3fº      %.3fº\n", mean_bond_angle, bond_angle_sd);

  puts("    Alpha carbon distances:");
  float mean_carbon_distance = total_carbon_distance / (n - 1);
  float total_carbon_distance_squared_diff = 0.0f;
  for(i = 0; i < current_carbon_distance; i++)
    total_carbon_distance_squared_diff += pow(carbon_distances[i] - mean_carbon_distance, 2);
  float carbon_distance_sd = sqrt(total_carbon_distance_squared_diff / (n - 1));
  printf("        CA_i - CA_i+1:           %.3fÅ        %.3fÅ\n\n", mean_carbon_distance, carbon_distance_sd);

  // Compute angles for residue 30
  AminoAcid *p29 = (AminoAcid *)array_get(polypeptide, 28);
  AminoAcid *p30 = (AminoAcid *)array_get(polypeptide, 29);
  AminoAcid *p31 = (AminoAcid *)array_get(polypeptide, 30);
  printf("Residue %d: %s (%d atoms)\n", p30->N->residue_sequence_number, p30->N->residue_name, amino_acid_num_atoms(p30));
  float phi = get_dihedral_angle(&p29->C->p, &p30->N->p, &p30->CA->p, &p30->C->p);
  float psi = get_dihedral_angle(&p30->N->p, &p30->CA->p, &p30->C->p, &p31->N->p);
  float omega = get_dihedral_angle(&p30->CA->p, &p30->C->p, &p31->N->p, &p31->CA->p);
  printf("    φ: %.3fº\n", phi);
  printf("    ψ: %.3fº\n", psi);
  printf("    ω: %.3fº\n\n", omega);

  // Set phi to 0, rotate all affected atoms
  //----------------------------------------
  
  // Define axis of rotation
  Vector axis = { p30->N->p.x - p30->CA->p.x, p30->N->p.y - p30->CA->p.y, p30->N->p.z - p30->CA->p.z };
  vector_normalize(&axis);

  // Update atom coordinates for residue 30
  // Skip atom 0 (the backbone nitrogen) and atom 1 (the alpha carbon)...start with atom 2 (the beta carbon)
  // Also skip the nitrogen's hydrogen
  for(i = 2; i < amino_acid_num_atoms(p30); i++)
  {
    Atom *a = (Atom *)array_get(p30->atoms, i);
    if(strcmp(a->atom_name, "H") != 0)
    {
      Vector v = { a->p.x - p30->CA->p.x, a->p.y - p30->CA->p.y, a->p.z - p30->CA->p.z };
      Vector *v_rot = vector_rotate(&axis, &v, 0-phi);
      a->p.x = p30->CA->p.x + v_rot->x;
      a->p.y = p30->CA->p.y + v_rot->y;
      a->p.z = p30->CA->p.z + v_rot->z;
      vector_delete(v_rot);
    }
  }
  
  // Update all atom coordinates for subsequent residues
  for(i = 30; i < array_size(polypeptide); i++)
  {
    AminoAcid *aa = (AminoAcid *)array_get(polypeptide, i);
    int j;
    for(j = 0; j < amino_acid_num_atoms(aa); j++)
    {
      Atom *a = (Atom *)array_get(aa->atoms, j);
      Vector v = { a->p.x - p30->CA->p.x, a->p.y - p30->CA->p.y, a->p.z - p30->CA->p.z };
      Vector *v_rot = vector_rotate(&axis, &v, 0-phi);
      a->p.x = p30->CA->p.x + v_rot->x;
      a->p.y = p30->CA->p.y + v_rot->y;
      a->p.z = p30->CA->p.z + v_rot->z;
      vector_delete(v_rot);
    }
  }
  
  // Now set psi to 0, rotate all affected atoms
  //----------------------------------------
  
  // Define new axis of rotation
  axis = (Vector){ p30->CA->p.x - p30->C->p.x, p30->CA->p.y - p30->C->p.y, p30->CA->p.z - p30->C->p.z };
  vector_normalize(&axis);
  
  // Update O coordinates for residue 30
  for(i = 0; i < amino_acid_num_atoms(p30); i++)
  {
    Atom *a = (Atom *)array_get(p30->atoms, i);
    if(strcmp(a->atom_name, "O") == 0)
    {
      Vector v = { a->p.x - p30->C->p.x, a->p.y - p30->C->p.y, a->p.z - p30->C->p.z };
      Vector *v_rot = vector_rotate(&axis, &v, 0-psi);
      a->p.x = p30->C->p.x + v_rot->x;
      a->p.y = p30->C->p.y + v_rot->y;
      a->p.z = p30->C->p.z + v_rot->z;
      vector_delete(v_rot);
      break;
    }
  }
  
  // Update all atom coordinates for subsequent residues
  for(i = 30; i < array_size(polypeptide); i++)
  {
    AminoAcid *aa = (AminoAcid *)array_get(polypeptide, i);
    int j;
    for(j = 0; j < amino_acid_num_atoms(aa); j++)
    {
      Atom *a = (Atom *)array_get(aa->atoms, j);
      Vector v = { a->p.x - p30->C->p.x, a->p.y - p30->C->p.y, a->p.z - p30->C->p.z };
      Vector *v_rot = vector_rotate(&axis, &v, 0-psi);
      a->p.x = p30->C->p.x + v_rot->x;
      a->p.y = p30->C->p.y + v_rot->y;
      a->p.z = p30->C->p.z + v_rot->z;
      vector_delete(v_rot);
    }
  }
  
  // Print out new coordinates
  instream = fopen(argv[1], "r");
  FILE *outstream = fopen(argv[2], "w");
  if(outstream == NULL)
  {
    fprintf(stderr, "Cannot open output file '%s'\n", argv[2]);
    return 1;
  }
  while(fgets(buffer, BUFFER_SIZE, instream))
  {
    if(strlen(buffer) > 4 && strncmp(buffer, "ATOM", 4) == 0)
      break;
    fputs(buffer, outstream);
  }
  for(i = 0; i < array_size(polypeptide); i++)
  {
    AminoAcid *aa = (AminoAcid *)array_get(polypeptide, i);
    int j;
    for(j = 0; j < amino_acid_num_atoms(aa); j++)
    {
      Atom *a = (Atom *)array_get(aa->atoms, j);
      atom_print(a, outstream);
    }
  }
  while(fgets(buffer, BUFFER_SIZE, instream))
  {
    if(strlen(buffer) < 4 || strncmp(buffer, "ATOM", 4) != 0)
      fputs(buffer, outstream);
  }
  fclose(instream);
  fclose(outstream);
  
  // Check for steric clashing
  float min_distance = 0.0f;
  Atom *a1 = NULL;
  Atom *a2 = NULL;
  for(i = 0; i < array_size(polypeptide); i++)
  {
    AminoAcid *p_i = (AminoAcid *)array_get(polypeptide, i);
    int j;
    for(j = i+1; j < array_size(polypeptide); j++)
    {
      AminoAcid *p_j = (AminoAcid *)array_get(polypeptide, j);            
      int k;
      for(k = 0; k < amino_acid_num_atoms(p_i); k++)
      {
        Atom *a_temp_i = (Atom *)array_get(p_i->atoms, k);
        int l;
        for(l = 0; l < amino_acid_num_atoms(p_j); l++)
        {
          Atom *a_temp_j = (Atom *)array_get(p_j->atoms, l);
          float test_distance = euclidean_distance(&a_temp_i->p, &a_temp_j->p);
          if(a1 == NULL || test_distance < min_distance)
          {
            min_distance = test_distance;
            a1 = a_temp_i;
            a2 = a_temp_j;
          }
//printf("~ test_distance=%.5f, min_distance=%.5f\n", test_distance, min_distance);
        }
      }
    }
  }
  puts("Steric clashing");
  printf("    Minimum atom distance: %.5fÅ (between %s of %s %d and %s of %s %d)\n\n", min_distance, a1->atom_name, a1->residue_name, a1->residue_sequence_number, a2->atom_name, a2->residue_name, a2->residue_sequence_number);

  puts("Extra credit--side chain torsion angles:");
  for(i = 0; i < array_size(polypeptide); i++)
  {
    AminoAcid *aa = (AminoAcid *)array_get(polypeptide, i);
    printf("    Residue %d (%s):\n", aa->N->residue_sequence_number, aa->N->residue_name);
    amino_acid_print_chis(aa, stdout);
  }

  // Free memory
  free(bond_lengths_N_CA);
  free(bond_lengths_CA_C);
  free(bond_lengths_C_N);
  free(bond_angles_N_CA_C);
  free(bond_angles_CA_C_N);
  free(bond_angles_C_N_CA);
  free(carbon_distances);
  for(i = 0; i < array_size(polypeptide); i++)
  {
    AminoAcid *aa = (AminoAcid *)array_get(polypeptide, i);
    amino_acid_delete(aa);
  }
  array_delete(polypeptide);

  return 0;
}
