#include "Atom.h"

double euclidean_distance(Point *p, Point *q)
{
  double sum_squared_distances = 0.0;

  sum_squared_distances += pow(q->x - p->x, 2);
  sum_squared_distances += pow(q->y - p->y, 2);
  sum_squared_distances += pow(q->z - p->z, 2);

  return sqrt(sum_squared_distances);
}

void atom_delete(Atom *a)
{
  int i;
  for(i = 0; i < array_size(a->surface); i++)
  {
    Point *p = (Point *)array_get(a->surface, i);
    free(p);
  }
  array_delete(a->surface);
  array_delete(a->accessible);
  free(a);
  a = NULL;
}

void atom_sample_surface(Atom *a, int n)
{
  int i;
  for(i = 0; i < n; i++)
  {
    Point *sampled_point = (Point *)malloc( sizeof(Point) );
    
    // Generate a random z value (representing a longitude)
    double z_range = a->radius * 2;
    double z_offset = rand() / (((double)RAND_MAX + 1) / z_range);
    sampled_point->z =  z_offset - a->radius;
    
    // Generate a random phi value (representing a latitude)
    // Calculate x and y coordinates from z and phi coordinates
    double phi = rand() / (((double)RAND_MAX + 1) / (2 * 3.14159265));
    double theta = asin(sampled_point->z / a->radius);
    sampled_point->x = a->radius * cos(theta) * cos(phi);
    sampled_point->y = a->radius * cos(theta) * sin(phi);

    // The point is centered at the origin...make sure to transform to center it
    // on the atomic center
    sampled_point->x += a->p.x;
    sampled_point->y += a->p.y;
    sampled_point->z += a->p.z;
 
    array_add(a->surface, sampled_point);
  }
}

Atom *atom_new(Point *p, char element, int rsn)
{
  Atom *a = (Atom *)malloc( sizeof(Atom) );

  a->p.x = p->x;
  a->p.y = p->y;
  a->p.z = p->z;
  a->element = element;
  a->residue_sequence_number = rsn;
  a->sasa = 0.0;

  switch(element)
  {
    case 'C':
      a->radius = VDWR_CARBON + PROBE_RADIUS;
      break;
    case 'H':
      a->radius = VDWR_HYDROGEN + PROBE_RADIUS;
      break;
    case 'N':
      a->radius = VDWR_NITROGEN + PROBE_RADIUS;
      break;
    case 'O':
      a->radius = VDWR_OXYGEN + PROBE_RADIUS;
      break;
    case 'S':
      a->radius = VDWR_SULFUR + PROBE_RADIUS;
      break;
    default:
      fprintf(stderr, "Error: unknown element '%c'\n", element);
      exit(1);
      break;
  }
  
  a->surface = array_new( sizeof(Point *) );
  a->accessible = array_new( sizeof(Point *) );
  atom_sample_surface(a, POINTS_PER_ATOM);
  
  return a;
}

Atom *parse_atom_from_pdb(const char *line)
{
  char buffer[16] = "";
  Point p;

  // Parse x,y,z coordinates
  strncpy(buffer, line + 30, 8);
  p.x = atof(buffer);
  strncpy(buffer, line + 38, 8);
  p.y = atof(buffer);
  strncpy(buffer, line + 46, 8);
  p.z = atof(buffer);

  // Parse residue sequence number
  int i = 22;
  while(line[i] == ' ')
    i++;
  int l = 25 - i + 1;
  strncpy(buffer, line + i, l);
  buffer[l] = '\0';
  int residue_sequence_number = atoi(buffer);

  // Parse element
  char element = line[77];
  
  return atom_new(&p, element, residue_sequence_number);
}
