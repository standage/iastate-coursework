/**
 * @author Daniel S. Standage
 *
 * This program implements a greedy algorithm for approximating the
 * minimum number of genome rearrangements between two related genomes.
 * A formalization of this problem is given by Jones and Pevzner in
 * 'An Introduction to Bioinformatics Algorithms', pp. 125-135.
 */

#include <stdio.h>

typedef int SyntenyBlock;
typedef int Direction;

/**
 * Data structure representing the endpoints (inclusive)
 * of a reversal permutation
 */
typedef struct
{
  int minIndex;
  int maxIndex;
} ReversalPermutation;

/**
 * Data structure representing the endpoints of a 
 * strip of ordered syteny blocks
 */
typedef struct
{
  int minIndex;
  int maxIndex;
  Direction direction;
} Strip;

/**
 * Determine the number of breakpoints in the given permutation
 *
 * @param int * pi      Pointer to a list of integers representing a permutation of synteny blocks
 * @param int length    The number of integers in the permutation (including "bogus" values on ends)
 * @returns             The number of breakpoints in the given permutation
 */
int b(SyntenyBlock * pi, int length, Strip * strips)
/*{
  int i;
  int breakpoints = 0;
  for(i = 0; i < length - 1; i++) 
  {
    if(pi[i] != pi[i+1] + 1 || pi[i] != pi[i+1] - 1)
    {
      breakpoints++;
    }
  } 
  return breakpoints;
}
int getStrips(SyntenyBlock * pi, int length, Strip * strips) */
{
  int numStrips = 0;
  int i;
  int stripStart = 0;
  for(i = 0; i < length - 1; i++)
  {
    if(pi[i] != pi[i+1] + 1 || pi[i] != pi[i+1] - 1)
    {
      Direction direction;
      if(stripStart == i || pi[i] == pi[i-1]-1)
      {
        direction = -1;
      }
      else
      {
        direction = 1;
      }
      strips[numStrips] = { stripStart, i, direction };
      numStrips++;
    }
  }
  return numStrips - 1;
}

/**
 *
 */
ReversalPermutation chooseRho(SyntenyBlock * pi, int length, Strip * strips)
{
  SyntenyBlock k = 0;
  int i;

  // Look for smallest element k in a decreasing strip
  for(i = 0; i < numStrips; i++)
  {
    Strip s = strips[i];
    if(s.direction == -1 && s.maxIndex < k)
    {
      k = s.maxIndex;
    }
  }

  // Determine the order of the indices of the k^th and
  // the (k-1)^th synteny blocks, set rho accordingly
  int indexOfk = indexOf(k);
  int indexOfkMinus1 = indexOf(k - 1);
  ReversalPermutation rho;
  if(indexOfk > indexOfkMinus1)
  {
    rho = { indexOfkMinus1, indexOfk };
  }
  else
  {
    rho = { indexOfk + 1, indexOfkMinus1 };
  }
  return rho;
}

/**
 * Determine whether there are any remaining decreasing strips
 *
 * @param   Strip * strips     Pointer to a list of strips
 * @param   int length         The maximum number of strips
 * @returns int                Whether there are any decreasing strips
 */
int hasDecreasingStrip(Strip * strips, int length)
{
  int i;
  for(i = 0; i < length; i++)
  {
    Strip s = strips[i];
    if(s.direction == -1)
    {
      return 1;
    }
  }
  return 0;
}

/**
 * Determine the index of the given synteny block
 *
 * @param   SyntenyBlock block     The synteny block for which we want the index
 * @returns int                    The index of the synteny block, or -1 if the block is not found
 */
int indexOf(SyntenyBlock * pi, int length, SyntenyBlock block)
{
  int i;
  for(i = 0; i < length; i++)
  {
    if(pi[i] == block)
    {
      return i;
    }
  }
  return -1;
}

/**
 * Perform in-place reversal of sublist of synteny blocks in the given rearrangement (permutation)
 *
 * @param SyntenyBlock * pi     Pointer to a list of synteny blocks that represents a permutation
 * @param int minIndex          Index of first integer to include in reversal
 * @param int maxIndex          Index of last integer to include in reversal
 */
void permute(int * pi, ReversalPermutation rho)
{
  int temp;
  while(rho.minIndex < rho.maxIndex)
  {
    temp = pi[rho.minIndex];
    pi[rho.minIndex] = pi[rho.maxIndex];
    pi[rho.maxIndex] = temp;
    rho.minIndex++;
    rho.maxIndex--;
  }
}

int main()
{
  int i;
  int pi[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  int length = 10;
  Strip strips[length];

  ReversalPermutation rho = {2, 5};
  permute(pi, rho);
  for(i = 0; i < 10; i++)
  {
    printf("i: %d\n", pi[i]);
  }
  printf("min: %d, max: %d\n", rho.minIndex, rho.maxIndex);

  return 0;
}
