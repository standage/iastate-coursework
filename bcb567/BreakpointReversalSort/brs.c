/**
 * @author Daniel S. Standage
 *
 * This program implements a greedy algorithm for approximating the
 * minimum number of genome rearrangements between two related genomes.
 * This algorithm is a relatively naive approach and ignores a lot of
 * biological reality associated with genome rearrangements. It is
 * used rather as a tool to teach the application of greedy algorithms
 * to problems in biology and bioinformatics. A formalization of this
 * problem is given by Jones and Pevzner in 'An Introduction to
 * Bioinformatics Algorithms', pp. 125-135.
 */

#include <stdio.h>

typedef int SyntenyBlock;

/**
 * Data structure representing the endpoints (inclusive)
 * of a reversal permutation
 */
typedef struct _ReversalPermutation_
{
  int minIndex;
  int maxIndex;
} ReversalPermutation;

/**
 * Data structure representing the endpoints of a 
 * strip of ordered syteny blocks
 */
typedef struct _Strip_
{
  int minIndex;
  int maxIndex;
  int direction;
} Strip;

/**
 * Function prototypes
 */
int b(SyntenyBlock * pi, int length, Strip * strips);
ReversalPermutation chooseRho(SyntenyBlock * pi, int length, Strip * strips, int numStrips);
int hasDecreasingStrip(Strip * strips, int numStrips);
int indexOf(SyntenyBlock * pi, int length, SyntenyBlock block);
void permute(int * pi, ReversalPermutation rho);

/**
 * Determine the number of breakpoints in the given permutation
 *
 * @param SyntenyBlock * pi     Pointer to a list of integers representing a permutation of synteny blocks
 * @param int length            The number of integers in the permutation (including "bogus" values on ends)
 * @param Strip * strips        Pointer to a list of strips in the permutation
 * @returns                     The number of breakpoints in the given permutation
 */
int b(SyntenyBlock * pi, int length, Strip * strips)
{
  int numStrips = 0;
  int i;
  int stripStart = 0;
  int direction;
  for(i = 0; i < length - 1; i++)
  {
    if(pi[i] != pi[i+1] + 1 && pi[i] != pi[i+1] - 1)
    {
      if((stripStart == i && stripStart != 0 && stripStart != length - 1) || pi[i] == pi[i-1]-1)
      {
        direction = -1;
      }
      else
      {
        direction = 1;
      }
      Strip s = { stripStart, i, direction };
      strips[numStrips] = s;
      numStrips++;
      stripStart = i+1;
    }
  }
  if(stripStart == length - 1 || pi[length-1] == pi[length-1] - 1)
  {
    direction = -1;
  }
  else
  {
    direction = 1;
  }
  Strip s = { stripStart, length - 1, direction };
  strips[numStrips] = s;
  numStrips++;
  return numStrips - 1;
}

/**
 * Choose the next reversal permutation
 *
 * @param   SyntenyBlock * pi       Pointer to a list of synteny block representing a permutation
 * @param   int length              Number of values in pi
 * @param   Strip * strips          Pointer to a list of strips
 * @param   int numStrips           The number of strips in the list
 * @returns ReversalPermutation     An optimal permutation
 */
ReversalPermutation chooseRho(SyntenyBlock * pi, int length, Strip * strips, int numStrips)
{
  SyntenyBlock k = length;
  int i;

  for(i = 0; i < numStrips; i++)
  {
    Strip s = strips[i]; 
    if(s.direction == -1 && pi[s.maxIndex] < k)
    {
      k = pi[s.maxIndex];
    }
  }

  if(k == length)
  {
    Strip s = strips[numStrips - 2];
    ReversalPermutation rho = { s.minIndex, s.maxIndex };
    return rho;
  }

  int indexOfk = indexOf(pi, length, k);
  int indexOfkMinus1 = indexOf(pi, length, k - 1);
  if(indexOfk > indexOfkMinus1)
  {
    ReversalPermutation rho = { indexOfkMinus1 + 1, indexOfk };
    return rho;
  }
  ReversalPermutation rho = { indexOfk + 1, indexOfkMinus1 };
  return rho;
}

/**
 * Determine whether there are any remaining decreasing strips
 * This function is not necessary since the chooseRho function
 * returns the correct permutation rho whether or not there are
 * decreasing strips
 *
 * @param   Strip * strips     Pointer to a list of strips
 * @param   int length         The maximum number of strips
 * @returns int                Whether there are any decreasing strips
 */
int hasDecreasingStrip(Strip * strips, int numStrips)
{
  int i;
  for(i = 0; i < numStrips; i++)
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
 * @param   SyntenyBlock * pi      Pointer to a list of synteny blocks
 * @param   int length             The number of synteny blocks in the permutation
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
void permute(SyntenyBlock * pi, ReversalPermutation rho)
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

/**
 * Print pi
 */
void printPi(SyntenyBlock * pi, int length)
{
  int i;
  for(i = 1; i < length - 1; i++)
  {
    if(i > 1) printf("-");
    printf("%d", pi[i]);
  }
}

/**
 * Implementation of the ImprovedBreakpointReversalSort algorithm
 * listed in the text
 */
int main()
{
  //SyntenyBlock pi[] = { 0, 8, 2, 7, 6, 5, 1, 4, 3, 9 };
  //int length = 10;
  SyntenyBlock pi[] = { 0,19,18,9,21,23,12,11,13,4,10,7,3,25,20,6,1,24,17,5,14,8,15,2,16,22,26 };
  int length = 27;
  Strip strips[length];

  int breakpoints = b(pi, length, strips);
  while(breakpoints > 0)
  {
    // I do not check whether there is a decreasing strip here.
    // If there is not, the chooseRho function will return the
    // appropriate permutation to reverse an increasing strip
    ReversalPermutation rho = chooseRho(pi, length, strips, breakpoints + 1);

    printf("Permuting: reversing from %d to %d\n", pi[rho.minIndex], pi[rho.maxIndex]);
    fputs("\tBefore: ", stdout);
    printPi(pi, length);

    permute(pi, rho);

    fputs("\n\tAfter:  ", stdout);
    printPi(pi, length);
    puts("");

    breakpoints = b(pi, length, strips);
  }
  return 0;
}
