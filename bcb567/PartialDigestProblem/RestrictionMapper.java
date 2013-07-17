/**
 * @author Daniel S. Standage
 *
 * This class implements a branch-and-bound algorithm for generating
 * a restriction map from a partial restriction enzyme digest. A
 * formalization of this problem is given by Jones and Pevzner in
 * 'An Introduction to Bioinformatics Algorithms', pp. 83-91.
 */

/**
 * Required libraries
 *
 * The 'mathCollection.jar' package can be obtained from
 * {@link http://mathcollection.berlios.de/}.
 * The documentation for this package is available at
 * {@link http://mathcollection.berlios.de/doc/mathCollection/package-summary.html}.
 */
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import mathCollection.*;

/**
 * Class implementation
 */
public class RestrictionMapper
{
 /**
  * Instance data
  *
  * Set<Set> solutions     A set of restriction maps that generate the specified
  *                        partial digest
  */
   private Set<Set> solutions;

 /**
  * Constructor
  */
  public RestrictionMapper()
  {
    this.solutions = new HashSet();
  }

  /**
   * Given a multiset of restriction fragment lengths L, determine
   * the set(s) of restriction sites X = {x_1, x_2, ...} that
   * generate the partial digest represented by L
   *
   * @param   Multiset partialDigest     A multiset of integers L, representing the
   *                                     lengths of fragments from a partial
   *                                     restriction enzyme digest
   * @returns Set<Set>                   A set X = {x_1, x_2, ...) where x_i is a
   *                                     set of integers representing the positions
   *                                     of restriction enzyme cut sites
   */
  public Set<Set> generateRestrictionMaps(Multiset partialDigest)
  {
    Integer width = this.getMaximumValue(partialDigest);
    partialDigest.remove(width);
    Set<Integer> restrictionMap = new HashSet();
    restrictionMap.add(new Integer(0));
    restrictionMap.add(width);
    this.place(partialDigest, restrictionMap, width);
    partialDigest.add(width);
    return solutions;
  }

  /**
   * @param   Integer      y     An integer
   * @param   Set<Integer> X     A set of integers
   * @returns Multiset           The multiset of distances between each element x in X and y
   */
  private Multiset delta(Integer y, Set<Integer> X)
  {
    Multiset d = new HashMultiset();
    Iterator<Integer> i = X.iterator();
    while(i.hasNext())
    {
      Integer x = i.next();
      int diff = java.lang.Math.abs(x.intValue() - y.intValue());
      d.add(new Integer(diff));
    }
    return d;
  }

  /**
   * @param   Multiset s     A multiset of integers
   * @returns Integer        The maximum integer value in set s
   */
  private Integer getMaximumValue(Multiset s)
  {
    Integer maxvalue = null;
    Iterator<Integer> i = s.iterator();
    while(i.hasNext())
    {
      Integer testvalue = i.next();
      if(maxvalue == null || testvalue.compareTo(maxvalue) > 0)
      {
        maxvalue = testvalue;
      }
    }
    return maxvalue;
  }

  /**
   * Recursive place method...see description in Jones and Pevzner
   */
  private void place(Multiset partialDigest, Set<Integer> restrictionMap, Integer width)
  {
    // Base case...bottom out
    if(partialDigest.size() == 0)
    {
      this.storeSolution(restrictionMap);
      return;
    }

    // Recurse left
    Integer y = this.getMaximumValue(partialDigest);
    Multiset delta = this.delta(y, restrictionMap);
    if(partialDigest.isSuperset(delta))
    {
      restrictionMap.add(y);
      Iterator<Integer> i = delta.iterator();
      while(i.hasNext())
      {
        Integer toRemove = i.next();
        partialDigest.remove(toRemove);
      }
      this.place(partialDigest, restrictionMap, width);
      restrictionMap.remove(y);
      i = delta.iterator();
      while(i.hasNext())
      {
        Integer toAdd = i.next();
        partialDigest.add(toAdd);
      }
    }

    // Recurse right
    Integer yprime = new Integer(width.intValue() - y.intValue());
    delta = this.delta(yprime, restrictionMap);
    if(partialDigest.isSuperset(delta))
    {
      restrictionMap.add(yprime);
      Iterator<Integer> i = delta.iterator();
      while(i.hasNext())
      {
        Integer toRemove = i.next();
        partialDigest.remove(toRemove);
      }
      this.place(partialDigest, restrictionMap, width);
      restrictionMap.remove(yprime);
      i = delta.iterator();
      while(i.hasNext())
      {
        Integer toAdd = i.next();
        partialDigest.add(toAdd);
      }
    }
  }

  /**
   * Store a restriction map that generates the specified partial digest
   * Since Java passes by pointer, deep copy is necessary
   */
  private void storeSolution(Set<Integer> restrictionMap)
  {
    Set<Integer> solution = new HashSet();
    Iterator<Integer> i = restrictionMap.iterator();
    while(i.hasNext())
    {
      Integer position = i.next();
      solution.add(new Integer(position.intValue()));
    }
    this.solutions.add(solution);
  }
}