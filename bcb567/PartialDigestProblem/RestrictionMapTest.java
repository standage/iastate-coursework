/**
 * Procedural solution to Partial Digest Problem
 */
import java.util.Iterator;
import java.util.Set;
import java.util.HashSet;
import mathCollection.*;

public class RestrictionMapTest
{
  private static Multiset L;
  private static Set X;

  public static void main(String [] args)
  {
    L = new HashMultiset();
    L.add(new Integer(2));
    L.add(new Integer(2));
    L.add(new Integer(3));
    L.add(new Integer(3));
    L.add(new Integer(4));
    L.add(new Integer(5));
    L.add(new Integer(6));
    L.add(new Integer(7));
    L.add(new Integer(8));
    L.add(new Integer(10));
    
    PartialDigest(L);
  }

  private static void PartialDigest(Multiset L)
  {
    Integer width = MaximumValue(L);
    L.remove(width);
    X = new HashSet();
    X.add(new Integer(0));
    X.add(width);
    Place(L, X, width);
  }
  
  private static void Place(Multiset L, Set X, Integer width)
  {
    if(L.size() == 0)
    {
      System.out.println("Success: "+ X);
      return;
    }
    
    Integer y = MaximumValue(L);
    Multiset d = delta(y, X);
    if(L.isSuperset(d))
    {
      X.add(y);
      Iterator<Integer> i = d.iterator();
      while(i.hasNext())
      {
        Integer toRemove = i.next();
        L.remove(toRemove);
      }
      Place(L, X, width);
      X.remove(y);
      i = d.iterator();
      while(i.hasNext())
      {
        Integer toAdd = i.next();
        L.add(toAdd);
      }
    }
    
    Integer yprime = new Integer(width.intValue() - y.intValue());
    d = delta(yprime, X);
    if(L.isSuperset(d))
    {
      X.add(yprime);
      Iterator<Integer> i = d.iterator();
      while(i.hasNext())
      {
        Integer toRemove = i.next();
        L.remove(toRemove);
      }
      Place(L, X, width);
      X.remove(yprime);
      i = d.iterator();
      while(i.hasNext())
      {
        Integer toAdd = i.next();
        L.add(toAdd);
      }
    }
  }
  
  private static int MaximumValue(Multiset set)
  {
    Integer maxvalue = null;
    Iterator<Integer> i = set.iterator();
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
  
  private static Multiset delta(Integer y, Set X)
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
}
