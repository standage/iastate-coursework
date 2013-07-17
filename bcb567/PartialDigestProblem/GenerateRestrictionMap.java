import java.util.Set;
import mathCollection.*;

public class GenerateRestrictionMap
{
  public static void main(String [] args)
  {
    Multiset L = new HashMultiset();
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
    
    RestrictionMapper mapper = new RestrictionMapper();
    Set<Set> maps = mapper.generateRestrictionMaps(L);
    System.out.println(maps);
  }
}
