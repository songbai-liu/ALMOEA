package jmetal.util.comparators;

import jmetal.core.Solution;

import java.util.Comparator;

public class ObjectiveDistanceCompartor  implements Comparator{   //将目标距离和从小到大进行排序
	 /**
	  * Compares two solutions.
	  * @param o1 Object representing the first <code>Solution</code>.
	  * @param o2 Object representing the second <code>Solution</code>.
	  * @return -1, or 0, or 1 if o1 is less than, equal, or greater than o2,
	  * respectively.
	  */
	  public int compare(Object o1, Object o2) {
	    if (o1==null)
	      return 1;
	    else if (o2 == null)
	      return -1;
	    
	    double distance1 = ((Solution)o1).getCrowdingDistance();
	    double distance2 = ((Solution)o2).getCrowdingDistance();
	    if (distance1 >  distance2)
	      return 1;
	    if (distance1 < distance2)
	      return - 1;
	    return 0;    
	  } // compare
}
