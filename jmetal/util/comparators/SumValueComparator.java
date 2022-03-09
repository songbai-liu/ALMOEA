package jmetal.util.comparators;

import java.util.Comparator;

import jmetal.core.Solution;

public class SumValueComparator implements Comparator{
	private boolean ascendingOrder_;


	public SumValueComparator() {
		ascendingOrder_ = true;
	} 

	public SumValueComparator(boolean descendingOrder) {
		if (descendingOrder)
			ascendingOrder_ = false;
		else
			ascendingOrder_ = true;
	} 

	public int compare(Object o1, Object o2) {
		if (o1 == null)
			return 1;
		else if (o2 == null)
			return -1;

		double f1 = ((Solution) o1).getSumValue();
		//double f1 = ((Solution) o1).getDistanceToIdealPoint();
		double f2 = ((Solution) o2).getSumValue();
		if (ascendingOrder_) {
			if (f1 < f2) {
				return -1;
			} else if (f1 > f2) {
				return 1;
			} else {
				return 0;
			}
		} else {
			if (f1 < f2) {
				return 1;
			} else if (f1 > f2) {
				return -1;
			} else {
				return 0;
			}
		}
	}
}
