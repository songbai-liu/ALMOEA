package jmetal.util;

public class PairData implements Comparable<PairData> {
	private int dataId;
	private double dataValue;

	public PairData(int id, double value) {
		this.dataId = id;
		this.dataValue = value;
	}

	public int getID(){
		return dataId;
	}
	
	public double getValue(){
		return dataValue;
	}
	

	public int compareTo(PairData sObj) {

		if (dataValue < sObj.dataValue)
			return -1;
		else if (dataValue > sObj.dataValue)
			return 1;
		else
			return 0;
	}
}
	