package jmetal.problems.LSMOP;
import jmetal.core.Problem;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;

import java.util.Random;

public abstract class LSMOP extends Problem {
	/**
	 * stores a epsilon default value
	 */
	protected int nk; //Number of subcomponents in each variable group
	protected double[] c;
	protected int[] sublen;//Number of Variables in each subcomponent
	protected int[] len;//Cumulative sum of lengths of variable groups
	protected Random random = new Random();
	
	/**
	 * Creates a LSMOP problem instance
	 * 
	 * @param numberOfVariables
	 *            Number of  
	 * @param numberOfObjectives
	 *            Number of objective functions
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public LSMOP(String solutionType, Integer numberOfVariables,
			Integer numberOfObjectives) {
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = numberOfObjectives;
		numberOfConstraints_ = 0;

		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		for (int var = 0; var < numberOfVariables; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for

		if (solutionType.compareTo("BinaryReal") == 0)
			solutionType_ = new BinaryRealSolutionType(this);
		else if (solutionType.compareTo("Real") == 0)
			solutionType_ = new RealSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType
					+ " invalid");
			System.exit(-1);
		}
	}
}
