//  CEC2009_UF6.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.problems.cec2009Competition;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem CEC2009_UF5
 */
public class UF6 extends Problem {
	int N_;
	double epsilon_;

	/**
	 * Constructor. Creates a default instance of problem CEC2009_UF6 (30
	 * decision variables)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public UF6(String solutionType) throws ClassNotFoundException {
		this(solutionType, 30, 2, 0.1); // 30 variables, N =10, epsilon = 0.1
	} // CEC2009_UF1

	/**
	 * Creates a new instance of problem CEC2009_UF6.
	 * 
	 * @param numberOfVariables
	 *            Number of variables.
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public UF6(String solutionType, Integer numberOfVariables, int N,
			double epsilon) {
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = 2;
		numberOfConstraints_ = 0;
		problemName_ = "UF6";

		N_ = N;
		epsilon_ = epsilon;

		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		lowerLimit_[0] = 0.0;
		upperLimit_[0] = 1.0;
		for (int var = 1; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -1.0;
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
	} // CEC2009_UF6

	/**
	 * Evaluates a solution.
	 * 
	 * @param solution
	 *            The solution to evaluate.
	 * @throws JMException
	 */
	public void evaluate(Solution solution) throws JMException {
		Variable[] decisionVariables = solution.getDecisionVariables();

		double[] x = new double[numberOfVariables_];
		for (int i = 0; i < numberOfVariables_; i++)
			x[i] = decisionVariables[i].getValue();

		int count1, count2;
		double prod1, prod2;
		double sum1, sum2, yj, hj, pj;
		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		prod1 = prod2 = 1.0;

		for (int j = 2; j <= numberOfVariables_; j++) {
			yj = x[j - 1]
					- Math.sin(6.0 * Math.PI * x[0] + j * Math.PI
							/ numberOfVariables_);
			pj = Math.cos(20.0 * yj * Math.PI / Math.sqrt(j));
			if (j % 2 == 0) {
				sum2 += yj * yj;
				prod2 *= pj;
				count2++;
			} else {
				sum1 += yj * yj;
				prod1 *= pj;
				count1++;
			}
		}
		hj = 2.0 * (0.5 / N_ + epsilon_) * Math.sin(2.0 * N_ * Math.PI * x[0]);
		if (hj < 0.0)
			hj = 0.0;

		solution.setObjective(0, x[0] + hj + 2.0
				* (4.0 * sum1 - 2.0 * prod1 + 2.0) / (double) count1);
		solution.setObjective(1, 1.0 - x[0] + hj + 2.0
				* (4.0 * sum2 - 2.0 * prod2 + 2.0) / (double) count2);
	} // evaluate
} // CEC2009_UF6
