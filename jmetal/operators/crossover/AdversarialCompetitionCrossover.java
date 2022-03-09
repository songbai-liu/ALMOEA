//  DifferentialEvolutionCrossover.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.operators.crossover;

import jmetal.core.Solution;
import jmetal.encodings.solutionType.ArrayRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.wrapper.XReal;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Differential evolution crossover operators Comments: - The operator receives
 * two parameters: the current individual and an array of three parent
 */
public class AdversarialCompetitionCrossover extends Crossover {
	/**
	 * DEFAULT_CR defines a default CR (crossover operation control) value
	 */
	private static final double DEFAULT_CR = 0.5;

	/**
	 * DEFAULT_F defines the default F (Scaling factor for mutation) value
	 */
	private static final double DEFAULT_F = 0.5;

	/**
	 * Valid solution types to apply this operator
	 */
	private static final List VALID_TYPES = Arrays.asList(
			RealSolutionType.class, ArrayRealSolutionType.class);

	private double CR_;
	private double F_;

	/**
	 * Constructor
	 */
	public AdversarialCompetitionCrossover(HashMap<String, Object> parameters) {
		super(parameters);

		CR_ = DEFAULT_CR;
		F_ = DEFAULT_F;

		if (parameters.get("CR") != null)
			CR_ = (Double) parameters.get("CR");
		if (parameters.get("F") != null)
			F_ = (Double) parameters.get("F");
	} // Constructor

	/**
	 * Executes the operation
	 * 
	 * @param object
	 *            An object containing an array of three parents
	 * @return An object containing the offSprings
	 */
	public Object execute(Object object) throws JMException {
		Object[] parameters = (Object[]) object;
		Solution current = (Solution) parameters[0];
		Solution[] parent = (Solution[]) parameters[1];
		double[] contribution = (double[]) parameters[2];
		double[] ub = (double[]) parameters[3];
		double[] lb = (double[]) parameters[4];
		int type = (int) parameters[5];
		Solution child;
		child = new Solution(current);
		XReal xChild = new XReal(child);
		int jrand;
		XReal xParent1 = new XReal(parent[0]);
		XReal xParent2 = new XReal(parent[1]);
		int numberOfVariables = xChild.getNumberOfDecisionVariables();
		List<Integer> index = new ArrayList<Integer>();
		for (int j = 0; j < numberOfVariables; j++) {
			index.add(j, j);
		}
		//double r1 = PseudoRandom.randDouble(0, 1);
		//double r2 = PseudoRandom.randDouble(0, 1);
		// STEP 4. DE Operator
		int rd1, rd2, rd;
		double value;
		if(type == 0 || type == 1){
			while(numberOfVariables > 0){
				rd1 = PseudoRandom.randInt(0, numberOfVariables-1);
				rd2 = PseudoRandom.randInt(0, numberOfVariables-1);
				while(rd2 == rd1){
					rd2 = PseudoRandom.randInt(0, numberOfVariables-1);
				}
				if(contribution[index.get(rd1)] > contribution[index.get(rd2)]){
					rd = index.get(rd1);
				}else{
					rd = index.get(rd2);
				}
				value = xChild.getValue(rd) + F_* (xParent1.getValue(rd) - xParent2.getValue(rd));
				
				if (value < lb[rd])
					value = lb[rd];
				if (value > ub[rd])
					value = ub[rd];
			
				xChild.setValue(rd, value);
				
				if(rd1 > rd2){
					index.remove(rd1);
					index.remove(rd2);
				}else{
					index.remove(rd2);
					index.remove(rd1);
				}
				numberOfVariables -= 2;
			}
		}else{
			while(numberOfVariables > 0){
				rd = numberOfVariables - 1;
				value = xChild.getValue(rd) + F_* (xParent1.getValue(rd) - xParent2.getValue(rd));
				
				if (value < xChild.getLowerBound(rd))
					value = xChild.getLowerBound(rd);
				if (value > xChild.getUpperBound(rd))
					value = xChild.getUpperBound(rd);
			
				xChild.setValue(rd, value);
				
				numberOfVariables -= 1;
			}
		}
		return child;
	}//execute
} // DifferentialEvolutionCrossover
