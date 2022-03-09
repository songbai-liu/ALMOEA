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

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Differential evolution crossover operators Comments: - The operator receives
 * two parameters: the current individual and an array of three parent
 */
public class DeepDECrossover extends Crossover {
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
	public DeepDECrossover(HashMap<String, Object> parameters) {
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
		Solution[][] parent = (Solution[][]) parameters[1];
		int[][] index = (int[][]) parameters[2];
		Solution child;
		child = new Solution(current);
		XReal xChild = new XReal(child);
        int layer = index.length;
        for(int i=0;i<layer;i++){
        	int jrand;
    		XReal xParent1 = new XReal(parent[i][0]);
    		XReal xParent2 = new XReal(parent[i][1]);
    		int numberOfVariables = index[i].length;
    		//double r1 = PseudoRandom.randDouble(0, 1);
    		//double r2 = PseudoRandom.randDouble(0, 1);
    		// STEP 4. DE Operator
    		for (int j = 0; j < numberOfVariables; j++) {
    			int d = index[i][j];
    			double value;
				value=PseudoRandom.randDouble(0, 1);
				value = xChild.getValue(d) + F_* (xParent1.getValue(d) - xParent2.getValue(d));
				if (value < xChild.getLowerBound(d))
					value = xChild.getLowerBound(d);
				if (value > xChild.getUpperBound(d))
					value = xChild.getUpperBound(d);
				xChild.setValue(d, value);
    		} // for
        }
		return child;
	}
} // DifferentialEvolutionCrossover
