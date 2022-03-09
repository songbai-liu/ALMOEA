//  IMF8.java
//
//  Author:
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.problems.IMF;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem DTLZ1
 */
public class IMF8 extends Problem {
	/**
	 * Creates a default DTLZ1 problem (7 variables and 3 objectives)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public IMF8(String solutionType) throws ClassNotFoundException {
		this(solutionType, 30, 3);
	} // DTLZ1

	/**
	 * Creates a DTLZ1 problem instance
	 * 
	 * @param numberOfVariables
	 *            Number of  
	 * @param numberOfObjectives
	 *            Number of objective functions
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public IMF8(String solutionType, Integer numberOfVariables,
			Integer numberOfObjectives) {
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = 3;
		numberOfConstraints_ = 0;
		problemName_ = "IMF8";

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

	/**
	 * Evaluates a solution
	 * 
	 * @param solution
	 *            The solution to evaluate
	 * @throws JMException
	 */
	public void evaluate(Solution solution) throws JMException {
		Variable[] decisionVariables = solution.getDecisionVariables();

		double[] x = new double[numberOfVariables_];
		double[] f = new double[3];

		for (int i = 0; i < numberOfVariables_; i++)
			x[i] = decisionVariables[i].getValue();

		double g = 0.0;
		double beta = 3.0;
		for (int i = 2; i < numberOfVariables_; i++)
			g += Math.pow(Math.pow(x[i], 1.0/(1.0+beta*(i+1)/numberOfVariables_))-x[0], 2);
		
		f[0] = (1+g)*Math.cos(0.5*x[0]*Math.PI)*Math.cos(0.5*x[1]*Math.PI);

	    f[1] = (1+g)*Math.cos(0.5*x[0]*Math.PI)*Math.sin(0.5*x[1]*Math.PI);
	    
	    f[2] = (1+g)*Math.sin(0.5*x[0]*Math.PI);

		for (int i = 0; i < numberOfObjectives_; i++)
			solution.setObjective(i, f[i]);
	} // evaluate

}
