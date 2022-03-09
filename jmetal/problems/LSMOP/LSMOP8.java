//  LSMOP8.java 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.problems.LSMOP;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem MaF1
 */
public class LSMOP8 extends LSMOP {
	/**
	 * Creates a default MaF1 problem (12 variables and 3 objectives)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public LSMOP8(String solutionType) throws ClassNotFoundException {
		this(solutionType, 300, 3);
	} // LSMOP8

	/**
	 * Creates a LSMOP8 problem instance
	 * 
	 * @param numberOfVariables
	 *            Number of  
	 * @param numberOfObjectives
	 *            Number of objective functions
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public LSMOP8(String solutionType, Integer numberOfVariables,
			Integer numberOfObjectives) {
		super(solutionType, numberOfVariables, numberOfObjectives);
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = numberOfObjectives;
		numberOfConstraints_ = 0;
		problemName_ = "LSMOP8";

		nk = 5;
		c = new double[numberOfObjectives_];
		c[0] = 3.8*0.1*(1.0-0.1);
		double sumC = c[0];
		for(int i = 1; i < numberOfObjectives_; i++){
			c[i] = 3.8*c[i-1]*(1.0-c[i-1]);
			sumC += c[i];
		}
		sublen = new int[numberOfObjectives_];
		for(int i = 0; i < numberOfObjectives_; i++){
			//sublen[i] = (int)Math.ceil(Math.round(numberOfVariables_*c[i]/sumC)/nk);
			//sublen[i] = (int)Math.floor(numberOfVariables_*c[i]/sumC/nk);
			sublen[i] = (int)Math.floor((numberOfVariables_-numberOfObjectives_+1)*c[i]/sumC/nk);
		}
		len = new int[numberOfObjectives_+1]; 
		len[0] = 0;
		for(int i = 0; i < numberOfObjectives_; i++){
			len[i+1] = len[i] + nk*sublen[i];
		}
		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		for (int var = 0; var < numberOfObjectives_ - 1; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for
		
		for (int var = numberOfObjectives_ - 1; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 10.0;
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
		Variable[] D = solution.getDecisionVariables();
		int dim = D.length;
		if(dim != numberOfVariables_) {
			System.out.println("LSMOP8: Recent Dim = "+ dim);
			System.out.println("LSMOP8: The Preset numberOfVariables = "+ numberOfVariables_);
		}
		double[] x = new double[dim];
		double[] f = new double[numberOfObjectives_];
		//Variable Linkage
		for (int i = 0; i < numberOfObjectives_-1; i++)
			x[i] = D[i].getValue();
		for(int i = numberOfObjectives_-1; i < dim; i++){
			x[i] = (1.0 + Math.cos(0.5*Math.PI*((i+1.0)/(dim))))*D[i].getValue() - 10.0*D[0].getValue();
		}
		
		double[] hx = new double[numberOfObjectives_]; 
		for (int i = 0; i < numberOfObjectives_; i++) {
			hx[i] = 1.0;
			for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
				hx[i] *= Math.cos(x[j] * 0.5 * Math.PI);
			if (i != 0) {
				int aux = numberOfObjectives_ - (i + 1);
				hx[i] *= Math.sin(x[aux] * 0.5 * Math.PI);
			} // if
		}// for
		
		double[] gx = new double[numberOfObjectives_];
		for (int i = 0; i < numberOfObjectives_; i++) {
			gx[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<nk;j++){
					double cg = 0.0;
					double cs = 1.0;
					int t = len[i]+numberOfObjectives_-1+j*sublen[i];
					for(;t<len[i]+numberOfObjectives_-1+(j+1)*sublen[i]-1;t++){
						cg += (x[t]*x[t])/4000;
						cs *= Math.cos(x[t]/Math.sqrt(t+1));
					}
					gx[i] += (cg - cs + 1)/sublen[i];
				}
			}else{
				for(int j=0;j<nk;j++){
					double cg = 0.0;
					int t = len[i]+numberOfObjectives_-1+j*sublen[i];
					for(;t<len[i]+numberOfObjectives_-1+(j+1)*sublen[i]-1;t++){
						cg+=x[t]*x[t];
					}
					cg = cg/sublen[i];
					gx[i] += cg;
				}
			}
			
			gx[i] = gx[i]/nk;
		}

		for (int i = 0; i < numberOfObjectives_; i++) {
			if(i != numberOfObjectives_-1){
				f[i] = hx[i]*(1+gx[i]+gx[i+1]);
			}else{
				f[i] = hx[i]*(1+gx[i]);
			}
			
		}// for

		for (int i = 0; i < numberOfObjectives_; i++)
			solution.setObjective(i,f[i]);
	} // evaluate

}//LSMOP8
