package jmetal.problems.LMF;
import java.io.IOException;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.util.JMException;

public class LMF8 extends LMF{
	
	public LMF8(String solutionType, int numberOfVariables, int numberOfObjectives) {
		super(solutionType, numberOfVariables, numberOfObjectives);
		problemName_ = "LMF8";
		gType_ = "LMF8";
		hType_ = "linear";
	}
	
	public LMF8(int numberOfVariables, int numberOfObjectives, int type1, int type2, int type3, int type4) {
		super("Real", numberOfVariables, numberOfObjectives, type1, type2, type3, type4);
		problemName_ = "LMF8";
		gType_ = "LMF8";
		hType_ = "linear";
	}
	
	public void evaluate(Solution solution) throws JMException {
		Variable[] D = solution.getDecisionVariables();
		int dim = D.length;
		double[] vars = new double[dim];
		for (int i = 0; i < dim; i++)
			vars[i] = D[i].getValue();
		
		double[][] xI = new double[numberOfObjectives_ - 1][];
		double[][][] xII = new double[numberOfObjectives_][][];

		for (int i = 0; i < numberOfObjectives_ - 1; i++){
			xI[i] = new double[gp[i].length];
			for(int j=0;j<gp[i].length;j++){
				xI[i][j] = vars[gp[i][j]];
			}
		}
		double[] y = new double[numberOfObjectives_-1];
		for(int i=0;i<numberOfObjectives_-1;i++){
			y[i] = sum_avg_abs(xI[i]);
		}

		for (int i = 0; i < numberOfObjectives_; i++){
			xII[i] = new double[dgd[i].length][];
			for(int j=0;j<dgd[i].length;j++){
				xII[i][j] = new double[dgd[i][j].length];
				for(int k=0;k<dgd[i][j].length;k++){
					double sb = vars[dgd[i][j][k]];
					int index = dgd[i][j][k];
					if(type_lk == 0) {
						xII[i][j][k] = (1.0 + (index+1.0)/(L))*sb - 10.0*y[0];
					}else if(type_lk == 1) {
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index+1.0)/(L))))*sb - 10.0*y[0];
					}else {
						if(i % 2 == 0) {
							xII[i][j][k] = (1.0 + (index+1.0)/(L))*sb - 10.0*y[0];
						}else {
							xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index+1.0)/(L))))*sb - 10.0*y[0];
						}
					}
					
				}
			}
		}
		double[] h = evalH(y);
		double[] g = evalG(y[0],xII);		
		for (int i = 0; i < numberOfObjectives_; i++) {
			if(type_fm == 0) {//addition model
				solution.setObjective(i, h[i] + g[i]);
			}else if(type_fm == 1) {//multiplication model
				solution.setObjective(i, h[i]*(1 + g[i]));
			}else {//Mixed model
				if(i % 2 == 0) {
					solution.setObjective(i, h[i]*(1 + g[i]));
				}else {
					solution.setObjective(i, h[i] + g[i]);
				}
			}
		}
	}
}
