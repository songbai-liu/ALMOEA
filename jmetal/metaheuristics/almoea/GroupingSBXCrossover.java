package jmetal.metaheuristics.almoea;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.wrapper.XReal;

public class GroupingSBXCrossover {
    private SolutionSet population;
    private int numVars;
    private int numObjs;
    private int[][] group_;
    private int groupSize_;
    private int populationSize;

    private double[] upBounds;
    private double[] lowBounds;

    private Problem problem;

    public GroupingSBXCrossover(SolutionSet population, Problem problem){
        this.population = population;
        this.problem = problem;
        this.numVars = problem.getNumberOfVariables();
        this.numObjs = problem.getNumberOfObjectives();
        upBounds = new double[numVars];
        lowBounds = new double[numVars];
        for(int var=0;var<numVars;var++) {
            upBounds[var] = problem.getUpperLimit(var);
            lowBounds[var] = problem.getLowerLimit(var);
        }
        groupSize_ = 4;
        populationSize = population.size();
    }
    
    public GroupingSBXCrossover(SolutionSet population, Problem problem, int groupSize){
        this.population = population;
        this.problem = problem;
        this.numVars = problem.getNumberOfVariables();
        this.numObjs = problem.getNumberOfObjectives();
        upBounds = new double[numVars];
        lowBounds = new double[numVars];
        for(int var=0;var<numVars;var++) {
            upBounds[var] = problem.getUpperLimit(var);
            lowBounds[var] = problem.getLowerLimit(var);
        }
        groupSize_ = groupSize;
        populationSize = population.size();
    }

    public SolutionSet doCrossover()throws JMException, ClassNotFoundException{
        SolutionSet offspringPopulation;
        linearGrouping();
        offspringPopulation = reproduction();
        return offspringPopulation;
    }

    public SolutionSet reproduction() throws JMException, ClassNotFoundException {
        SolutionSet offspringPopulation = new SolutionSet(populationSize);
        Solution[] parents = new Solution[2];
        for (int i = 0; i < populationSize; i++) {
    		int rdInt1 = PseudoRandom.randInt(0, populationSize-1);
    		int rdInt2 = PseudoRandom.randInt(0, populationSize-1);
    		while(rdInt2 == rdInt1) {
    			rdInt2 = PseudoRandom.randInt(0, populationSize-1);
    		}
    		parents[0] = population.get(rdInt1);
    		parents[1] = population.get(rdInt2);
            XReal[] xParents = new XReal[2];
			for(int p=0; p<2; p++) {
				xParents[p] = new XReal(parents[p]);
			}
			Solution child = new Solution(population.get(i));
			double learningRate = PseudoRandom.randDouble(0.0, 1.0);
			//encode
			double[][] encodedParents = encode(xParents);
            //search in the compressed space
			double[] newEncode = searchInTransferedSpace(encodedParents,learningRate);
			//decode
			decode(newEncode,child);
            doMutation(0.001,child);
            problem.evaluate(child);
            problem.evaluateConstraints(child);
            child.setindex(1);
            offspringPopulation.add(child);
        }//for i

        return offspringPopulation;
    }

    public void doMutation(double probability, Solution solution)
            throws JMException {
        double rnd, delta1, delta2, mut_pow, deltaq;
        double y, yl, yu, val, xy;
        XReal x = new XReal(solution);
        for (int var = 0; var < solution.numberOfVariables(); var++) {
            if (PseudoRandom.randDouble() <= probability) {
                y = x.getValue(var);
                yl = x.getLowerBound(var);
                yu = x.getUpperBound(var);
                delta1 = (y - yl) / (yu - yl);
                delta2 = (yu - y) / (yu - yl);
                rnd = PseudoRandom.randDouble();
                mut_pow = 1.0 / (20 + 1.0);
                if (rnd <= 0.5) {
                    xy = 1.0 - delta1;
                    val = 2.0 * rnd + (1.0 - 2.0 * rnd)
                            * (Math.pow(xy, (20 + 1.0)));
                    deltaq = java.lang.Math.pow(val, mut_pow) - 1.0;
                } else {
                    xy = 1.0 - delta2;
                    val = 2.0
                            * (1.0 - rnd)
                            + 2.0
                            * (rnd - 0.5)
                            * (java.lang.Math.pow(xy,
                            (20 + 1.0)));
                    deltaq = 1.0 - (java.lang.Math.pow(val, mut_pow));
                }
                y = y + deltaq * (yu - yl);
                if (y < yl)
                    y = yl;
                if (y > yu)
                    y = yu;
                x.setValue(var, y);
            }
        } // for

    } // doMutation

    public double[][] encode(XReal[] parents) throws JMException {
		int len = parents.length;
		double[][] encodedParents = new double[len][groupSize_];
		for(int l=0;l<len;l++) {
			double[] encodedSolution = new double[numVars];
			encodedSolution = normalizedVariable(parents[l]);
			for(int g=0;g<groupSize_;g++) {
				encodedParents[l][g] = 0;
				for(int var=0;var<group_[g].length;var++) {
					encodedParents[l][g] += encodedSolution[group_[g][var]];
				}
				encodedParents[l][g] = encodedParents[l][g]/group_[g].length;
			}
		}
		return encodedParents; 
	}
	
	public double[] searchInTransferedSpace(double[][] encodedParents, double learningRate) {
		double[][] newEncode = new double[2][groupSize_];
		double rand;
		double EPS = 1.0e-14;
		double distributionIndex_ = 20.0;
		double y1, y2, yL, yu;
		double c1, c2;
		double alpha, beta, betaq;
		double valueX1, valueX2;
		for(int j=0;j<groupSize_;j++) {
			valueX1 = encodedParents[0][j];
			valueX2 = encodedParents[1][j];
			if (PseudoRandom.randDouble() <= 0.5) {
				if (java.lang.Math.abs(valueX1 - valueX2) > EPS) {

					if (valueX1 < valueX2) {
						y1 = valueX1;
						y2 = valueX2;
					} else {
						y1 = valueX2;
						y2 = valueX1;
					} // if

					yL = lowBounds[j];
					yu = upBounds[j];
					rand = PseudoRandom.randDouble();
					beta = 1.0 + (2.0 * (y1 - yL) / (y2 - y1));
					alpha = 2.0 - java.lang.Math.pow(beta,
							-(distributionIndex_ + 1.0));

					if (rand <= (1.0 / alpha)) {
						betaq = java.lang.Math.pow((rand * alpha),
								(1.0 / (distributionIndex_ + 1.0)));
					} else {
						betaq = java.lang.Math.pow((1.0 / (2.0 - rand
								* alpha)),
								(1.0 / (distributionIndex_ + 1.0)));
					} // if

					c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
					beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
					alpha = 2.0 - java.lang.Math.pow(beta,
							-(distributionIndex_ + 1.0));

					if (rand <= (1.0 / alpha)) {
						betaq = java.lang.Math.pow((rand * alpha),
								(1.0 / (distributionIndex_ + 1.0)));
					} else {
						betaq = java.lang.Math.pow((1.0 / (2.0 - rand
								* alpha)),
								(1.0 / (distributionIndex_ + 1.0)));
					} // if

					c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

					if (c1 < yL)
						c1 = yL;

					if (c2 < yL)
						c2 = yL;

					if (c1 > yu)
						c1 = yu;

					if (c2 > yu)
						c2 = yu;

					if (PseudoRandom.randDouble() <= 0.5) {
						newEncode[0][j] = c2;
						newEncode[1][j] = c1;
					} else {
						newEncode[0][j] = c1;
						newEncode[1][j] = c2;
					} // if
				} else {
					newEncode[0][j] = valueX1;
					newEncode[1][j] = valueX2;
				} // if
			} else {
				newEncode[0][j] = valueX2;
				newEncode[1][j] = valueX1;
			} // if

			if(newEncode[0][j] <= 0) { newEncode[0][j] = 0.00001; } 
			if(newEncode[0][j] >= 1) { newEncode[0][j] = 0.99999; }
		}
		return newEncode[0];
	}
	
	public void decode(double[] newEncode, Solution child) throws JMException {
		XReal offspring = new XReal(child);
		for(int g=0;g<groupSize_;g++) {
			double sum = 0;
			for(int var=0;var<group_[g].length;var++) {
				int cVar = group_[g][var];
				sum += (offspring.getValue(cVar)-lowBounds[cVar])/(upBounds[cVar] - lowBounds[cVar]);
			}
			if(sum == 0) sum = 0.00001;
			for(int var=0;var<group_[g].length;var++) {
				int cVar = group_[g][var];
				double normalizedValue = (offspring.getValue(cVar)-lowBounds[cVar])/(upBounds[cVar] - lowBounds[cVar]);
				double decodedValue = (normalizedValue/sum) * newEncode[g] * group_[g].length;
				decodedValue = decodedValue*(upBounds[cVar] - lowBounds[cVar]) + lowBounds[cVar];
				if(decodedValue < lowBounds[cVar]) {
					decodedValue = lowBounds[cVar];
				}
				if(decodedValue > upBounds[cVar]) {
					decodedValue = upBounds[cVar];
				}
				offspring.setValue(cVar, decodedValue);
			}
		}
	}
	
	public double[] normalizedVariable(XReal solution) throws JMException {
		double[] x = new double[numVars];
		for (int i = 0; i < numVars; i++) 
			x[i] = solution.getValue(i);
		for (int i = 0; i < numVars; i++) {
			x[i] = (x[i] - lowBounds[i])  / (upBounds[i] - lowBounds[i]);
		}
		return x;
	}
	
    public void linearGrouping(){
		int gSize_ = numVars/groupSize_;
		group_ = new int[groupSize_][];
		for(int g=0;g<groupSize_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVars-(groupSize_-1)*gSize_;//the variable size of the last group
		group_[groupSize_-1] = new int[lSize_];
		
		int t = 0;
		for(int g=0;g<groupSize_-1;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = t;
				t++;
			}
		}
		//assign variable to the last group
		for(int m=0;m<lSize_;m++){
			group_[groupSize_-1][m] = t;
			t++;
		}
	}

}
