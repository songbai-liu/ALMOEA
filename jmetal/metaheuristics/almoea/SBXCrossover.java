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

public class SBXCrossover {
    private SolutionSet population;
    private int numVars;
    private int numObjs;


    private double[] upBounds;
    private double[] lowBounds;

    private double[] zIdeal;
    private double[] zNadir;

    private Problem problem;
    static Distance distance = new Distance();

    public SBXCrossover(SolutionSet population, Problem problem){
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
        zIdeal = new double[numObjs];
        zNadir = new double[numObjs];
        for(int m=0;m<numObjs;m++){
            zIdeal[m] = Double.MAX_VALUE;
            zNadir[m] = Double.MIN_VALUE;
        }
    }

    public SolutionSet doCrossover()throws JMException, ClassNotFoundException{
        int size = population.size();
        SolutionSet offspringPopulation;
        offspringPopulation = reproduction();

        return offspringPopulation;
    }

    public SolutionSet reproduction() throws JMException, ClassNotFoundException {
    	double EPS = 1.0e-14;
    	double probability = 1.0;
    	double distributionIndex_ = 20.0;
        int populationSize = population.size();
        SolutionSet offspringPopulation = new SolutionSet(populationSize);
        Solution[] parents = new Solution[2];
        Solution[] offSpring;
        for (int i = 0; i < population.size(); i++) {  
            int rdInt1 = PseudoRandom.randInt(0, population.size()-1);
            int rdInt2 = PseudoRandom.randInt(0, population.size()-1);
            while(rdInt1 == rdInt2) {
                rdInt2 = PseudoRandom.randInt(0, population.size()-1);
            }
            parents[0] = population.get(rdInt1);
            parents[1] = population.get(rdInt2);
            offSpring = new Solution[2];
    		offSpring[0] = new Solution(parents[0]);
    		offSpring[1] = new Solution(parents[1]);
    		int j;
    		double rand;
    		double y1, y2, yL, yu;
    		double c1, c2;
    		double alpha, beta, betaq;
    		double valueX1, valueX2;
    		XReal x1 = new XReal(parents[0]);
    		XReal x2 = new XReal(parents[1]);
    		XReal offs1 = new XReal(offSpring[0]);
    		XReal offs2 = new XReal(offSpring[1]);

    		int numberOfVariables = x1.getNumberOfDecisionVariables();

    		if (PseudoRandom.randDouble() <= probability) {
    			for (j = 0; j < numberOfVariables; j++) {
    				valueX1 = x1.getValue(j);
    				valueX2 = x2.getValue(j);
    				if (PseudoRandom.randDouble() <= 0.5) {
    					if (java.lang.Math.abs(valueX1 - valueX2) > EPS) {

    						if (valueX1 < valueX2) {
    							y1 = valueX1;
    							y2 = valueX2;
    						} else {
    							y1 = valueX2;
    							y2 = valueX1;
    						} // if

    						yL = x1.getLowerBound(j);
    						yu = x1.getUpperBound(j);
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
    							offs1.setValue(j, c2);
    							offs2.setValue(j, c1);
    						} else {
    							offs1.setValue(j, c1);
    							offs2.setValue(j, c2);
    						} // if
    					} else {
    						offs1.setValue(j, valueX1);
    						offs2.setValue(j, valueX2);
    					} // if
    				} else {
    					offs1.setValue(j, valueX2);
    					offs2.setValue(j, valueX1);
    				} // if
    			} // if
    		} // if
            doMutation(0.001,offSpring[0]);
            problem.evaluate(offSpring[0]);
            problem.evaluateConstraints(offSpring[0]);
            offSpring[0].setindex(1);
            offspringPopulation.add(offSpring[0]);
        }//for i

        return offspringPopulation;
    }

    public SolutionSet getStSolutionSet(SolutionSet ss,int size) {
        Ranking ranking = new NondominatedRanking(ss);
        int remain = size;
        int index = 0;
        SolutionSet front = null;
        SolutionSet mgPopulation = new SolutionSet();
        front = ranking.getSubfront(index);
        while ((remain > 0) && (remain >= front.size())) {

            for (int k = 0; k < front.size(); k++) {
                mgPopulation.add(front.get(k));
            } // for
            // Decrement remain
            remain = remain - front.size();
            // Obtain the next front
            index++;
            if (remain > 0) {
                front = ranking.getSubfront(index);
            } // if
        }
        if (remain > 0) { // front contains individuals to insert
            for (int k = 0; k < front.size(); k++) {
                mgPopulation.add(front.get(k));
            }
        }
        return mgPopulation;
    }

    public void estimateIdealNadirPoint(SolutionSet solSet){
        for(int p=0;p<solSet.size();p++){
            for(int i=0;i<numObjs;i++){
                if(solSet.get(p).getObjective(i) < zIdeal[i])
                    zIdeal[i] = solSet.get(p).getObjective(i);
                if(solSet.get(p).getObjective(i) > zNadir[i])
                    zNadir[i] = solSet.get(p).getObjective(i);
            }
        }
    }

    public void normalization(SolutionSet solSet){
        double value;
        for(int p=0;p<solSet.size();p++){
            double sum = 0.0;
            for(int i=0;i<numObjs;i++){
                value = (solSet.get(p).getObjective(i)-zIdeal[i])/
                        (zNadir[i]-zIdeal[i]);
                solSet.get(p).setNormalizedObjective(i, value);
                sum += value;
            }
            solSet.get(p).setSumValue(sum);
        }
    }

    public SolutionSet getLabel(SolutionSet[] input){
        int size = input[1].size();
        SolutionSet labelSet = new SolutionSet(size);
        for(int i=0;i<size;i++) {
            Solution sol = input[1].get(i);
            double minDis = distance(sol, input[0].get(0));
            int minIndex = 0;
            for(int j=1;j<input[0].size();j++) {
                double dis = distance(sol, input[0].get(j));
                if(dis < minDis) {
                    minDis = dis;
                    minIndex = j;
                }
            }
            labelSet.add(input[0].get(minIndex));
        }
        return labelSet;
    }

    public double distance(Solution s1, Solution s2) {
        double dis = 0.0;
        double sum1 = s1.getSumValue();
        double sum2 = s2.getSumValue();
        double v1,v2;
        for(int i=0;i<numObjs;i++) {
            v1 = s1.getNormalizedObjective(i)/sum1;
            v2 = s2.getNormalizedObjective(i)/sum2;
            dis += (v1 - v2) * (v1 - v2);
        }
        dis = Math.sqrt(dis);
        return dis;
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

    public SolutionSet[] getClassifiedSet(){
        int populationSize = population.size();
        SolutionSet[] st = new SolutionSet[2];
        NondominatedRanking ranking = new NondominatedRanking(population);
        int size = populationSize/2;
        st[0] = new SolutionSet(size);
        st[1] = new SolutionSet(size);
        int remain = size;
        int index = 0;
        SolutionSet front = null;
        // Obtain the next front
        front = ranking.getSubfront(index);

        while ((remain > 0) && (remain >= front.size())) {
            // Assign crowding distance to individuals
            distance.crowdingDistanceAssignment(front,
                    problem.getNumberOfObjectives());
            // Add the individuals of this front
            for (int k = 0; k < front.size(); k++) {
                st[0].add(front.get(k));
            } // for

            // Decrement remain
            remain = remain - front.size();

            // Obtain the next front
            index++;
            if (remain > 0) {
                front = ranking.getSubfront(index);
            } // if
        } // while

        // Remain is less than front(index).size, insert only the best one
        if (remain > 0) { // front contains individuals to insert
            distance.crowdingDistanceAssignment(front,
                    problem.getNumberOfObjectives());
            front.sort(new CrowdingComparator());
            for (int k = 0; k < remain; k++) {
                st[0].add(front.get(k));
            } // for
            for (int k = remain; k < front.size(); k++) {
                st[1].add(front.get(k));
            } // for
            remain = 0;
            index++;
        } // if

        for(int i=index;i<ranking.getNumberOfSubfronts();i++){
            front = ranking.getSubfront(i);
            for (int k = 0; k < front.size(); k++) {
                st[1].add(front.get(k));
            } // for
        }
        return st;
    }
    
    public SolutionSet[] getClassifiedSet1(){
        int populationSize = population.size();
        SolutionSet[] st = new SolutionSet[2];
        NondominatedRanking ranking = new NondominatedRanking(population);
        SolutionSet front0 = ranking.getSubfront(0);
        st[0] = front0;
        st[1] = new SolutionSet(populationSize-front0.size());
        for(int i=1;i<ranking.getNumberOfSubfronts();i++){
            SolutionSet front = ranking.getSubfront(i);
            for (int k = 0; k < front.size(); k++) {
                st[1].add(front.get(k));
            } // for
        } 
        return st;
    }

}
