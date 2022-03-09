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

public class AGSMetaCrossover {
    private SolutionSet population;
    private MLPModel MLP;
    private int numVars;
    private int numObjs;

    static double learningRate = 0.1;
    static int K = 10;
    static int numLayers = 3;

    private double[] upBounds;
    private double[] lowBounds;

    private double[] zIdeal;
    private double[] zNadir;

    private Problem problem;
    static Distance distance = new Distance();

    public AGSMetaCrossover(SolutionSet population, Problem problem){
        this.population = population;
        this.problem = problem;
        this.numVars = problem.getNumberOfVariables();
        this.numObjs = problem.getNumberOfObjectives();
        MLP = new MLPModel(numVars, K, numLayers, learningRate);
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
    
    public AGSMetaCrossover(SolutionSet population, Problem problem, MLPModel currentMLPModel){
        this.population = population;
        this.problem = problem;
        this.numVars = problem.getNumberOfVariables();
        this.numObjs = problem.getNumberOfObjectives();
        MLP = currentMLPModel;
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

    public SolutionSet doCrossover(int type)throws JMException, ClassNotFoundException{
        int size = population.size();
        SolutionSet offspringPopulation;
        /*Data Processing*/
        SolutionSet st = getStSolutionSet(population, (int)(size*0.7));
        estimateIdealNadirPoint(st);
        //estimateIdealNadirPoint(population);
        normalization(population);
        /*Classification*/
        //SolutionSet[] trainingSet = new PopulationClassification(population, numObjs, lamda).classification();
        SolutionSet[] trainingSet = getClassifiedSet();
        /*Training*/
        SolutionSet labelSet = getLabel(trainingSet);
        MLP.getTrainingModel(trainingSet[1],labelSet,upBounds, lowBounds);

        offspringPopulation = reproduction(type);

        return offspringPopulation;
    }

    public SolutionSet reproduction(int type) throws JMException, ClassNotFoundException {
        int populationSize = population.size();
        SolutionSet offspringPopulation = new SolutionSet(populationSize);
        Solution[] parents = new Solution[3];
        Solution offSpring;
        for (int i = 0; i < (populationSize); i++) {
            offSpring = new Solution(population.get(i));
            int rdInt1 = PseudoRandom.randInt(0, populationSize-1);
            while(rdInt1 == i) {
                rdInt1 = PseudoRandom.randInt(0, populationSize-1);
            }
            int rdInt2 = PseudoRandom.randInt(0, populationSize-1);
            while(rdInt1 == rdInt2 || rdInt2 == i) {
                rdInt2 = PseudoRandom.randInt(0, populationSize-1);
            }
            parents[0] = population.get(rdInt1);
            parents[1] = population.get(rdInt2);
            parents[2] = population.get(i);
            XReal xChild = new XReal(offSpring);
            XReal xParent0 = new XReal(parents[0]);
            XReal xParent1 = new XReal(parents[1]);
            XReal xParent2 = new XReal(parents[2]);
            XReal xCurrent = new XReal(population.get(i));
            int jrand = PseudoRandom.randInt(0, numVars - 1);

            double[] xLearning = learningViaDAE(xChild);

            for (int j = 0; j < numVars; j++) {
                if (PseudoRandom.randDouble(0, 1) < 1.0 || j == jrand) {
                    double value;

                    double r1 = PseudoRandom.randDouble(0, 1);
                    double r2 = PseudoRandom.randDouble(0, 1);
                    if(type == 0){
                        value = xParent2.getValue(j) +r1*(xParent0.getValue(j) - xParent1.getValue(j))
                                + r2*(-xLearning[j] + xParent2.getValue(j));
                    }else{
                        value = xParent2.getValue(j) + r1*(xParent0.getValue(j) - xParent1.getValue(j));
                    }
                    if (value < lowBounds[j])
                        value = lowBounds[j];
                    if (value > upBounds[j])
                        value = upBounds[j];
                    xChild.setValue(j, value);
                } else {
                	double value;
                	if (PseudoRandom.randDouble(0, 1) < 0.5) {
                        double r1 = PseudoRandom.randDouble(0, 1);
                        value = xParent2.getValue(j) + r1*(xParent0.getValue(j) - xParent1.getValue(j));
                        if (value < lowBounds[j])
                            value = lowBounds[j];
                        if (value > upBounds[j])
                            value = upBounds[j];
                	}else {
                		value = xCurrent.getValue(j);
                	}
                	xChild.setValue(j, value);
                } // else
            } // for j
            doMutation(0.01,offSpring);
            problem.evaluate(offSpring);
            problem.evaluateConstraints(offSpring);
            offspringPopulation.add(offSpring);
        }//for i
        return offspringPopulation;
    }
    
    public MLPModel getCurrentMLPModel() {
    	return this.MLP;
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

    public double[] learningViaDAE(XReal xsol) throws JMException, ClassNotFoundException {
        double[] input = new double[numVars];
        for(int i=0;i<numVars;i++) {
            input[i] = (xsol.getValue(i) - lowBounds[i])/(upBounds[i] - lowBounds[i]);
        }
        double[] output = MLP.baseMLP.computeOut(input);
        double value;
        for(int var=0;var<numVars;var++) {
            value = output[var]*(upBounds[var] - lowBounds[var]) + lowBounds[var];
            if(value < lowBounds[var]) {
                value = lowBounds[var];
            }
            if(value > upBounds[var]) {
                value = upBounds[var];
            }
            output[var] = value;
        }
        return output;
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
            double dis2Ideal = 0.0;
            for(int i=0;i<numObjs;i++){
                value = (solSet.get(p).getObjective(i)-zIdeal[i])/
                        (zNadir[i]-zIdeal[i]);
                solSet.get(p).setNormalizedObjective(i, value);
                sum += value;
                dis2Ideal += value*value;
            }
            dis2Ideal = Math.sqrt(dis2Ideal);
            solSet.get(p).setSumValue(sum);
            solSet.get(p).setDistanceToIdealPoint(dis2Ideal);
        }
    }

    public SolutionSet getLabel(SolutionSet[] input){
        int size = input[1].size();
        SolutionSet labelSet = new SolutionSet(size);
        for(int i=0;i<size;i++) {
        	int rdInt = PseudoRandom.randInt(0, input[0].size() - 1);
        	if(PseudoRandom.randDouble(0, 1) < 0.0) {
        		labelSet.add(input[0].get(rdInt));
        	}else {
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

}
