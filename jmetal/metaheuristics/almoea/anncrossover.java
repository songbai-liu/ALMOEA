package jmetal.metaheuristics.almoea;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class anncrossover {
    private SolutionSet population;
    private MLPModel DAE;
    private int numVars;
    private int numObjs;

    static double learningRate = 0.1;
    static int reD = 10;
    static int numLayers = 3;
    static int div1_ = 49;
    static int div2_ = 0;

    private double[] upBounds;
    private double[] lowBounds;

    private double[] zIdeal;
    private double[] zNadir;
    private double[][] lamda;
    private VectorGenerator vg;

    private Problem problem;

    public anncrossover(SolutionSet population, Problem problem){
        this.population = population;
        this.problem = problem;
        this.numVars = problem.getNumberOfVariables();
        this.numObjs = problem.getNumberOfObjectives();
        DAE = new MLPModel(numVars, reD, numLayers, learningRate);
        vg = new TwoLevelWeightVectorGenerator(div1_, div2_,numObjs);
        lamda = vg.getVectors();
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
        /*Data Processing*/
        estimateIdealNadirPoint(population);
        normalization(population);
        /*Classification*/
        SolutionSet[] trainingSet = new PopulationClassification(population, numObjs, lamda).classification();
        /*Training*/
        SolutionSet labelSet = getLabel(trainingSet);
        DAE.getTrainingModel(trainingSet[1],labelSet,upBounds, lowBounds);

        offspringPopulation = reproduction();

        return offspringPopulation;
    }

    public SolutionSet reproduction() throws JMException, ClassNotFoundException {
        int populationSize = population.size();
        SolutionSet offspringPopulation = new SolutionSet(populationSize);
        Solution[] parents = new Solution[3];
        Solution offSpring;
        double CR = 1.0;
        double F = 0.5;
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
                if (PseudoRandom.randDouble(0, 1) < CR || j == jrand) {
                    double value;
                    value=PseudoRandom.randDouble(0, 1);
                    value = xParent2.getValue(j) + PseudoRandom.randDouble(0, 1)*(xParent0.getValue(j) - xParent1.getValue(j))
                            + PseudoRandom.randDouble(0, 1)*(-xLearning[j] + xParent2.getValue(j));
                    //value = xParent2.getValue(j) + 0.5*(-xLearning[j] + xParent2.getValue(j));
                    if (value < lowBounds[j])
                        value = lowBounds[j];
                    if (value > upBounds[j])
                        value = upBounds[j];
                    xChild.setValue(j, value);
                } else {
                    double value;
                    value = xCurrent.getValue(j);
                    xChild.setValue(j, value);
                } // else
            } // for j

            //doMutation(0.05,offSpring);
            problem.evaluate(offSpring);
            problem.evaluateConstraints(offSpring);
            offspringPopulation.add(offSpring);
        }//for i
        return offspringPopulation;
    }

    public double[] learningViaDAE(XReal xsol) throws JMException, ClassNotFoundException {
        double[] input = new double[numVars];
        for(int i=0;i<numVars;i++) {
            input[i] = (xsol.getValue(i) - lowBounds[i])/(upBounds[i] - lowBounds[i]);
        }
        double[] output = DAE.baseMLP.computeOut(input);
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



}
