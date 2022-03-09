package jmetal.metaheuristics.almoea;

import java.util.ArrayList;
import java.util.List;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.almoea.BaseMLP.IHierarchicalClustering;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.comparators.SumValueComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.wrapper.XReal;

public class AGSBalanceableCrossover {
    private SolutionSet population;
    private MLPModel cMLP;
    private MLPModel dMLP;
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

    public AGSBalanceableCrossover(SolutionSet population, Problem problem){
        this.population = population;
        this.problem = problem;
        this.numVars = problem.getNumberOfVariables();
        this.numObjs = problem.getNumberOfObjectives();
        cMLP = new MLPModel(numVars, K, numLayers, learningRate);
        dMLP = new MLPModel(numVars, K, numLayers, learningRate);
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
    
    public AGSBalanceableCrossover(SolutionSet population, Problem problem, MLPModel cMLP, MLPModel dMLP){
        this.population = population;
        this.problem = problem;
        this.numVars = problem.getNumberOfVariables();
        this.numObjs = problem.getNumberOfObjectives();
        this.cMLP = cMLP;
        this.dMLP = dMLP;
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
        SolutionSet st = getStSolutionSet(population, (int)(size*0.85));
        estimateIdealNadirPoint(st);
        normalization(population);
        /*Classification*/
        SolutionSet[] cTrainingSet = getConClassifiedSet();
        SolutionSet[] dTrainingSet = getDivClassifiedSet();
        /*Training*/
        //SolutionSet cLabelSet = getConLabel(cTrainingSet);
        cMLP.getTrainingModel(cTrainingSet[1],cTrainingSet[2],upBounds, lowBounds);
        //SolutionSet dLabelSet = getConLabel(dTrainingSet);
        dMLP.getTrainingModel(dTrainingSet[1],dTrainingSet[2],upBounds, lowBounds);
        offspringPopulation = reproduction(cTrainingSet[0],dTrainingSet[0]);

        return offspringPopulation;
    }
    
    public SolutionSet reproduction(SolutionSet CS, SolutionSet DS) throws JMException, ClassNotFoundException {
        int populationSize = population.size();
        double[] xCLearning,xDLearning;
        SolutionSet offspringPopulation = new SolutionSet(populationSize);
        Solution[] parents = new Solution[4];
        Solution offSpring;
        double CR = 1.0;
        double F = 0.33;
        for (int i = 0; i < populationSize; i++) {
            offSpring = new Solution(population.get(i));
            int rdInt1 = PseudoRandom.randInt(0, populationSize-1);
            while(rdInt1 == i) {
            	rdInt1 = PseudoRandom.randInt(0, populationSize-1);
            }
            
            int rdInt2 = PseudoRandom.randInt(0, populationSize-1);
            while(rdInt2 == rdInt1 || rdInt2 == i) {
            	rdInt2 = PseudoRandom.randInt(0, populationSize-1);
            }
            parents[0] = population.get(rdInt1);
            parents[1] = population.get(rdInt2);

            XReal xChild = new XReal(offSpring);
            XReal xParent0 = new XReal(parents[0]);
            XReal xParent1 = new XReal(parents[1]);
            XReal xCurrent = new XReal(population.get(i));
            
            if(offSpring.getClusterID() == 0 && offSpring.getindex() == 0){//Case 1
            	xCLearning = learningViaDAE(xChild,cMLP);
       		 	xDLearning = learningViaDAE(xChild,dMLP);
       		 	for (int j = 0; j < numVars; j++) {
            		double value;
                    value = PseudoRandom.randDouble(0, 1);
                    value = xCurrent.getValue(j) + 0*(-xCLearning[j] + xCurrent.getValue(j)) 
                    		+ F*(-xDLearning[j] + xCurrent.getValue(j))
                    		+ F*(xParent0.getValue(j) - xParent1.getValue(j));
                    if (value < lowBounds[j])
                        value = lowBounds[j];
                    if (value > upBounds[j])
                        value = upBounds[j];
                    xChild.setValue(j, value);
            	}
            }else if(offSpring.getClusterID() == 1 && offSpring.getindex() == 1) {//case 2
            	int rd = PseudoRandom.randInt(0, DS.size()-1);
            	XReal xD = new XReal(DS.get(rd));
            	rd = PseudoRandom.randInt(0, CS.size()-1);
            	XReal xC = new XReal(CS.get(rd));
            	for (int j = 0; j < numVars; j++) {
            		double value;
                    value = PseudoRandom.randDouble(0, 1);
                    value = xCurrent.getValue(j) + F*(xD.getValue(j) - xCurrent.getValue(j)) 
                    		+ F*(xC.getValue(j) - xCurrent.getValue(j))
                    		+ F*(xParent0.getValue(j) - xParent1.getValue(j));
                    if (value < lowBounds[j])
                        value = lowBounds[j];
                    if (value > upBounds[j])
                        value = upBounds[j];
                    xChild.setValue(j, value);
            	}
            }else if(offSpring.getClusterID() == 0 && offSpring.getindex() == 1) {//Case 3
            	xCLearning = learningViaDAE(xChild,cMLP);
            	int rd = PseudoRandom.randInt(0, DS.size()-1);
            	XReal xD = new XReal(DS.get(rd));
            	for (int j = 0; j < numVars; j++) {
            		double value;
                    value = PseudoRandom.randDouble(0, 1);
                    value = xCurrent.getValue(j) + F*(xD.getValue(j) - xCurrent.getValue(j)) 
                    		+ F*(-xCLearning[j] + xCurrent.getValue(j))
                    		+ F*(xParent0.getValue(j) - xParent1.getValue(j));
                    if (value < lowBounds[j])
                        value = lowBounds[j];
                    if (value > upBounds[j])
                        value = upBounds[j];
                    xChild.setValue(j, value);
            	}
            }else {//Case 4
            	xDLearning = learningViaDAE(xChild,dMLP);
            	int rd = PseudoRandom.randInt(0, CS.size()-1);
            	XReal xC = new XReal(CS.get(rd));
            	for (int j = 0; j < numVars; j++) {
            		double value;
                    value = PseudoRandom.randDouble(0, 1);
                    value = xCurrent.getValue(j) + F*(xC.getValue(j) - xCurrent.getValue(j)) 
                    		+ F*(-xDLearning[j] + xCurrent.getValue(j))
                    		+ F*(xParent0.getValue(j) - xParent1.getValue(j));
                    if (value < lowBounds[j])
                        value = lowBounds[j];
                    if (value > upBounds[j])
                        value = upBounds[j];
                    xChild.setValue(j, value);
            	}
            }
            doMutation(0.01,offSpring);
            problem.evaluate(offSpring);
            problem.evaluateConstraints(offSpring);
            offspringPopulation.add(offSpring);
        }//for i

        return offspringPopulation;
    }

    public MLPModel getConvergenceMLPModel() {
    	return this.cMLP;
    }
    
    public MLPModel getDiversityMLPModel() {
    	return this.dMLP;
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

    public double[] learningViaDAE(XReal xsol, MLPModel MLP) throws JMException, ClassNotFoundException {
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
            double disToIdeal = 0.0;
            for(int i=0;i<numObjs;i++){
                value = (solSet.get(p).getObjective(i)-zIdeal[i])/
                        (zNadir[i]-zIdeal[i]);
                solSet.get(p).setNormalizedObjective(i, value);
                sum += value;
                disToIdeal += value*value;
            }
            disToIdeal = Math.sqrt(disToIdeal);
            solSet.get(p).setSumValue(sum);
            solSet.get(p).setFitness(sum-numObjs);
            solSet.get(p).setDistanceToIdealPoint(disToIdeal);
            for(int i=0;i<numObjs;i++){
            	solSet.get(p).setIthTranslatedObjective(i, solSet.get(p).getNormalizedObjective(i)/sum);
            }
        }
    }

    public SolutionSet getConLabel(SolutionSet[] input){
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

    public SolutionSet[] getConClassifiedSet(){
    	int populationSize = population.size();
        SolutionSet[] st = new SolutionSet[3];
        st[0] = new SolutionSet();
        st[1] = new SolutionSet();
        st[2] = new SolutionSet();
        List<SolutionSet> list = new <SolutionSet>ArrayList();
		for(int i=0;i<populationSize;i++){
			SolutionSet sols = new SolutionSet();
			sols.add(population.get(i));
		    list.add(sols);  
		}
		list = new IHierarchicalClustering(list).clusteringAnalysis(populationSize);
		for(int i=0;i<list.size();i++) {
			if(list.get(i).size() > 1) {
				list.get(i).sort(new FitnessComparator(false));
				list.get(i).get(list.get(i).size()-1).setClusterID(0);
				st[0].add(list.get(i).get(list.get(i).size()-1));
				for(int j=0; j<list.get(i).size()-1; j++) {
					st[1].add(list.get(i).get(j));
					list.get(i).get(j).setClusterID(1);
					int rd = PseudoRandom.randInt(j+1, list.get(i).size() - 1);
					st[2].add(list.get(i).get(rd));
				}
			}else {
				list.get(i).get(0).setClusterID(0);
				st[0].add(list.get(i).get(0));
			}
		}
        return st;
    }
    
    public SolutionSet[] getDivClassifiedSet(){
    	int populationSize = population.size();
        SolutionSet[] st = new SolutionSet[3];
        st[0] = new SolutionSet();
        st[1] = new SolutionSet();
        st[2] = new SolutionSet();
        List<SolutionSet> list = new <SolutionSet>ArrayList();
		for(int i=0;i<populationSize;i++){
			SolutionSet sols = new SolutionSet();
			sols.add(population.get(i));
		    list.add(sols);  
		}
		list = new HierarchicalClustering(list).clusteringAnalysis(populationSize);
		for(int i=0;i<list.size();i++) {
			int size = list.get(i).size();
			if(size > 1) {
				list.get(i).sort(new SumValueComparator(false));
				list.get(i).get(size-1).setindex(0);
				st[0].add(list.get(i).get(size-1));
				for(int j=0; j<size-1; j++) {
					st[1].add(list.get(i).get(j));
					list.get(i).get(j).setindex(1);
					int rd = PseudoRandom.randInt(j+1, size - 1);
					st[2].add(list.get(i).get(rd));
				}
			}else {
				list.get(i).get(0).setindex(0);
				st[0].add(list.get(i).get(0));
			}
		}
        return st;
    }

}
