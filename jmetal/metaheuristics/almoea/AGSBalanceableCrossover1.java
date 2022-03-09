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

public class AGSBalanceableCrossover1 {
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

    public AGSBalanceableCrossover1(SolutionSet population, Problem problem){
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
    
    public AGSBalanceableCrossover1(SolutionSet population, Problem problem, MLPModel cMLP, MLPModel dMLP){
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
        offspringPopulation = reproduction1(cTrainingSet[0],dTrainingSet[0]);

        return offspringPopulation;
    }
    
    public SolutionSet reproduction(SolutionSet CS, SolutionSet DS) throws JMException, ClassNotFoundException {
        int populationSize = population.size();
        double[] xCLearning,xDLearning;
        SolutionSet offspringPopulation = new SolutionSet(populationSize);
        Solution[] parents = new Solution[4];
        Solution offSpring;
        for (int i = 0; i < populationSize; i++) {
            offSpring = new Solution(population.get(i));
            XReal xChild = new XReal(offSpring);
            XReal xCurrent = new XReal(population.get(i));
            if(PseudoRandom.randDouble(0, 1) < 0.5) {
            	XReal xCS = new XReal(CS.get(PseudoRandom.randInt(0, CS.size() -1)));
        		XReal xDS = new XReal(DS.get(PseudoRandom.randInt(0, DS.size() -1)));
            	for (int j = 0; j < numVars; j++) {
            		double value = PseudoRandom.randDouble(0, 1);
            		value = xCurrent.getValue(j) + 0.5*(xDS.getValue(j) - xCS.getValue(j));
            		if (value < lowBounds[j])
                        value = lowBounds[j];
                    if (value > upBounds[j])
                        value = upBounds[j];
                    xChild.setValue(j, value);
            	}
            }else {
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

                
                 XReal xParent0 = new XReal(parents[0]);
                 XReal xParent1 = new XReal(parents[1]);
                 
                 xCLearning = learningViaDAE(xChild,cMLP);
        		 	xDLearning = learningViaDAE(xChild,dMLP);
        		 	for (int j = 0; j < numVars; j++) {
             		double value = PseudoRandom.randDouble(0, 1);
             		double r1 = PseudoRandom.randDouble(0.2, 0.6);
             		double r2 = PseudoRandom.randDouble(0.2, 0.6);
             		double r3 = PseudoRandom.randDouble(0.2, 0.6);
             		if(PseudoRandom.randDouble(0, 1) < 0.5) {
             			r1 = 0;
             		}
             		if(PseudoRandom.randDouble(0, 1) < 0.5) {
             			r2 = 0;
             		}
                     value = xCurrent.getValue(j) + 0.5*(-xCLearning[j] + xCurrent.getValue(j)) 
                     		+ 0.5*(-xDLearning[j] + xCurrent.getValue(j))
                     		+ 0.5*(xParent0.getValue(j) - xParent1.getValue(j));
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
    
    public SolutionSet reproduction1(SolutionSet CS, SolutionSet DS) throws JMException, ClassNotFoundException {
        int populationSize = population.size();
        double[] xCLearning,xDLearning;
        SolutionSet offspringPopulation = new SolutionSet(populationSize);
        Solution[] parents = new Solution[2];
        Solution offSpring;
        for (int i = 0; i < populationSize; i++) {
            offSpring = new Solution(population.get(i));
            XReal xChild = new XReal(offSpring);
            XReal xCurrent = new XReal(population.get(i));

    		int crd1 = PseudoRandom.randInt(0, CS.size()-1);
            int crd2 = PseudoRandom.randInt(0, CS.size()-1);
            while(crd2 == crd1) {
            	crd2 = PseudoRandom.randInt(0, CS.size()-1);
            }
            XReal xCP1 = new XReal(CS.get(crd1));
    		XReal xCP2 = new XReal(CS.get(crd2));
    		
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
            XReal xParent0 = new XReal(parents[0]);
            XReal xParent1 = new XReal(parents[1]);
    		
    		int drd1 = PseudoRandom.randInt(0, DS.size()-1);
            int drd2 = PseudoRandom.randInt(0, DS.size()-1);
            while(drd2 == drd1) {
            	drd2 = PseudoRandom.randInt(0, DS.size()-1);
            }
            XReal xDP1 = new XReal(DS.get(drd1));
    		XReal xDP2 = new XReal(DS.get(drd2));
            
            xCLearning = learningViaDAE(xChild,cMLP);
   		 	xDLearning = learningViaDAE(xChild,dMLP);
   		 
    		for (int j = 0; j < numVars; j++) {
        		double value = PseudoRandom.randDouble(0, 1);
        		if(PseudoRandom.randDouble(0, 1) < 0.5) {
        			value = xCurrent.getValue(j) + 0.5*(xCP1.getValue(j) - xCP2.getValue(j))
        					+ 0.5*(xDP1.getValue(j) - xDP2.getValue(j));
        		}else {
        			value = xCurrent.getValue(j) + 0.25*(-xDLearning[j] + xCurrent.getValue(j))
        					+ 0.25*(-xCLearning[j] + xCurrent.getValue(j))
                     		+ 0.25*(xParent0.getValue(j) - xParent1.getValue(j));
        		}
        		
        		if (value < lowBounds[j])
                    value = lowBounds[j];
                if (value > upBounds[j])
                    value = upBounds[j];
                xChild.setValue(j, value);
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
            double disToNadir = 0.0;
            for(int i=0;i<numObjs;i++){
                value = (solSet.get(p).getObjective(i)-zIdeal[i])/
                        (zNadir[i]-zIdeal[i]);
                solSet.get(p).setNormalizedObjective(i, value);
                sum += value;
                disToIdeal += value*value;
                disToNadir += (value-1)*(value-1);
            }
            
            disToIdeal = Math.sqrt(disToIdeal);
            disToNadir = Math.sqrt(disToNadir);
            
            solSet.get(p).setSumValue(sum);
            //solSet.get(p).setFitness(sum-numObjs);
            solSet.get(p).setFitness(disToNadir);
            solSet.get(p).setDistanceToIdealPoint(disToIdeal);
            for(int i=0;i<numObjs;i++){
            	solSet.get(p).setIthTranslatedObjective(i, solSet.get(p).getNormalizedObjective(i)/sum);
            }
        }
    }

    public SolutionSet getConLabel(SolutionSet st0, SolutionSet st1){
        int size = st1.size();
        SolutionSet labelSet = new SolutionSet(size);
        for(int i=0;i<size;i++) {
        	int rdInt = PseudoRandom.randInt(0, st0.size() - 1);
        	if(PseudoRandom.randDouble(0, 1) < 0.0) {
        		labelSet.add(st0.get(rdInt));
        	}else {
        		Solution sol = st1.get(i);
                double minDis = distance(sol, st0.get(0));
                int minIndex = 0;
                for(int j=1;j<st0.size();j++) {
                    double dis = distance(sol, st0.get(j));
                    if(dis < minDis) {
                        minDis = dis;
                        minIndex = j;
                    }
                }
                labelSet.add(st0.get(minIndex));
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

    public SolutionSet[] getInvertedDivClassifiedSet(){
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
				list.get(i).sort(new FitnessComparator());
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
    
    public SolutionSet[] getConClassifiedSet(){
        SolutionSet[] st = new SolutionSet[3];
        st[0] = new SolutionSet();
        st[1] = new SolutionSet();
        NondominatedRanking ranking = new NondominatedRanking(population);
        int remain = population.size()/2;
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
        st[2] = getConLabel(st[0], st[1]);
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
				list.get(i).sort(new SumValueComparator(true));
				list.get(i).get(size-1).setindex(0);
				st[0].add(list.get(i).get(size-1));
				for(int j=0; j<size-1; j++) {
					st[1].add(list.get(i).get(j));
					list.get(i).get(j).setindex(1);
					//int rd = PseudoRandom.randInt(j+1, size - 1);
					//st[2].add(list.get(i).get(rd));
					int rd = PseudoRandom.randInt(0, list.size() - 1);
					if(rd == i) {rd = PseudoRandom.randInt(0, list.size() - 1);}
					int rnd = PseudoRandom.randInt(0, list.get(rd).size() - 1);
					st[2].add(list.get(rd).get(rnd));
				}
			}else {
				list.get(i).get(0).setindex(0);
				st[0].add(list.get(i).get(0));
			}
		}
        return st;
    }

}
