package jmetal.metaheuristics.almoea;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

//  AMOEAD.java
//
//  Authors:
//       Songbai Liu <songbai209@qq.com>
//       Jun Li <lijun@szu.edu.cn>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.


import jmetal.core.*;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class AMOEADV4 extends Algorithm {

    private int populationSize_;

    private int H_;
    /**
     * Stores the population
     */
    private SolutionSet population;
    /**
     * Z vector (ideal point)
     */
    double[] z_;
    double[] N_;
    /**
     * Lambda vectors
     */
    // Vector<Vector<Double>> lambda_ ;
    double[][] lambda_;

    int evaluations_;
    
    static double learningRate = 0.1;
    static int K = 10;
    static int numLayers = 3;
    private MLPModel MLP;
    
    double[] ur_;//update rate of the environmental selection during the evolutionary process
    double[] nr_;//Nodominated ratio of the population during the evolutionary process
    
    double[] upBounds;
    double[] lowBounds;
    int numVars;

    /**
     * Constructor
     *
     * @param problem
     *            Problem to solve
     */
    public AMOEADV4(Problem problem) {
        super(problem);
    } // MOEAD

    public SolutionSet execute() throws JMException, ClassNotFoundException {
        int maxEvaluations;
        
        MLP = new MLPModel(problem_.getNumberOfVariables(), K, numLayers, learningRate);

        evaluations_ = 0;
        maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations"))
                .intValue();
        populationSize_ = ((Integer) this.getInputParameter("populationSize"))
                .intValue();

        population = new SolutionSet(populationSize_);

        z_ = new double[problem_.getNumberOfObjectives()];
        N_ = new double[problem_.getNumberOfObjectives()];
        lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];
        /* generate two-layer weight vectors */
		VectorGenerator vg;
		if(problem_.getNumberOfObjectives() == 2) {
			vg = new TwoLevelWeightVectorGenerator(99, 0,problem_.getNumberOfObjectives());
		}else {
			vg = new TwoLevelWeightVectorGenerator(16, 0,problem_.getNumberOfObjectives());
		}
		lambda_ = vg.getVectors();

        // STEP 1. Initialization
        initPopulation();
        initIdealAndNadirPoint();
        numVars = problem_.getNumberOfVariables();
        int sss = maxEvaluations/populationSize_;
        ur_ = new double[sss-1];
        nr_ = new double[sss-1];
        upBounds = new double[numVars];
        lowBounds = new double[numVars];
        for(int var=0;var<numVars;var++) {
            upBounds[var] = problem_.getUpperLimit(var);
            lowBounds[var] = problem_.getLowerLimit(var);
        }
        
        int t = 0;
        // STEP 2. Update
        do {
            SolutionSet offspringPopulation;
            NondominatedRanking rank = new NondominatedRanking(population);
            nr_[t] = (double)rank.getSubfront(0).size()/population.size();
            
            if(evaluations_ < 0.15*maxEvaluations){
            	//InnovizedRepair crossover = new InnovizedRepair(population,problem_);
            	AGSCompetitiveCrossover crossover = new AGSCompetitiveCrossover(population,problem_,MLP);
                offspringPopulation = crossover.doCrossover();
                MLP = crossover.getCurrentMLPModel();
            }else{
            	if(nr_[t] < 0.2) {
            		CompetitiveDECrossover crossover = new CompetitiveDECrossover(population,problem_);
                    offspringPopulation = crossover.doCrossover();
                }else if (nr_[t] > 0.2 && ur_[t] < 0.1){
                    offspringPopulation = samplingNearbyParetoSet();
                }else {
                	CompetitiveDECrossover1 crossover = new CompetitiveDECrossover1(population,problem_);
                    offspringPopulation = crossover.doCrossover();
                }
            	
            }
            SolutionSet union = ((SolutionSet) population).union(offspringPopulation);
            SolutionSet st = getStSolutionSet(union,(int)(populationSize_*(1.0 + 0.5*evaluations_/maxEvaluations)));
            estimateIdealNadirPoint(st);
            normalization(st);
            population.clear();
            population = new DecompositionBasedSelection(st, problem_.getNumberOfObjectives(), lambda_).environmentalSelection();
            
            double sum = 0;
            for(int i=0;i<population.size();i++) {
            	sum += population.get(i).getindex();
            	population.get(i).setindex(0);
            }
            ur_[t] = sum/populationSize_;
            t++;
            
            evaluations_ += populationSize_;
            
        } while (evaluations_ < maxEvaluations);

        //printMetric("UR_"+ problem_.getNumberOfObjectives()+"Obj_"+problem_.getName()+
                //"_" + numVars + "D" +".txt", ur_);
        //printMetric("NR_"+ problem_.getNumberOfObjectives()+"Obj_"+problem_.getName()+
                //"_" + numVars + "D" +".txt", nr_);
        NondominatedRanking ranking = new NondominatedRanking(population);

        return ranking.getSubfront(0);
    }
    
    public static void printMetric(String path,double[] indicator){
        try {
            /* Open the file */
            FileOutputStream fos   = new FileOutputStream(path)     ;
            OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
            BufferedWriter bw      = new BufferedWriter(osw)        ;
            for (int i = 0; i < indicator.length; i++) {
                bw.write(indicator[i]+" ");
                bw.newLine();
            }

            /* Close the file */
            bw.close();
        }catch (IOException e) {
            Configuration.logger_.severe("Error acceding to the file");
            e.printStackTrace();
        }
    } // printGD

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

    public void initPopulation() throws JMException, ClassNotFoundException {
        for (int i = 0; i < populationSize_; i++) {
            Solution newSolution = new Solution(problem_);
            problem_.evaluate(newSolution);
            evaluations_++;
            newSolution.setindex(0);
            population.add(newSolution);
        } // for
       
    } // initPopulation
    
    
    public void resamplingPopulation() throws JMException, ClassNotFoundException {
    	SolutionSet newSet = new SolutionSet(populationSize_);
        for (int i = 0; i < populationSize_; i++) {
            Solution newSolution = new Solution(problem_);
			/*
			 * XReal xChild = new XReal(newSolution); for (int j = 0; j < numVars; j++) {
			 * double value; value=PseudoRandom.randDouble((lowBounds[j] -
			 * xChild.getValue(j))/2, (upBounds[j] - xChild.getValue(j))/2); value =
			 * xChild.getValue(j) + value; xChild.setValue(j, value); }
			 */
            problem_.evaluate(newSolution);
            //evaluations_++;
            newSolution.setindex(0);
            newSet.add(newSolution);
        } // for
        
        SolutionSet union = ((SolutionSet) population).union(newSet);
        //SolutionSet st = getStSolutionSet(union,(int)(populationSize_*1.5));
        estimateIdealNadirPoint(union);
        normalization(union);
        population.clear();
        population = new DecompositionBasedSelection(union, problem_.getNumberOfObjectives(), lambda_).environmentalSelection();
       
    } // initPopulation
    
    public SolutionSet samplingNearbyParetoSet() throws JMException {
    	SolutionSet newSet = new SolutionSet(populationSize_);
    	NondominatedRanking rank = new NondominatedRanking(population);
    	SolutionSet front = rank.getSubfront(0);
    	int size = front.size();
    	
    	for(int i=0;i<populationSize_;i++) {
    		int index = PseudoRandom.randInt(0, size-1);
    		Solution newSolution = new Solution(front.get(index));
    		doMutation(0.001, newSolution);
    		problem_.evaluate(newSolution);
    		newSolution.setindex(1);
            newSet.add(newSolution);
    	}

    	return newSet;
    }

    /**
     *
     */
    void initIdealAndNadirPoint() throws JMException, ClassNotFoundException {
        for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
            z_[i] = 1.0e+30;
            N_[i] = -1.0e+30;
        } // for
    } // initIdealPoint
    
    public void estimateIdealNadirPoint(SolutionSet solSet){
        for(int p=0;p<solSet.size();p++){
            for(int i=0;i<problem_.getNumberOfObjectives();i++){
                if(solSet.get(p).getObjective(i) < z_[i])
                    z_[i] = solSet.get(p).getObjective(i);
                if(solSet.get(p).getObjective(i) > N_[i])
                	N_[i] = solSet.get(p).getObjective(i);
            }
        }
    }

    public void normalization(SolutionSet solSet){
        double value;
        for(int p=0;p<solSet.size();p++){
            double sum = 0.0;
            for(int i=0;i< problem_.getNumberOfObjectives();i++){
                value = (solSet.get(p).getObjective(i)-z_[i])/ (N_[i]-z_[i]);
                //value = (solSet.get(p).getObjective(i)-z_[i]);
                solSet.get(p).setNormalizedObjective(i, value);
                sum += value;
            }
            solSet.get(p).setSumValue(sum);
        }
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

} // AMOEAD

