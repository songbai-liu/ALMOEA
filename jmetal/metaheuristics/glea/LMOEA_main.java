//  LMPSO_main.java
//
// Author:
//     Songbai Liu <Songbai209@qq.com>
// Copyright (c) 2019 Songbai Liu
//
// This Program is free software: you can redistribute it and/or modify 
// it under the terms of the GNU Lesser General Public License as published
// by the free software foundation, either version 4 of license, or any 
// later version (at your option).

package jmetal.metaheuristics.glea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.nsgaII.NSGAII_ANN;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.DTLZ.*;
import jmetal.problems.IMF.*;
import jmetal.problems.LMF.*;
import jmetal.problems.LSMOP.*;
import jmetal.problems.M2M.*;

import jmetal.problems.WFG.*;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

/**
 * Class to configure and execute the NSGA-II algorithm.
 * 
 * Besides the classic NSGA-II, a steady-state version (ssNSGAII) is also
 * included (See: J.J. Durillo, A.J. Nebro, F. Luna and E. Alba "On the Effect
 * of the Steady-State Selection Scheme in Multi-Objective Genetic Algorithms"
 * 5th International Conference, EMO 2009, pp: 183-197. April 2009)
 */

public class LMOEA_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**

	 *            Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName
	 *             paretoFrontFile
	 */
	public static void printGD(String path,double[] GD){
	    try {
	    	/* Open the file */
	    	FileOutputStream fos   = new FileOutputStream(path)     ;
	    	OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
	    	BufferedWriter bw      = new BufferedWriter(osw)        ;             
	    	for (int i = 0; i < GD.length; i++) {  
	    		bw.write(GD[i]+" ");
	    		bw.newLine();     
	    	}
	    	/* Close the file */
	    	bw.close();
	    }catch (IOException e) {
	    	Configuration.logger_.severe("Error acceding to the file");
	    	e.printStackTrace();
	    }       
	} // print
	
	public static void printave(String path,double aveIGD,double varianceIGD,double aveHypervolume,double varianceHV){
	    try {
	    	/* Open the file */
	    	FileOutputStream fos   = new FileOutputStream(path) ;
	    	OutputStreamWriter osw = new OutputStreamWriter(fos) ;
	    	BufferedWriter bw      = new BufferedWriter(osw)  ;            
	        bw.write(aveIGD+" ");
	        bw.newLine(); 
	        bw.write(varianceIGD+" ");
	        bw.newLine();
	        bw.write(aveHypervolume+" ");
	        bw.newLine();  
	        bw.write(varianceHV+" ");
	        bw.newLine();
	        /* Close the file */
		    bw.close();
		 }catch (IOException e) {
		      Configuration.logger_.severe("Error acceding to the file");
		      e.printStackTrace();
		 }       
	} // print
	
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAII_main.log");
		logger_.addHandler(fileHandler_);
		int n = 3000;
		int m = 2;
		int num_samples = 5000;
		if(m == 3) {
			num_samples = 10000;
		}
		int tp1 = 2;//type of formulation model, 0 for addition, 1 for multiplication, and 2 for mixed/hybrid
		int tp2 = 2;//type of variable linkage, 0 for linear linkage, 1 for nonlinear linkage, and 2 for mixed/hybrid
		int tp3 = 1;//type of deep grouping on variables, 0 for even grouping, 1 for nonuniform grouping
		int tp4 = 1;//type of contribution on variables, 0 for balanced contribution, 1 for unbalanced contribution
		Problem problem=null; // The problem to solve
	    Algorithm algorithm; // The algorithm to use
	    Operator crossover; // Crossover operator
	    Operator mutation; // Mutation operator
	    Operator selection; // Selection operator
	    HashMap parameters; // Operator parameters
	    QualityIndicator indicators; // Object to get quality indicators
	    indicators = null;
		for(int fun=31;fun<=38;fun++){
			int runtimes=20;
			double[] GDarray=new double[runtimes];
			double[] IGDarray=new double[runtimes];
			double[] spreadarray=new double[runtimes];
			double[] Hypervolume=new double[runtimes];
			long Execution_time=0;
		    if (args.length == 1) {
		    	Object[] params = { "Real" };
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    } // if
		    else if (args.length == 2) {
		    	Object[] params = { "Real" };
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    	indicators = new QualityIndicator(problem, args[1]);
		    } // if
		    else { // Default problem
				if(fun == 31) {
					problem = new LMF1("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "concave_" + m +"D" +".pf");
				}
				if(fun == 32) {
					problem = new LMF2("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "convex_" + m +"D" +".pf");
				}
				if(fun == 33) {
					problem = new LMF3("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "linear_" + m +"D" +".pf");
				}
				if(fun == 34) {
					problem = new LMF4("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "concave_" + m +"D" +".pf");
				}
				if(fun == 35) {
					problem = new LMF5("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "convex_" + m +"D" +".pf");
				}
				if(fun == 36) {
					problem = new LMF6("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "linear_" + m +"D" +".pf");
				}
				if(fun == 37) {
					problem = new LMF7("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "concave_" + m +"D" +".pf");
				}
				if(fun == 38) {
					problem = new LMF8("Real",n,m);
					indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "linear_" + m +"D" +".pf");
				}
		    } // else
		   for( int i=0;i<runtimes;i++)
		   {

		   		 algorithm= new GLEA_new(problem);
		    	// Algorithm parameters
		    	if(m == 2) {
	    			algorithm.setInputParameter("populationSize", 100);
		    		algorithm.setInputParameter("maxGenerations", 100);
		    		algorithm.setInputParameter("div1", 99);
		    		algorithm.setInputParameter("div2", 0);
		    		algorithm.setInputParameter("div3", 49);
		    		algorithm.setInputParameter("div4", 0);
	    		}else if(m == 3) {
	    			algorithm.setInputParameter("populationSize", 153);
		    		algorithm.setInputParameter("maxGenerations", 100);
		    		algorithm.setInputParameter("div1", 16);
		    		algorithm.setInputParameter("div2", 0);
		    		algorithm.setInputParameter("div3", 10);
		    		algorithm.setInputParameter("div4", 0);
	    		}
		    	// Mutation and Crossover for Real codification
		    	// Crossover operator		
		    	parameters = new HashMap();
		    	parameters.put("probability", 1.0);
		    	parameters.put("distributionIndex", 20.0);
		    	crossover = CrossoverFactory.getCrossoverOperator(
		    			"SBXCrossover", parameters);
		    	parameters = new HashMap();
		    	parameters.put("probability", 1.0/problem.getNumberOfVariables());
		    	parameters.put("distributionIndex", 20.0);
		    	mutation = MutationFactory.getMutationOperator("PolynomialMutation",
		    			parameters);

		    	// Selection Operator
		    	parameters = null;
		    	selection = SelectionFactory.getSelectionOperator("BinaryTournament",
		    			parameters);
		    	// Add the operators to the algorithm
		    	algorithm.addOperator("crossover", crossover);
		    	algorithm.addOperator("mutation", mutation);
		    	algorithm.addOperator("selection", selection);
		    	// Add the indicator object to the algorithm
		    	algorithm.setInputParameter("indicators", indicators);
		    	// Execute the Algorithm
		    	long initTime = System.currentTimeMillis();
		    	SolutionSet population = algorithm.execute();
		    	Execution_time+=(System.currentTimeMillis() - initTime);
		    	// Result messages
		    	population.printObjectivesToFile("GLEA_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1)+".txt" );
		    	IGDarray[i]=indicators.getIGD1(population);
		    }
		    
		    //printGD("GLEAP6_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + n + "D_IGD" + tp1 + tp2 + tp3 + tp4 + ".txt",IGDarray);
		    double sumIGD=0;
		    for(int i=0;i<runtimes;i++){
		    	sumIGD+=IGDarray[i];
		    }	  	  
		    logger_.info("Total execution time: " + Execution_time + "ms");
		    System.out.println("IGD_"+ problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ 
				  "_" + problem.getNumberOfVariables() + "D" + " = "+sumIGD/runtimes);
		} //main
	}
} // NSGAII_main
