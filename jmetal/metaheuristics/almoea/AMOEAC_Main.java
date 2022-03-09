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

package jmetal.metaheuristics.almoea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.IMF.*;
import jmetal.problems.LMF.*;
import jmetal.problems.LSMOP.*;
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

public class AMOEAC_Main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	//* @param args
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
	    	FileOutputStream fos   = new FileOutputStream(path);//java�?�囦欢�?�撳嚭娴侊�?�?�涘缓�?�囦欢娴�
	    	OutputStreamWriter osw = new OutputStreamWriter(fos);//OutputStreamWriter�?�瓧绗︽祦閫氬悜瀛楄妭娴佺殑妗ユ 
	    	BufferedWriter bw      = new BufferedWriter(osw);//缂撳啿�?��               
	    	for (int i = 0; i < GD.length; i++) {  
	    		bw.write(GD[i]+" ");//�??�?埌缂撳啿�?��
	    		bw.newLine(); //鎹㈣       
	    	}
	    	/*Close the file*/
	    	bw.close();
	    }catch (IOException e) {
	    	Configuration.logger_.severe("Error acceding to the file");
	    	e.printStackTrace();
	    }       
	} // printGD
	
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
	} // printave
	
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAII_main.log");
		logger_.addHandler(fileHandler_);
		int n = 1000;
		int m = 2;
		int num_samples = 5000;
		if(m == 3) {
			num_samples = 10000;
		}
		Problem problem=null; // The problem to solve
	    Algorithm algorithm; // The algorithm to use
	    Operator crossover; // Crossover operator
	    Operator mutation; // Mutation operator
	    Operator selection; // Selection operator
	    HashMap parameters; // Operator parameters
	    QualityIndicator indicators; // Object to get quality indicators
	    indicators = null;
		for(int fun=31;fun<=38;fun++){
			int runtimes=1;
			double[] GDarray=new double[runtimes];
			double[] IGDarray=new double[runtimes];
			double[] spreadarray=new double[runtimes];
			double[] Hypervolume=new double[runtimes];
			long Execution_time=0;
		    if (args.length == 1) {
		    	Object[] params = {"Real"};
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    } // if
		    else if (args.length == 2) {
		    	Object[] params = { "Real" };
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    	indicators = new QualityIndicator(problem, args[1]);
		    } // if
		    else { // Default problem
		    	if(fun==1){
		    		problem = new LSMOP1("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==2){
		    		problem = new LSMOP2("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==3){
		    		problem = new LSMOP3("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==4){
		    		problem = new LSMOP4("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==5){
		    		problem = new LSMOP5("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==6){
		    		problem = new LSMOP6("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==7){
		    		problem = new LSMOP7("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==8){
		    		problem = new LSMOP8("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==9){
		    		problem = new LSMOP9("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
			    }
		    	if(fun==10) {
		    		m=2;
		    		problem = new IMF1("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun==11) {
		    		m=2;
		    		problem = new IMF2("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun==12) {
		    		m=2;
		    		problem = new IMF3("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun==13) {
		    		m=3;
		    		problem = new IMF4("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + 10000 +".txt");
		    	}
		    	if(fun==14) {
		    		m=2;
		    		problem = new IMF5("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun==15) {
		    		m=2;
		    		problem = new IMF6("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun==16) {
		    		m=2;
		    		problem = new IMF7("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun==17) {
		    		m=3;
		    		problem = new IMF8("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + 10000 +".txt");
		    	}
		    	if(fun==18) {
		    		m=2;
		    		problem = new IMF9("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun==19) {
		    		m=2;
		    		problem = new IMF10("Real",n,m);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\" + problem.getName() + "_" + m +"D_" + num_samples +".txt");
		    	}
		    	if(fun == 20) {
		    		problem = new UF1("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\" + problem.getName() + "_" + m +"D_" + 1000 +".txt");
		    	}
		    	if(fun == 21) {
		    		problem = new UF2("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\" + problem.getName() + "_" + m +"D_" + 1000 +".txt");
		    	}
		    	if(fun == 22) {
		    		problem = new UF3("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\" + problem.getName() + "_" + m +"D_" + 1000 +".txt");
		    	}
		    	if(fun == 23) {
		    		problem = new UF4("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\" + problem.getName() + "_" + m +"D_" + 1000 +".txt");
		    	}
		    	if(fun == 24) {
		    		problem = new UF5("Real",n, 10, 0.1);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\" + problem.getName() + "_" + m +"D_" + 21 +".txt");
		    	}
		    	if(fun == 25) {
		    		problem = new UF6("Real",n, 2, 0.1);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\" + problem.getName() + "_" + m +"D_" + 668 +".txt");
		    	}
		    	if(fun == 26) {
		    		problem = new UF7("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\" + problem.getName() + "_" + m +"D_" + 1000 +".txt");
		    	}
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
		    	if(fun == 39) {
		    		problem = new LMF9("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "inverted_concave_" + 3 +"D" +".pf");
		    	}
		    	if(fun == 40) {
		    		problem = new LMF10("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "inverted_linear_" + 3 +"D" +".pf");
		    	}
		    	if(fun == 41) {
		    		problem = new LMF11("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "inverted_concave_" + 3 +"D" +".pf");
		    	}
		    	if(fun == 42) {
		    		problem = new LMF12("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "inverted_linear_" + 3 +"D" +".pf");
		    	}
		    } // else
		    m = problem.getNumberOfObjectives();
		    for(int i=0;i<runtimes;i++){
		    	algorithm = new AMOEAC(problem);
		    	// Algorithm parameters
		    	if(m == 2) {
		    		algorithm.setInputParameter("populationSize", 100);
		    		algorithm.setInputParameter("maxEvaluations", 100*100);
		    	}else {
		    		algorithm.setInputParameter("populationSize", 153);
		    		algorithm.setInputParameter("maxEvaluations", 153*1000);
		    	}
		    	
		    	// Add the indicator object to the algorithm
		    	algorithm.setInputParameter("indicators", indicators);
		    	// Execute the Algorithm
		    	long initTime = System.currentTimeMillis();
		    	SolutionSet population = algorithm.execute();
		    	Execution_time+=(System.currentTimeMillis() - initTime);
		    	// Result messages
		    	population.printObjectivesToFile("AMOEAC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1)+".txt" );
		    	IGDarray[i]=indicators.getIGD1(population);
		    }
		    double sumIGD=0;
		    for(int i=0;i<runtimes;i++){
		    	sumIGD+=IGDarray[i];
		    }	  	  
		    logger_.info("Total execution time: " + Execution_time + "ms");
		    System.out.println("IGD_"+ problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ 
				  "_" + problem.getNumberOfVariables() + "D" + " = "+sumIGD/runtimes);
		} //main
	}
} // AMOEAC_main
