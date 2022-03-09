//  MOEAD_main.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.moead;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.problems.Schaffer;
import jmetal.problems.LMF.LMF1;
import jmetal.problems.LMF.LMF2;
import jmetal.problems.LMF.LMF3;
import jmetal.problems.LMF.LMF4;
import jmetal.problems.LMF.LMF5;
import jmetal.problems.LMF.LMF6;
import jmetal.problems.LMF.LMF7;
import jmetal.problems.LMF.LMF8;
import jmetal.problems.LSMOP.LSMOP1;
import jmetal.problems.LSMOP.LSMOP2;
import jmetal.problems.LSMOP.LSMOP3;
import jmetal.problems.LSMOP.LSMOP4;
import jmetal.problems.LSMOP.LSMOP5;
import jmetal.problems.LSMOP.LSMOP6;
import jmetal.problems.LSMOP.LSMOP7;
import jmetal.problems.LSMOP.LSMOP8;
import jmetal.problems.LSMOP.LSMOP9;
import jmetal.problems.cec2009Competition.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
import jmetal.qualityIndicator.util.MetricsUtil;
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
 * This class executes the algorithm described in: H. Li and Q. Zhang,
 * "Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D
 * and NSGA-II". IEEE Trans on Evolutionary Computation, vol. 12, no 2, pp
 * 284-302, April/2009.
 */
public class computeIGD_main_2obj {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object
	
	public static jmetal.qualityIndicator.util.MetricsUtil utils_;

	/**
	// * @param args
	 *            Command line arguments. The first (optional) argument
	 *            specifies the problem to solve.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options - jmetal.metaheuristics.moead.MOEAD_main
	 *             - jmetal.metaheuristics.moead.MOEAD_main problemName -
	 *             jmetal.metaheuristics.moead.MOEAD_main problemName
	 *             ParetoFrontFile
	 * @throws ClassNotFoundException
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
	  } // printGD
	
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		
		// Logger object and file to store log messages
				logger_ = Configuration.logger_;
				fileHandler_ = new FileHandler("MOEAD.log");
				logger_.addHandler(fileHandler_);
		int n = 5000;
		int m = 2;
		int runtimes=10;
		for(int fun=31;fun<=38;fun++){
		double[] IGDarray=new double[runtimes];
		double[] HVarray=new double[runtimes];
		long Execuion_time=0;
		
		Problem problem=null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator

		QualityIndicator indicators; // Object to get quality indicators

		HashMap parameters; // Operator parameters

		indicators = null;
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
			if(fun==1){
			      problem = new LSMOP1("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP1_2D_5000.txt");
			    	}
			if(fun==2){
			      problem = new LSMOP2("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP2_2D_5000.txt");
			    	}
			if(fun==3){
			      problem = new LSMOP3("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP3_2D_5000.txt");
			    	}
			if(fun==4){
			      problem = new LSMOP4("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP4_2D_5000.txt");
			    	}
			if(fun==5){
			      problem = new LSMOP5("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP5_2D_5000.txt");
			    	}
			if(fun==6){
			      problem = new LSMOP6("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP6_2D_5000.txt");
			    	}
			if(fun==7){
			      problem = new LSMOP7("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP7_2D_5000.txt");
			    	}
			if(fun==8){
			      problem = new LSMOP8("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP8_2D_5000.txt");
			    	}
			if(fun==9){
			      problem = new LSMOP9("Real",n,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP9_2D_5000.txt");
			    	}
			if(fun==10){
			      problem = new UF1("Real",n);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF1_500.txt");
			    	}
			if(fun==11){
			      problem = new UF2("Real",n);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF2_500.txt");
			    	}
			if(fun==12){
			      problem = new UF3("Real",n);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF3_500.txt");
			    	}
			if(fun==13){
			      problem = new UF4("Real",n);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF4_500.txt");
			    	}
			if(fun==14){
			      problem = new UF5("Real",n,10,0.1);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF5_21.txt");
			    	}
			if(fun==15){
			      problem = new UF6("Real",n,2,0.1);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF6_668.txt");
			    	}
			if(fun==16){
			      problem = new UF7("Real",n);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF7_500.txt");
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
		  
		} // else
		utils_ = new MetricsUtil();
		String[] algorithms = new String[22];
		algorithms[0] = "LMOCSO";
		algorithms[1] = "LSMOF";
		algorithms[2] = "WOFNSGA";
		algorithms[3] = "MOEADVA";
		algorithms[4] = "WOF";
		algorithms[5] = "DGEA";
		algorithms[6] = "NSGAII";
		algorithms[7] = "NSGAII_ANN";
		algorithms[8] = "SMSEMOA";
		algorithms[9] = "SMSEMOA_LDE";
		algorithms[10] = "MOEAD_DE";
		algorithms[11] = "GLEA";
		algorithms[12] = "AMOEADV2";
		algorithms[13] = "AMOEADV3";
		algorithms[14] = "AMOEADIn";

		for(int a=14;a<=14;a++){
			for(int i=0;i<runtimes;i++){//D:\ResearchWork\2021\2021_P8\Review1\AMOEAD
			   String path = "D:\\ResearchWork\\2021\\2021_P8\\Review1\\AMOEAD\\" + algorithms[a];
			//   path += "\\learing-DE\\m2\\"+algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName() + "_" + n +"D_run"+(i+1)+".txt";
				//String path = "D:\\2021_P1\\population\\FE_10000\\"+ algorithms[a];
				//path +="\\"+ algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName() + "_" + n +"D_run"+(i+1)+".txt";
			//	path +="\\"+ algorithms[a]+"_"+problem.getName()+"_M"+problem.getNumberOfObjectives()+"_D"+problem.getNumberOfVariables()+"_"+(i+1)+".txt";
			 path += "\\"+algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName() + "_" + n +"D_run"+(i+1)+".txt";
			   double[][] population = utils_.readFront(path); 
		       IGDarray[i] = indicators.getIGD1(population);
			}
		//	printGD("D:\\2021_P1\\IGD\\"+algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + n + "D_IGD.txt",IGDarray);
			printGD(algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+ problem.getName() + "_" + n + "D_IGD.txt",IGDarray);
			double sumIGD = 0;
			for(int i=0;i<runtimes;i++){
			  sumIGD+=IGDarray[i];
			}
			//logger_.info("Total execution time: " + Execuion_time + "ms");
			System.out.println("avrHV-fun"+fun+" = "+sumIGD/runtimes);
		}
	}//for
	} // main
} // MOEAD_main
