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
import jmetal.problems.LSMOP.*;
import jmetal.problems.ProblemFactory;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.*;
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
 * This class executes the algorithm described in: H. Li and Q. Zhang,
 * "Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D
 * and NSGA-II". IEEE Trans on Evolutionary Computation, vol. 12, no 2, pp
 * 284-302, April/2009.
 */
public class MOEAD_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

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
		int n =10000;
		int m = 3;
		int num_samples = 5000;
		if(m == 3) {
			num_samples = 10000;
		}


		//double[] GDarray=new double[runtimes];

		//double[] spreadarray=new double[runtimes];
		//double[] Hypervolume=new double[runtimes];

		
		//for(int i=0;i<runtimes;i++){
		Problem problem=null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator

		QualityIndicator indicators=null; // Object to get quality indicators

		HashMap parameters; // Operator parameters

		indicators = null;
		for(int fun=1;fun<=9;fun++){
			int runtimes=20;
			double[] IGDarray=new double[runtimes];
			long Execuion_time=0;
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
		} // else
			m = problem.getNumberOfObjectives();
			for(int i=0;i<runtimes;i++){
		algorithm = new MOEAD(problem);
		// algorithm = new MOEAD_DRA(problem);
				if(m == 2) {
					algorithm.setInputParameter("populationSize", 100);
					algorithm.setInputParameter("maxEvaluations", 100*100);
					//algorithm.setInputParameter("div1", 49);
					//algorithm.setInputParameter("div2", 0);
					//algorithm.setInputParameter("div3", 99);
				//	algorithm.setInputParameter("div4", 0);
				}else if(m == 3) {
					algorithm.setInputParameter("populationSize", 153);
					algorithm.setInputParameter("maxEvaluations", 100*153);
				//	algorithm.setInputParameter("div1", 11);
					//algorithm.setInputParameter("div2", 0);
					//algorithm.setInputParameter("div3", 16);
				//	algorithm.setInputParameter("div4", 0);
				}
		// Algorithm parameters

		// Directory with the files containing the weight vectors used in
		// Q. Zhang, W. Liu, and H Li, The Performance of a New Version of
		// MOEA/D
		// on CEC09 Unconstrained MOP Test Instances Working Report CES-491,
		// School
		// of CS & EE, University of Essex, 02/2009.
		// http://dces.essex.ac.uk/staff/qzhang/MOEAcompetition/CEC09final/code/ZhangMOEADcode/moead0305.rar
	//	algorithm.setInputParameter("dataDirectory",
		//		"E:\\sbMaOP\\SbMaOP\\Weights");

		//algorithm.setInputParameter("finalSize", 200); // used by MOEAD_DRA

		algorithm.setInputParameter("T", 20);
		algorithm.setInputParameter("delta", 0.9);
		algorithm.setInputParameter("nr", 2);

		// Crossover operator
		parameters = new HashMap();
		parameters.put("CR", 1.0);
		parameters.put("F", 0.5);
		crossover = CrossoverFactory.getCrossoverOperator(
				"DifferentialEvolutionCrossover", parameters);

		// Mutation operator
		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation",
				parameters);

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);

		// Execute the Algorithm
		long initTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		//long estimatedTime = System.currentTimeMillis() - initTime;

		Execuion_time+=(System.currentTimeMillis() - initTime);
		population.printObjectivesToFile("MOEAD_DE_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1)+".txt" );
		IGDarray[i]=indicators.getIGD1(population);
	    //spreadarray[i]=indicators.getSpread(population);
	   // Hypervolume[i]=indicators.getHypervolume(population);
		}
		//printGD("FUN"+fun+"IGD",IGDarray);
		//double sumGD=0;
		double sumIGD=0;
		//double sumSP=0;
		//double sumHV = 0;
		for(int i=0;i<runtimes;i++){
		 // sumGD+=GDarray[i];
		  sumIGD+=IGDarray[i];
		 // sumSP+=spreadarray[i];
		  //sumHV+=Hypervolume[i];
		}
		logger_.info("Total execution time: " + Execuion_time + "ms");
		//System.out.println("avrGD-fun"+fun+" = "+sumGD/runtimes);
		//System.out.println("avrSP-fun"+fun+" = "+sumSP/runtimes);
			System.out.println("IGD_"+ problem.getNumberOfObjectives()+"Obj_"+problem.getName()+
					"_" + problem.getNumberOfVariables() + "D" + " = "+sumIGD/runtimes);
		//System.out.println("avrHV-fun"+fun+" = "+sumHV/runtimes);
		}
	} // main
} // MOEAD_main
