//MOEAD_main.java
//
//Author:
//   Antonio J. Nebro <antonio@lcc.uma.es>
//   Juan J. Durillo <durillo@lcc.uma.es>
//
//Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU Lesser General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU Lesser General Public License for more details.
//
//You should have received a copy of the GNU Lesser General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
import jmetal.problems.LZ09.*;
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
public class My_MOEAD_main {
public static Logger logger_; // Logger object
public static FileHandler fileHandler_; // FileHandler object

/**
 * @param args
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
	for(int fun=32;fun<=34;fun++){	
	int runtimes=30;
	double[] GDarray=new double[runtimes];
	double[] IGDarray=new double[runtimes];
	double[] spreadarray=new double[runtimes];
	double[] Hypervolume=new double[runtimes];
	long Execuion_time=0;
	
	for(int i=0;i<runtimes;i++){
	Problem problem=null; // The problem to solve
	Algorithm algorithm; // The algorithm to use
	Operator crossover,crossover2; // Crossover operator
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
		// problem = new Kursawe("Real", 3);
		if(fun==1){
	  	      problem = new ZDT1("Real");
	  	      
	  	      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT1_501.txt" ) ;
	  	    	}//problem = new WFG1("Real");
	  	if(fun==2){
	  	      problem = new ZDT2("Real");
	  	      
	  	      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT2_501.txt" ) ;
	  	    	}//problem = new WFG1("Real");
	  	if(fun==3){
	  	      problem = new ZDT3("Real");
	  	      
	  	      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT3_269.txt" ) ;
	  	    	}
	  	if(fun==4){
		      problem = new ZDT4("Real",10);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT4_501.txt" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==5){
		      problem = new ZDT6("Real",10);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT6_774.txt" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==6){
		      problem = new DTLZ1("Real",10,3);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ1.pf" ) ;
		    	}
	  	if(fun==7){
		      problem = new DTLZ2("Real",10,3);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ2.pf" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==8){
		      problem = new DTLZ3("Real",10,3);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ3.pf" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==9){
		      problem = new DTLZ4("Real",10,3);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ4.pf" ) ;
		    	}
	  	if(fun==10){
		      problem = new DTLZ5("Real",10,3);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ5.txt" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==11){
		      problem = new DTLZ6("Real",10,3);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ6.txt" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==12){
		      problem = new DTLZ7("Real",10,3);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ7.pf" ) ;
		    	}
		if(fun==13){
	  	      problem = new WFG1("Real",8,2,2);
	  	      
	  	      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG1_605.txt" ) ;
	  	    	}//problem = new WFG1("Real");
	  	if(fun==14){
	  	      problem = new WFG2("Real",8,2,2);
	  	      
	  	      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG2_111.txt" ) ;
	  	    	}//problem = new WFG1("Real");
	  	if(fun==15){
	  	      problem = new WFG3("Real",8,2,2);
	  	      
	  	      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG3_right.txt" ) ;
	  	    	}
		if(fun==16){
		      problem = new WFG4("Real",8,2,2);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG4_1181.txt" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==17){
		      problem = new WFG5("Real",8,2,2);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG5_694.txt" ) ;
		    	}
	  	if(fun==18){
	  	      problem = new WFG6("Real",8,2,2);
	  	      
	  	      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\wfg6_166.txt" ) ;
	  	    	}//problem = new WFG1("Real");
	  	if(fun==19){
		      problem = new WFG7("Real",8,2,2);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG7_2435.txt" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==20){
		      problem = new WFG8("Real",8,2,2);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG8_201.txt" );
		    	}
		if(fun==21){
		      problem = new WFG9("Real",8,2,2);
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG9_2591.txt" ) ;
		    	}//problem = new WFG1("Real");
	  	if(fun==22){
		      problem = new Fonseca("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\Fonseca.pf" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==23){
		      problem = new Kursawe("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\Kursawe.pf" ) ;
		    	}
		if(fun==24){
		      problem = new Schaffer("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\Schaffer.pf" ) ;
		    	}//problem = new WFG1("Real");
		if(fun==25){
		      problem = new UF1("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF1.dat" ) ;
		    	}
		if(fun==26){
		      problem = new UF2("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF2.dat" ) ;
		    	}
		if(fun==27){
		      problem = new UF3("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF3.dat" ) ;
		    	}
		if(fun==28){
		      problem = new UF4("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF4.dat" ) ;
		    	}
		if(fun==29){
		      problem = new UF5("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF5.dat" ) ;
		    	}
		if(fun==30){
		      problem = new UF6("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF6.dat" ) ;
		    	}
		if(fun==31){
		      problem = new UF7("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF7.dat" ) ;
		    	}
		if(fun==32){
		      problem = new UF8("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF8.DAT" ) ;
		    	}
		if(fun==33){
		      problem = new UF9("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF9.DAT" ) ;
		    	}
		if(fun==34){
		      problem = new UF10("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF10.DAT" ) ;
		    	}
		if(fun==35){
		      problem = new LZ09_F1("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F1.pf" ) ;
		    	}
		if(fun==36){
		      problem = new LZ09_F2("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F2.pf" ) ;
		    	}
		if(fun==37){
		      problem = new LZ09_F3("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F3.pf" ) ;
		    	}
		if(fun==38){
		      problem = new LZ09_F4("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F4.pf" ) ;
		    	}
		if(fun==39){
		      problem = new LZ09_F5("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F5.pf" ) ;
		    	}
		if(fun==40){
		      problem = new LZ09_F6("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F6.pf" ) ;
		    	}
		if(fun==41){
		      problem = new LZ09_F7("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F7.pf" ) ;
		    	}
		if(fun==42){
		      problem = new LZ09_F8("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F8.pf" ) ;
		    	}
		if(fun==43){
		      problem = new LZ09_F9("Real");
		      
		      indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\LZ09_F9.pf" ) ;
		    	}
		// problem = new Kursawe("BinaryReal", 3);
		// problem = new Water("Real");
		// problem = new ZDT1("ArrayReal", 100);
		// problem = new ConstrEx("Real");
		// problem = new DTLZ1("Real");
		// problem = new OKA2("Real") ;
	} // else

	algorithm = new My_MOEAD(problem);
	// algorithm = new MOEAD_DRA(problem);My_MOEAD

	// Algorithm parameters
	algorithm.setInputParameter("populationSize", 1000);
	algorithm.setInputParameter("maxEvaluations", 300000);

	// Directory with the files containing the weight vectors used in
	// Q. Zhang, W. Liu, and H Li, The Performance of a New Version of
	// MOEA/D
	// on CEC09 Unconstrained MOP Test Instances Working Report CES-491,
	// School
	// of CS & EE, University of Essex, 02/2009.
	// http://dces.essex.ac.uk/staff/qzhang/MOEAcompetition/CEC09final/code/ZhangMOEADcode/moead0305.rar
	//algorithm.setInputParameter("dataDirectory","E:\\new_multiobjective\\jMetal\\Weights\\W2D_600.dat");
	algorithm.setInputParameter("dataDirectory","E:\\new_multiobjective\\jMetal\\Weights\\W3D_1000.dat");

	algorithm.setInputParameter("finalSize", 100); // used by MOEAD_DRA

	algorithm.setInputParameter("T", 20);
	algorithm.setInputParameter("delta", 0.9);
	algorithm.setInputParameter("nr", 2);

	// Crossover operator
	parameters = new HashMap();
	parameters.put("CR", 1.0);
	parameters.put("F", 0.5);
	crossover = CrossoverFactory.getCrossoverOperator(
			"DifferentialEvolutionCrossover", parameters);
	
	parameters = new HashMap(); 
	parameters.put("probability", 1.0);
	parameters.put("distributionIndex", 20.0); 
	crossover2 = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

	// Mutation operator
	parameters = new HashMap();
	parameters.put("probability", 1.0 / problem.getNumberOfVariables());
	parameters.put("distributionIndex", 20.0);
	mutation = MutationFactory.getMutationOperator("PolynomialMutation",
			parameters);

	algorithm.addOperator("crossover", crossover);
	algorithm.addOperator("crossover2", crossover2);
	algorithm.addOperator("mutation", mutation);

	// Execute the Algorithm
	long initTime = System.currentTimeMillis();
	SolutionSet population = algorithm.execute();
	//long estimatedTime = System.currentTimeMillis() - initTime;

	Execuion_time+=(System.currentTimeMillis() - initTime);
	
	// Result messages
	//logger_.info("Total execution time: " + estimatedTime + "ms");
	//logger_.info("Objectives values have been writen to file FUN");
	//population.printObjectivesToFile("FUN");
	//logger_.info("Variables values have been writen to file VAR");
	//population.printVariablesToFile("VAR");
	//population.printObjectivesToFile("FUN"+(i+1));

	//if (indicators != null) {
		// logger_.info("Quality indicators") ;
		// logger_.info("Hypervolume: " +
		// indicators.getHypervolume(population)) ;
		// logger_.info("EPSILON    : " + indicators.getEpsilon(population))
		// ;
		// logger_.info("GD         : " + indicators.getGD(population)) ;
		//logger_.info("IGD        : " + indicators.getCECIGD(population));
		// logger_.info("Spread     : " + indicators.getSpread(population))
		// ;
	//} // if
	population.printObjectivesToFile("FUN"+fun+"T"+(i+1));
	//population.printVariablesToFile("VAR"+fun+"T"+(i+1));
	GDarray[i]=indicators.getGD(population);
    IGDarray[i]=indicators.getCECIGD(population);
    spreadarray[i]=indicators.getSpread(population);
    Hypervolume[i]=indicators.getHypervolume(population);
	}
	printGD("FUN"+fun+"IGD",IGDarray);
	double sumGD=0;
	double sumIGD=0;
	double sumSP=0;
	for(int i=0;i<runtimes;i++){
	  sumGD+=GDarray[i];
	  sumIGD+=IGDarray[i];
	  sumSP+=spreadarray[i];
	}
	logger_.info("Total execution time: " + Execuion_time + "ms");
	System.out.println(sumGD/runtimes);
	System.out.println(sumSP/runtimes);
	System.out.println(sumIGD/runtimes);
	}
} // main
} // MOEAD_main
