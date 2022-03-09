//  SMPSO_main.java
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

package jmetal.metaheuristics.smpso;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.mutation.Mutation;
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
import jmetal.problems.LSMOP.LSMOP1;
import jmetal.problems.LSMOP.LSMOP2;
import jmetal.problems.LSMOP.LSMOP3;
import jmetal.problems.LSMOP.LSMOP4;
import jmetal.problems.LSMOP.LSMOP5;
import jmetal.problems.LSMOP.LSMOP6;
import jmetal.problems.LSMOP.LSMOP7;
import jmetal.problems.LSMOP.LSMOP8;
import jmetal.problems.LSMOP.LSMOP9;
import jmetal.problems.M2M.MOP1;
import jmetal.problems.M2M.MOP2;
import jmetal.problems.M2M.MOP3;
import jmetal.problems.M2M.MOP4;
import jmetal.problems.M2M.MOP5;
import jmetal.problems.M2M.MOP6;
import jmetal.problems.M2M.MOP7;
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
import jmetal.problems.cec2009Competition.UF1;
import jmetal.problems.cec2009Competition.UF10;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;
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
 * This class executes the SMPSO algorithm described in: A.J. Nebro, J.J.
 * Durillo, J. Garcia-Nieto, C.A. Coello Coello, F. Luna and E. Alba
 * "SMPSO: A New PSO-based Metaheuristic for Multi-objective Optimization". IEEE
 * Symposium on Computational Intelligence in Multicriteria Decision-Making
 * (MCDM 2009), pp: 66-73. March 2009
 */
public class SMPSO_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args
	 *            Command line arguments. The first (optional) argument
	 *            specifies the problem to solve.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options - jmetal.metaheuristics.smpso.SMPSO_main
	 *             - jmetal.metaheuristics.smpso.SMPSO_main problemName -
	 *             jmetal.metaheuristics.smpso.SMPSO_main problemName
	 *             ParetoFrontFile
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
	
	public static void main(String[] args) throws JMException, IOException,
			ClassNotFoundException {

		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("SMPSO_main.log");
		logger_.addHandler(fileHandler_);

		for(int fun=29;fun<=37;fun++){
			int runtimes=1;
			//double[] GDarray=new double[runtimes];
			  double[] IGDarray=new double[runtimes];
			  //double[] spreadarray=new double[runtimes];
			  //double[] Hypervolume=new double[runtimes];
			  long Execution_time=0;

		for(int i=0;i<runtimes;i++){
		Problem problem=null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Mutation mutation; // "Turbulence" operator

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
		  	      problem = new ZDT1("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT1_501.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==2){
		  	      problem = new ZDT2("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT2_501.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==3){
		  	      problem = new ZDT3("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT3_269.txt" ) ;
		  	    	}
		  	if(fun==4){
			      problem = new ZDT4("Real",10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT4_501.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==5){
			      problem = new ZDT6("Real",10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT6_774.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==6){
				//problem = new DTLZ1("Real",10,2);
			      problem = new DTLZ1("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ1_M10.dat" ) ;
			    	}
		  	if(fun==7){
			      //problem = new DTLZ2("Real",10,2);
		  		problem = new DTLZ2("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==8){
			      problem = new DTLZ3("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==9){
			      problem = new DTLZ4("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}
		  	if(fun==10){
			      problem = new DTLZ5("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ5_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==11){
			      problem = new DTLZ6("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ5_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==12){
			      problem = new DTLZ7("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ7_M10.dat" ) ;
			    	}
			if(fun==13){
		  	      problem = new WFG1("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG1_605.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==14){
		  	      problem = new WFG2("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG2_111.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==15){
		  	      problem = new WFG3("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG3_right.txt" ) ;
		  	    	}
			if(fun==16){
			      problem = new WFG4("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG4_1181.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==17){
			      problem = new WFG5("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG5_694.txt" ) ;
			    	}
		  	if(fun==18){
		  	      problem = new WFG6("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\wfg6_166.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==19){
			      problem = new WFG7("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG7_2435.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==20){
			      problem = new WFG8("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG8_201.txt" );
			    	}
			if(fun==21){
			      problem = new WFG9("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG9_2591.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==29){
			      problem = new LSMOP1("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP1_2D_5000.txt");
			    	}
			if(fun==30){
			      problem = new LSMOP2("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP2_2D_5000.txt");
			    	}
			if(fun==31){
			      problem = new LSMOP3("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP3_2D_5000.txt");
			    	}
			if(fun==32){
			      problem = new LSMOP4("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP4_2D_5000.txt");
			    	}
			if(fun==33){
			      problem = new LSMOP5("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP5_2D_5000.txt");
			    	}
			if(fun==34){
			      problem = new LSMOP6("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP6_2D_5000.txt");
			    	}
			if(fun==35){
			      problem = new LSMOP7("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP7_2D_5000.txt");
			    	}
			if(fun==36){
			      problem = new LSMOP8("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP8_2D_5000.txt");
			    	}
			if(fun==37){
			      problem = new LSMOP9("Real",200,2);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP9_2D_5000.txt");
			    	}
			if(fun==38){
			      problem = new UF1("Real",100);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF1_500.txt");
			    	}
			if(fun==39){
			      problem = new UF2("Real",100);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF2_500.txt");
			    	}
			if(fun==40){
			      problem = new UF3("Real",100);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF3_500.txt");
			    	}
			if(fun==41){
			      problem = new UF4("Real",100);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF4_500.txt");
			    	}
			if(fun==42){
			      problem = new UF5("Real",100,2,0.1);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF5_21.txt");
			    	}
			if(fun==43){
			      problem = new UF6("Real",100,2,0.1);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF6_668.txt");
			    	}
			if(fun==44){
			      problem = new UF7("Real",100);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF7_500.txt");
			    	}
			if(fun==45){
			      problem = new UF8("Real");;
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF8_500.txt");
			    	}
			if(fun==46){
			      problem = new UF9("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==47){
			      problem = new UF10("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==48){
			      problem = new MOP1("Real");
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==49){
			      problem = new MOP2("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==50){
			      problem = new MOP3("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==51){
			      problem = new MOP4("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==52){
			      problem = new MOP5("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==53){
			      problem = new MOP6("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==54){
			      problem = new MOP7("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			
		} // else

		algorithm = new SMPSO(problem);

		// Algorithm parameters
		if(fun>=29 && fun <= 37) {
			algorithm.setInputParameter("swarmSize", 300);
			algorithm.setInputParameter("archiveSize", 300);
			algorithm.setInputParameter("maxIterations", 10000);
			algorithm.setInputParameter("div1", 299);
			algorithm.setInputParameter("div2", 0);
			algorithm.setInputParameter("T", 20);
		}else if(fun >= 38 && fun <= 44) {
			algorithm.setInputParameter("swarmSize", 600);
			algorithm.setInputParameter("archiveSize", 600);
			algorithm.setInputParameter("maxIterations", 500);
		}else if(fun >= 42 && fun <= 50){
			algorithm.setInputParameter("swarmSize", 300);
			algorithm.setInputParameter("archiveSize", 300);
			algorithm.setInputParameter("maxIterations", 4000);
			algorithm.setInputParameter("div1", 23);
			algorithm.setInputParameter("div2", 0);
			algorithm.setInputParameter("T", 20);
		}
		

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation",
				parameters);

		algorithm.addOperator("mutation", mutation);

		// Execute the Algorithm
		long initTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		//long estimatedTime = System.currentTimeMillis() - initTime;
		Execution_time+=(System.currentTimeMillis() - initTime);

		// Result messages
		/*logger_.info("Total execution time: " + estimatedTime + "ms");
		logger_.info("Objectives values have been writen to file FUN");
		population.printObjectivesToFile("FUN");
		logger_.info("Variables values have been writen to file VAR");
		population.printVariablesToFile("VAR");

		if (indicators != null) {
			logger_.info("Quality indicators");
			logger_.info("Hypervolume: "
					+ indicators.getHypervolume(population));
			logger_.info("GD         : " + indicators.getGD(population));
			logger_.info("IGD        : " + indicators.getIGD(population));
			logger_.info("Spread     : " + indicators.getSpread(population));
			logger_.info("Epsilon    : " + indicators.getEpsilon(population));
		} // if*/
		//population.printObjectivesToFile("Run"+ i + "-FUN-" + problem.getName()
				//+ "-SMPSO-FUN");
		//GDarray[i]=indicators.getCECIGD(population);
	   IGDarray[i]=indicators.getIGD1(population);
	   // spreadarray[i]=indicators.getEpsilon(population);
	   // Hypervolume[i]=indicators.getHypervolume(population);
		//population.printObjectivesToFile("SMPSO_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_run"+(i+1)+".txt" );
		}
	//printGD("FUN"+fun+"IGD",GDarray);
		//printGD("FUN"+fun+"IGD",IGDarray);
		//printGD("FUN"+fun+"IGD",IGDarray);
		//printDiversity("FUN"+fun+"Diversity",spreadarray);
		//printHypervolume("FUN"+fun+"Hypervolume",Hypervolume);
		//double sumGD=0;
		double sumIGD=0;
		//double sumSP=0;
		//double sumHypervolume=0.0;
		for(int i=0;i<runtimes;i++){
		 // sumGD+=GDarray[i];
		  sumIGD+=IGDarray[i];
		  //sumSP+=spreadarray[i];
		 // sumHypervolume+=Hypervolume[i];
		}
		logger_.info("Total execution time: " + Execution_time + "ms");
		//System.out.println("avrGD-fun"+fun+" = "+sumGD/runtimes);
		//System.out.println("avrSP-fun"+fun+" = "+sumSP/runtimes);
		System.out.println("avrIGD-fun"+fun+"= "+sumIGD/runtimes);
		//System.out.println("avrHV-fun"+fun+" = "+sumHypervolume/runtimes);
		}
	} // main
} // SMPSO_main
