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
import jmetal.problems.IMF.IMF1;
import jmetal.problems.IMF.IMF10;
import jmetal.problems.IMF.IMF2;
import jmetal.problems.IMF.IMF3;
import jmetal.problems.IMF.IMF4;
import jmetal.problems.IMF.IMF5;
import jmetal.problems.IMF.IMF6;
import jmetal.problems.IMF.IMF7;
import jmetal.problems.IMF.IMF8;
import jmetal.problems.IMF.IMF9;
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
public class computeIGD_main_3obj {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object
	
	public static jmetal.qualityIndicator.util.MetricsUtil utils_;

	/**
	 //* @param args
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
	int n = 1000;
	int runtimes=10;
	for(int fun=1;fun<=9;fun++){
	
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
			     problem = new LSMOP1("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP1_3D_10000.txt");
			    	}
			if(fun==2){
			     problem = new LSMOP2("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP2_3D_10000.txt");
			    	}
			if(fun==3){
			     problem = new LSMOP3("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP3_3D_10000.txt");
			    	}
			if(fun==4){
			     problem = new LSMOP4("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP4_3D_10000.txt");
			    	}
			if(fun==5){
			     problem = new LSMOP5("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP5_3D_10000.txt");
			    	}
			if(fun==6){
			     problem = new LSMOP6("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP6_3D_10000.txt");
			    	}
			if(fun==7){
			     problem = new LSMOP7("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP7_3D_10000.txt");
			    	}
			if(fun==8){
			     problem = new LSMOP8("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP8_3D_10000.txt");
			    	}
			if(fun==9){
			     problem = new LSMOP9("Real",n,3);
			     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP9_3D_10000.txt");
			}
			if(fun==10){
	    		problem = new IMF1("Real",n,2);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF1_2D_5000.txt");
		    }
	    	if(fun==11){
	    		problem = new IMF2("Real",n,2);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF2_2D_5000.txt");
		    }
	    	if(fun==12){
	    		problem = new IMF3("Real",n,2);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF3_2D_5000.txt");
		    }
	    	if(fun==13){
	    		problem = new IMF4("Real",n,3);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF4_3D_10000.txt");
		    }
	    	if(fun==14){
	    		problem = new IMF5("Real",n,2);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF5_2D_5000.txt");
		    }
	    	if(fun==15){
	    		problem = new IMF6("Real",n,2);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF6_2D_5000.txt");
		    }
	    	if(fun==16){
	    		problem = new IMF7("Real",n,2);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF7_2D_5000.txt");
		    }
	    	if(fun==17){
	    		problem = new IMF8("Real",n,3);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF8_3D_10000.txt");
		    }
	    	if(fun==18){
	    		problem = new IMF9("Real",n,2);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF9_2D_5000.txt");
		    }

		} // else
		utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
		String[] algorithms = new String[25];
		algorithms[1] = "LMOCSO";
		algorithms[2] = "LSMOF";
		algorithms[3] = "MLDE";
		algorithms[4] = "WOF";
		algorithms[5] = "IMMOEA";
		algorithms[6] = "NSGAII";
		algorithms[7] = "NSGAII_ANN";
		algorithms[8] = "SMSEMOA";
		algorithms[9] = "AMOEADV2";
		algorithms[10] = "AMOEADV3";
		algorithms[11] = "MOEAD_LDE";
		algorithms[13] = "SWTEA_FM2";
		algorithms[14] = "SWTEA_FM5";
		algorithms[15] = "SWTEA_FM10";
		algorithms[16] = "SWTEA_FM20";
		algorithms[17] = "SWTEA_FM50";
		algorithms[18] = "SWTEA_lg";
		
		algorithms[19] = "LVIDE_Inverse";
		algorithms[20] = "LVIDE_PS";
		algorithms[21] = "LVIDE_rndRank";
		algorithms[22] = "LVIDE_Alter";
		algorithms[23] = "LVIDE_Rand";
		algorithms[24] = "LMOEA";
		for(int a=9;a<=10;a++){
			for(int i=0;i<runtimes;i++){
			   //String path = "D:\\ResearchWork\\2021\\2021_SWTEA\\Experimental Studies\\SWTEA_Variants\\" + algorithms[a];
			//	String path = "D:\\2021_P1\\population\\FE_10000\\" + algorithms[a];
				//   path += "\\learing-DE\\m2\\"+algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName() + "_" + n +"D_run"+(i+1)+".txt";
			//	path +="\\"+ algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName() + "_" + n +"D_run"+(i+1)+".txt";

				String path = "D:\\ResearchWork\\2021\\2021_P8\\Review1\\AMOEAD\\" + algorithms[a];
				path +="\\"+ algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName() + "_" + n +"D_run"+(i+1)+".txt";
				//	path +="\\"+ algorithms[a]+"_"+problem.getName()+"_M"+problem.getNumberOfObjectives()+"_D"+problem.getNumberOfVariables()+"_"+(i+1)+".txt";
				double[][] population = utils_.readFront(path);
				IGDarray[i] = indicators.getIGD1(population);
			}
			printGD(algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+ problem.getName() + "_" + n + "D_IGD.txt",IGDarray);
			double sumIGD = 0;
			for(int i=0;i<runtimes;i++){
			  sumIGD+=IGDarray[i];
			}
			logger_.info("Total execution time: " + Execuion_time + "ms");
			System.out.println("avrIGD-fun"+fun+" = "+sumIGD/runtimes);
		}
	 }//for
	} // main
} // MOEAD_main
