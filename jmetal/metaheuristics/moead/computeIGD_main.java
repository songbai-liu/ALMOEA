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
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
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
public class computeIGD_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object
	
	public static jmetal.qualityIndicator.util.MetricsUtil utils_;

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
		int[] N = new int[5];
		N[0] = 10;
		N[1] = 50;
		N[2] = 100;
		N[3] = 500;
		N[4] = 1000;
		for(int sb=0;sb<5;sb++){
			int n = N[sb];
			int runtimes=12;
			int t = 1;
			if(sb==0){
				t=10;
			}
			for(int fun=t;fun<=57;fun++){	
			
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
			    	if(fun==19){
			    		problem = new IMF10("Real",n,2);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF10_2D_5000.txt");
				    }
			    	if(fun==20){
			    		problem = new UF1("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF1_2D_1000.txt");
				    }
			    	if(fun==21){
			    		problem = new UF2("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF2_2D_1000.txt");
				    }
			    	if(fun==22){
			    		problem = new UF3("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF3_2D_1000.txt");
				    }
			    	if(fun==23){
			    		problem = new UF4("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF4_2D_1000.txt");
				    }
			    	if(fun==24){
			    		problem = new UF5("Real",n,10,0.1);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF5_21.txt");
				    }
			    	if(fun==25){
			    		problem = new UF6("Real",n,2,0.1);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF6_668.txt");
				    }
			    	if(fun==26){
			    		problem = new UF7("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF7_2D_1000.txt");
				    }
			    	if(fun==27){
			    		problem = new UF8("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF8_3D_10000.txt");
				    }
			    	if(fun==28){
			    		problem = new UF9("Real",n,0.1);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF9_3D_10000.txt");
				    }
			    	if(fun==29){
			    		problem = new UF10("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF10_3D_10000.txt");
				    }
			    	if(fun==30){
					     problem = new WFG1("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG1_3D_10000.txt");
					    	}
					if(fun==31){
						 problem = new WFG2("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG2_3D_10000.txt");
					    	}
					if(fun==32){
						 problem = new WFG3("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG3_3D_10000.txt");
					    	}
					if(fun==33){
						 problem = new WFG4("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG9_3D_10000.txt");
					    	}
					if(fun==34){
						 problem = new WFG5("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG9_3D_10000.txt");
					    	}
					if(fun==35){
						 problem = new WFG6("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG9_3D_10000.txt");
					    	}
					if(fun==36){
						 problem = new WFG7("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG9_3D_10000.txt");
					    	}
					if(fun==37){
						 problem = new WFG8("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG9_3D_10000.txt");
					    	}
					if(fun==38){
						 problem = new WFG9("Real",4,n-4,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\WFG\\WFG9_2D_5000.txt");
					}
					if(fun==39){
					     problem = new DTLZ1("Real",n,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\DTLZ1_3D_10000.txt");
					}
					if(fun==40){
						 problem = new DTLZ2("Real",n,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\DTLZ2_3D_10000.txt");
					    	}
					if(fun==41){
						 problem = new DTLZ3("Real",n,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\DTLZ2_3D_10000.txt");
					}
					if(fun==42){
						 problem = new DTLZ4("Real",n,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\DTLZ2_3D_10000.txt");
					    	}
					if(fun==43){
						 problem = new DTLZ5("Real",n,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\DTLZ5_3D_10000.txt");
					}
					if(fun==44){
						 problem = new DTLZ6("Real",n,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\DTLZ5_3D_10000.txt");
					}
					if(fun==45){
						 problem = new DTLZ7("Real",n,3);
					     indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\DTLZ7_3D_10000.txt");
					}
					if(fun==46){
			    		problem = new MOP1("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\MOP\\MOP1_2D_5000.txt");
				    }
			    	if(fun==47){
			    		problem = new MOP2("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\MOP\\MOP2_2D_5000.txt");
				    }
			    	if(fun==48){
			    		problem = new MOP3("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\MOP\\MOP3_2D_5000.txt");
				    }
			    	if(fun==49){
			    		problem = new MOP4("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\MOP\\MOP4_2D_5000.txt");
				    }
			    	if(fun==50){
			    		problem = new MOP5("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\MOP\\MOP5_2D_5000.txt");
				    }
			    	if(fun==51){
			    		problem = new MOP6("Real",n,3);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\MOP\\MOP6_3D_5000.txt");
				    }
			    	if(fun==52){
			    		problem = new MOP7("Real",n,3);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\MOP\\MOP7_3D_5000.txt");
				    }
			    	if(fun==53){
			    		problem = new ZDT1("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\ZDT\\ZDT1_2D.txt");
				    }
			    	if(fun==54){
			    		problem = new ZDT2("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\ZDT\\ZDT2_2D.txt");
				    }
			    	if(fun==55){
			    		problem = new ZDT3("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\ZDT\\ZDT3_2D.txt");
				    }
			    	if(fun==56){
			    		problem = new ZDT4("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\ZDT\\ZDT4_2D.txt");
				    }
			    	if(fun==57){
			    		problem = new ZDT6("Real",n);
			    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\ZDT\\ZDT6_2D.txt");
				    }
			    
				} // else
				utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
				String[] algorithms = new String[22];
				algorithms[0] = "MOEADDE";
				algorithms[1] = "LMOCSO";
				algorithms[2] = "LSMOF";
				algorithms[3] = "MOEADVA";
				algorithms[4] = "WOF";
				algorithms[5] = "IMMOEA";
				algorithms[6] = "NSGAII";
				algorithms[7] = "MOEAD";
				algorithms[8] = "IBEA";
				algorithms[9] = "SMPSO";
				algorithms[10] = "DGEA";
				algorithms[11] = "MOEADM2M";
				
				for(int a=0;a<=11;a++){
					for(int i=0;i<runtimes;i++){
					   String path = "D:\\CityU\\2021\\BenchmarkProblems\\ExperimentResults\\10MFE\\" + algorithms[a];
					   if(fun>=10 && fun<=19){//IMF1 to IMF10
						   path += "\\"+algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+ "IMMOEA_F" + (fun-9) + "_" + n +"D_run"+(i+1)+".txt";
					   }else if(fun>=46 && fun<=52){//MOP1 to MOP7
						   path += "\\"+algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+ "MOEADM2M_F" + (fun-45) + "_" + n +"D_run"+(i+1)+".txt";
					   }else{
						   path += "\\"+algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+ problem.getName() + "_" + n +"D_run"+(i+1)+".txt";
					   }
					   
					  
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
		}
	
	} // main
} // MOEAD_main
