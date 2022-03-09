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
import jmetal.problems.LMF.*;
import jmetal.problems.LSMOP.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
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
public class computeIGD {
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
		int n = 100;
		int m = 2;
		int runtimes=5;
		int tp1;//type of formulation model, 0 for addition, 1 for multiplication, and 2 for mixed/hybrid
		int tp2;//type of variable linkage, 0 for linear linkage, 1 for nonlinear linkage, and 2 for mixed/hybrid
		int tp3;//type of deep grouping on variables, 0 for even grouping, 1 for nonuniform grouping
		int tp4;//type of contribution on variables, 0 for balanced contribution, 1 for unbalanced contribution
		for(int fun=9;fun<=12;fun++){	
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
				if(fun == 1) {
	    		problem = new LMF1("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "concave_" + m +"D" +".pf");
	    	}
	    	if(fun == 2) {
	    		problem = new LMF2("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "convex_" + m +"D" +".pf");
	    	}
	    	if(fun == 3) {
	    		problem = new LMF3("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "linear_" + m +"D" +".pf");
	    	}
	    	if(fun == 4) {
	    		problem = new LMF4("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "concave_" + m +"D" +".pf");
	    	}
	    	if(fun == 5) {
	    		problem = new LMF5("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "convex_" + m +"D" +".pf");
	    	}
	    	if(fun == 6) {
	    		problem = new LMF6("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "linear_" + m +"D" +".pf");
	    	}
	    	if(fun == 7) {
	    		problem = new LMF7("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "concave_" + m +"D" +".pf");
	    	}
	    	if(fun == 8) {
	    		problem = new LMF8("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + "linear_" + m +"D" +".pf");
	    	}
	    	if(fun == 9) {
	    		m = 3;
	    		problem = new LMF9("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + problem.getHtype() + "_" + m +"D" +".pf");
	    	}
	    	if(fun == 10) {
	    		m = 3;
	    		problem = new LMF10("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + problem.getHtype() + "_" + m +"D" +".pf");
	    	}
	    	if(fun == 11) {
	    		m = 3;
	    		problem = new LMF11("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + problem.getHtype() + "_" + m +"D" +".pf");
	    	}
	    	if(fun == 12) {
	    		m = 3;
	    		problem = new LMF12("Real",n,m);
	    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LMF\\" + problem.getHtype() + "_" + m +"D" +".pf");
	    	}
		} // else
		utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
		String[] algorithms = new String[22];
		algorithms[0] = "DGEA";
		algorithms[1] = "LMOCSO";
		algorithms[2] = "LSMOF";
		algorithms[3] = "MOEADVA";
		algorithms[4] = "WOF";
	    tp1 = 2;
	    tp2 = 2;
	    tp3 = 1;
	    tp4 = 1;
	    int turn = fun + 10;
		for(int a=0;a<=4;a++){
			for(int i=0;i<runtimes;i++){
			   String path = "D:\\ResearchWork\\2021\\2021_P4_LMOP\\Review1\\Analysis_LMF\\" + tp1 + tp2 + tp3 + tp4;
			   path += "\\" + algorithms[a] + "\\" +algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+ "LMOP" + turn + "_" + n +"D_run"+(i+1)+".txt";
			   double[][] population = utils_.readFront(path); 
		       IGDarray[i] = indicators.getIGD1(population);
			}
			printGD(algorithms[a]+"_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + n + "D_IGD" + tp1 + tp2 + tp3 + tp4 + ".txt",IGDarray);
			double sumIGD = 0;
			for(int i=0;i<runtimes;i++){
			  sumIGD+=IGDarray[i];
			}
			logger_.info("Total execution time: " + Execuion_time + "ms");
			System.out.println("avrHV-fun"+fun+" = "+sumIGD/runtimes);
		}
	}//for
	} // main
} // MOEAD_main
