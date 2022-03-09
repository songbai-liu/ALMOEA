package jmetal.metaheuristics.almoea;

//  AMOEAD_main.java

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.problems.LSMOP.*;
import jmetal.problems.ProblemFactory;
import jmetal.problems.DTLZ.*;
import jmetal.problems.LMF.LMF1;
import jmetal.problems.LMF.LMF10;
import jmetal.problems.LMF.LMF11;
import jmetal.problems.LMF.LMF12;
import jmetal.problems.LMF.LMF2;
import jmetal.problems.LMF.LMF3;
import jmetal.problems.LMF.LMF4;
import jmetal.problems.LMF.LMF5;
import jmetal.problems.LMF.LMF6;
import jmetal.problems.LMF.LMF7;
import jmetal.problems.LMF.LMF8;
import jmetal.problems.LMF.LMF9;
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

public class AMOEAD_Main {
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
        int n = 5000;
        int m = 2;
        int num_samples = 5000;
        if(m == 3) {
            num_samples = 10000;
        }

        Problem problem=null; // The problem to solve
        Algorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator

        QualityIndicator indicators=null; // Object to get quality indicators

        HashMap parameters; // Operator parameters

        indicators = null;
        for(int fun=1;fun<=9;fun++){
            int runtimes=1;
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
                if(fun==11){
                    problem = new DTLZ1("Real",n,m);
                    indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\" + "DTLZ1" + "_" + m +"D"+".txt");
                }
                if(fun==12){
                    problem = new DTLZ2("Real",n,m);
                    indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\" + "DTLZ2" + "_" + m +"D"+".txt");
                }
                if(fun==13){
                    problem = new DTLZ3("Real",n,m);
                    indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\" + "DTLZ2" + "_" + m +"D"+".txt");
                }
                if(fun==14){
                    problem = new DTLZ4("Real",n,m);
                    indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\DTLZ\\" + "DTLZ2" + "_" + m +"D"+".txt");
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
                //algorithm = new AMOEAD_Hybrid(problem);
                algorithm = new AMOEAD_Innov(problem);
                if(m == 2) {
                    algorithm.setInputParameter("populationSize", 100);
                    algorithm.setInputParameter("maxEvaluations", 100*1000);
                }else if(m == 3) {
                    algorithm.setInputParameter("populationSize", 153);
                    algorithm.setInputParameter("maxEvaluations", 100*153);
                }

                // Execute the Algorithm
                long initTime = System.currentTimeMillis();
                SolutionSet population = algorithm.execute();
                //long estimatedTime = System.currentTimeMillis() - initTime;

                Execuion_time+=(System.currentTimeMillis() - initTime);
                //population.printObjectivesToFile("AMOEADIn_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1)+".txt" );
                IGDarray[i] = indicators.getIGD1(population);
            }
            double sumIGD=0;
            for(int i=0;i<runtimes;i++){
                sumIGD+=IGDarray[i];
            }

            System.out.println("IGD_"+ problem.getNumberOfObjectives()+"Obj_"+problem.getName()+
                    "_" + problem.getNumberOfVariables() + "D" + " = "+sumIGD/runtimes);
        }
    } // main
} // AMOEAD_Main
