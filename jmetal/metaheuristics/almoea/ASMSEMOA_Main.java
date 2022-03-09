package jmetal.metaheuristics.almoea;

//  ASMSEMOA_main.java

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.LSMOP.*;
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

public class ASMSEMOA_Main {
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
            } // else
            m = problem.getNumberOfObjectives();
            for(int i=0;i<runtimes;i++){
                algorithm = new ASMSEMOA(problem);
                // algorithm = new MOEAD_DRA(problem);
                if(m == 2) {
                    algorithm.setInputParameter("populationSize", 100);
                    algorithm.setInputParameter("maxEvaluations", 100*100);
                }else if(m == 3) {
                    algorithm.setInputParameter("populationSize", 153);
                    algorithm.setInputParameter("maxEvaluations", 100*153);
                }

                algorithm.setInputParameter("offset", 100.0);

                // Mutation and Crossover for Real codification
                parameters = new HashMap();
                parameters.put("probability", 1.0);
                parameters.put("distributionIndex", 20.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
                        parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / problem.getNumberOfVariables());
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
                // Execute the Algorithm
                long initTime = System.currentTimeMillis();
                SolutionSet population = algorithm.execute();
                //long estimatedTime = System.currentTimeMillis() - initTime;

                Execuion_time+=(System.currentTimeMillis() - initTime);
               population.printObjectivesToFile("ASMSEMOA_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1)+".txt" );
                IGDarray[i]=indicators.getIGD1(population);
            }
            double sumIGD=0;
            for(int i=0;i<runtimes;i++){
                sumIGD+=IGDarray[i];
            }
            logger_.info("Total execution time: " + Execuion_time + "ms");
            System.out.println("IGD_"+ problem.getNumberOfObjectives()+"Obj_"+problem.getName()+
                    "_" + problem.getNumberOfVariables() + "D" + " = "+sumIGD/runtimes);
        }
    } // main
} // ASMSEMOA_Main
