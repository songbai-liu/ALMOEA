package jmetal.metaheuristics.nsgaII;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.ranking.NondominatedRanking;

/**
 * Implementation of NSGA-II. This implementation of NSGA-II makes use of a
 * QualityIndicator object to obtained the convergence speed of the algorithm.
 * This version is used in the paper: A.J. Nebro, J.J. Durillo, C.A. Coello
 * Coello, F. Luna, E. Alba
 * "A Study of Convergence Speed in Multi-Objective Metaheuristics." To be
 * presented in: PPSN'08. Dortmund. September 2008.
 */

public class NSGADE extends Algorithm {
	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public NSGADE(Problem problem) {
		super(problem);
	} // NSGAII

	
	/**
	 * Runs the NSGA-II algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int populationSize;
		int maxEvaluations;
		int evaluations;
		double[] IGD;
		QualityIndicator indicators; // QualityIndicator object
		int requiredEvaluations; // Use in the example of use of the
		// indicators object (see below)
		int neighborsize=20;
		int[][] neighborhood;

		SolutionSet population;
		SolutionSet offspringPopulation;
		SolutionSet union;

		Operator mutationOperator;
		Operator crossoverOperator;
		Operator selectionOperator;

		Distance distance = new Distance();
		neighborhood=new int[problem_.getNumberOfObjectives()][neighborsize];

		// Read the parameters
		populationSize = ((Integer) getInputParameter("swarmSize"))
				.intValue();
		maxEvaluations = ((Integer) getInputParameter("maxIterations"))
				.intValue();
		indicators = (QualityIndicator) getInputParameter("indicators");

		// Initialize the variables
		population = new SolutionSet(populationSize);
		evaluations = 0;

		requiredEvaluations = 0;

		// Read the operators
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");
		Distance distance_ = new Distance();
		IGD =  new double[maxEvaluations/100];
		// Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			//evaluations++;
			population.add(newSolution);
		} // for
		
		//findneighbor();
        int t = 0;
		// Generations
		while (evaluations < maxEvaluations) {
			if(evaluations % 100 == 0){
				IGD[t] = indicators.getIGD1(population);
				t++;
			}
			// Create the offSpring solutionSet
			offspringPopulation = new SolutionSet(populationSize);
			Solution[] parents = new Solution[4];
			for (int i = 0; i < populationSize; i++) {
				if (evaluations < maxEvaluations) {
					// obtain parents
					parents[2] = (Solution) selectionOperator
							.execute(population);
					parents[0] = (Solution) selectionOperator
							.execute(population);
					parents[1] = (Solution) selectionOperator
							.execute(population);
					double mindistance=distance_.distanceBetweenSolutions(parents[0], parents[2]);
					double mindistance2=distance_.distanceBetweenSolutions(parents[1], parents[2]);
					if(mindistance>mindistance2)
						parents[0]=parents[1];
					parents[3] = (Solution) selectionOperator
							.execute(population);
					parents[1] = (Solution) selectionOperator
							.execute(population);
					mindistance=distance_.distanceBetweenSolutions(parents[3], parents[2]);
					mindistance2=distance_.distanceBetweenSolutions(parents[1], parents[2]);
					if(mindistance<mindistance2)
						parents[1]=parents[3];
					Solution child = (Solution) crossoverOperator.execute(new Object[] {
							parents[2], parents});
					mutationOperator.execute(child);
					//mutationOperator.execute(offSpring[1]);
					problem_.evaluate(child);
					problem_.evaluateConstraints(child);
					offspringPopulation.add(child);
				} // if
			} // for

			// Create the solutionSet union of solutionSet and offSpring
			union = ((SolutionSet) population).union(offspringPopulation);

			// Ranking the union
			NondominatedRanking ranking = new NondominatedRanking(union);

			int remain = populationSize;
			int index = 0;
			SolutionSet front = null;
			population.clear();

			// Obtain the next front
			front = ranking.getSubfront(index);

			while ((remain > 0) && (remain >= front.size())) {
				// Assign crowding distance to individuals
				distance.crowdingDistanceAssignment(front,
						problem_.getNumberOfObjectives());
				// Add the individuals of this front
				for (int k = 0; k < front.size(); k++) {
					population.add(front.get(k));
				} // for

				// Decrement remain
				remain = remain - front.size();

				// Obtain the next front
				index++;
				if (remain > 0) {
					front = ranking.getSubfront(index);
				} // if
			} // while

			// Remain is less than front(index).size, insert only the best one
			if (remain > 0) { // front contains individuals to insert
				distance.crowdingDistanceAssignment(front,
						problem_.getNumberOfObjectives());
				front.sort(new CrowdingComparator());
				for (int k = 0; k < remain; k++) {
					population.add(front.get(k));
				} // for

				remain = 0;
			} // if
			evaluations += 1;

			// This piece of code shows how to use the indicator object into the
			// code
			// of NSGA-II. In particular, it finds the number of evaluations
			// required
			// by the algorithm to obtain a Pareto front with a hypervolume
			// higher
			// than the hypervolume of the true Pareto front.
			/*
			 * if ((indicators != null) && (requiredEvaluations == 0)) { double HV =
			 * indicators.getHypervolume(population); if (HV >= (0.98 *
			 * indicators.getTrueParetoFrontHypervolume())) { requiredEvaluations =
			 * evaluations; } // if } // if
			 */		
		} // while
		printGD("NSGADE_"+problem_.getNumberOfObjectives()+"Obj_"+problem_.getName()+ "_" + problem_.getNumberOfVariables() + "D_IGD"+".txt", IGD);
		// Return as output parameter the required evaluations
		//setOutputParameter("evaluations", requiredEvaluations);

		// Return the first non-dominated front
		NondominatedRanking ranking = new NondominatedRanking(population);
		//ranking.getSubfront(0).printFeasibleFUN("FUN_NSGAII");

		return ranking.getSubfront(0);
	} // execute
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");//写到缓冲区
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
} // NSGA-II

