//  ANSGAII.java
//
//  Authors:
//       Songbai Liu <songbai209@qq.com>
//       Jun Li <lijun@szu.edu.cn>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
package jmetal.metaheuristics.almoea;

import jmetal.core.*;
import jmetal.metaheuristics.almoea.AGSMetaCrossover;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.ranking.NondominatedRanking;

/**
 * Implementation of ANSGA-II.
 */

public class ANSGAII extends Algorithm {
	
	static double learningRate = 0.1;
    static int K = 10;
    static int numLayers = 3;
    private MLPModel MLP;
    /**
     * Constructor
     *
     * @param problem
     *            Problem to solve
     */
    public ANSGAII(Problem problem) {
        super(problem);
    } // NSGAII

    /**
     * Runs the ANSGA-II algorithm.
     *
     * @return a <code>SolutionSet</code> that is a set of non dominated
     *         solutions as a result of the algorithm execution
     * @throws JMException
     */
    public SolutionSet execute() throws JMException, ClassNotFoundException {
    	
    	MLP = new MLPModel(problem_.getNumberOfVariables(), K, numLayers, learningRate);
    	
        int populationSize;
        int maxEvaluations;
        int evaluations;

        QualityIndicator indicators; // QualityIndicator object
        int requiredEvaluations; // Use in the example of use of the

        SolutionSet population;
        SolutionSet offspringPopulation;
        SolutionSet union;

        Operator mutationOperator;
        Operator crossoverOperator;
        Operator selectionOperator;

        Distance distance = new Distance();

        // Read the parameters
        populationSize = ((Integer) getInputParameter("populationSize"))
                .intValue();
        maxEvaluations = ((Integer) getInputParameter("maxEvaluations"))
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

        // Create the initial solutionSet
        Solution newSolution;
        for (int i = 0; i < populationSize; i++) {
            newSolution = new Solution(problem_);
            problem_.evaluate(newSolution);
            problem_.evaluateConstraints(newSolution);
            //evaluations++;
            population.add(newSolution);
        } // for

        while (evaluations < maxEvaluations) {

            // Create the offSpring solutionSet
          //  offspringPopulation = new anncrossover1(population,problem_).doCrossover();
        	AGSMetaCrossover crossover = new AGSMetaCrossover(population,problem_,MLP);
            if(evaluations < 0.95*maxEvaluations){
                offspringPopulation = crossover.doCrossover(0);
            }else{
                offspringPopulation = crossover.doCrossover(1);
            }
            MLP = crossover.getCurrentMLPModel();
            
            evaluations += populationSize;
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
        } // while


        NondominatedRanking ranking = new NondominatedRanking(population);

        return ranking.getSubfront(0);
    } // execute
} // ANSGA-II
