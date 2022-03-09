package jmetal.metaheuristics.moead;

//  MOEAD.java
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


import jmetal.core.*;
import jmetal.metaheuristics.almoea.AGSMetaCrossover;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;

public class MOEAD_ANN1 extends Algorithm {

    private int populationSize_;

    private int H_;
    /**
     * Stores the population
     */
    private SolutionSet population;
    /**
     * Z vector (ideal point)
     */
    double[] z_;
    double[] N_;
    /**
     * Lambda vectors
     */
    // Vector<Vector<Double>> lambda_ ;
    double[][] lambda_;

    int evaluations_;

    /**
     * Constructor
     *
     * @param problem
     *            Problem to solve
     */
    public MOEAD_ANN1(Problem problem) {
        super(problem);
    } // MOEAD

    public SolutionSet execute() throws JMException, ClassNotFoundException {
        int maxEvaluations;

        evaluations_ = 0;
        maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations"))
                .intValue();
        populationSize_ = ((Integer) this.getInputParameter("populationSize"))
                .intValue();

        population = new SolutionSet(populationSize_);

        z_ = new double[problem_.getNumberOfObjectives()];
        N_ = new double[problem_.getNumberOfObjectives()];
        // lambda_ = new Vector(problem_.getNumberOfObjectives()) ;
        lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];

        H_= 16;

        // STEP 1. Initialization
        // STEP 1.1. Compute euclidean distances between weight vectors and find
        // T
        initUniformWeight();

        // STEP 1.2. Initialize population
        initPopulation();

        // STEP 1.3. Initialize z_
        initIdealAndNadirPoint();

        // STEP 2. Update
        do {
            SolutionSet offspringPopulation;
            if(PseudoRandom.randDouble(0.0,1.0) <0.5){
                offspringPopulation = new AGSMetaCrossover(population,problem_).doCrossover(0);
            }else{
                offspringPopulation = new AGSMetaCrossover(population,problem_).doCrossover(1);
            }

            SolutionSet union = ((SolutionSet) population).union(offspringPopulation);
            SolutionSet st = getStSolutionSet(union,populationSize_);
            for(int i = 0; i < populationSize_; i++){
                updateReference(st.get(i));
            };
            normalization(st);
            population.clear();
            population = new DecompositionBasedSelection(st, problem_.getNumberOfObjectives(), lambda_).environmentalSelection();
            evaluations_ += populationSize_;
        } while (evaluations_ < maxEvaluations);

        return population;
    }

    public SolutionSet getStSolutionSet(SolutionSet ss,int size) {
        Ranking ranking = new NondominatedRanking(ss);
        int remain = size;
        int index = 0;
        SolutionSet front = null;
        SolutionSet mgPopulation = new SolutionSet();
        front = ranking.getSubfront(index);
        while ((remain > 0) && (remain >= front.size())) {

            for (int k = 0; k < front.size(); k++) {
                mgPopulation.add(front.get(k));
            } // for
            // Decrement remain
            remain = remain - front.size();
            // Obtain the next front
            index++;
            if (remain > 0) {
                front = ranking.getSubfront(index);
            } // if
        }
        if (remain > 0) { // front contains individuals to insert
            for (int k = 0; k < front.size(); k++) {
                mgPopulation.add(front.get(k));
            }
        }
        return mgPopulation;
    }

    public void initUniformWeight() { // init lambda vectors
        int nw = 0;
        if (problem_.getNumberOfObjectives() == 2) {
            for (int n = 0; n < populationSize_; n++) {
                double a = 1.0 * n / (populationSize_ - 1);
                lambda_[n][0] = a;
                lambda_[n][1] = 1 - a;
                nw++;
            } // for
        } // if
        else if(problem_.getNumberOfObjectives() == 3) {
            int i, j;
            for (i = 0; i <= H_; i++) {
                for (j = 0; j <= H_; j++) {
                    if (i + j <= H_) {
                        lambda_[nw][0] = (double) (1.0 * i) / H_;
                        lambda_[nw][1] = (double) (1.0 * j) / H_;
                        lambda_[nw][2] = (double) (1.0 * (H_ - i - j) / H_);
                        nw++;
                    } // if
                } // for
            } // for
            if (nw != populationSize_) {
                System.out.println(nw + "---" + (populationSize_));
                System.out.println("ERROR: population size <> #weights");
                System.exit(0);
            }
        }else{
            System.out.println("This is not a Multiobjective Optimization Problem!");
            System.exit(0);
        }
    } // initUniformWeight


    /**
     *
     */
    public void initPopulation() throws JMException, ClassNotFoundException {
        for (int i = 0; i < populationSize_; i++) {
            Solution newSolution = new Solution(problem_);
            problem_.evaluate(newSolution);
            evaluations_++;
            population.add(newSolution);
        } // for
    } // initPopulation

    /**
     *
     */
    void initIdealAndNadirPoint() throws JMException, ClassNotFoundException {
        for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
            z_[i] = 1.0e+30;
            N_[i] = -1.0e+30;
        } // for
        for (int i = 0; i < populationSize_; i++) {
            updateReference(population.get(i));
        } // for
    } // initIdealPoint

    public void normalization(SolutionSet solSet){
        double value;
        for(int p=0;p<solSet.size();p++){
            double sum = 0.0;
            for(int i=0;i< problem_.getNumberOfObjectives();i++){
                //value = (solSet.get(p).getObjective(i)-z_[i])/ (N_[i]-z_[i]);
                value = (solSet.get(p).getObjective(i)-z_[i]);
                solSet.get(p).setNormalizedObjective(i, value);
                sum += value;
            }
            solSet.get(p).setSumValue(sum);
        }
    }

    /**
     *
     * @param individual
     */
    void updateReference(Solution individual) {
        for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
            if (individual.getObjective(n) < z_[n]) {
                z_[n] = individual.getObjective(n);
            }
            if (individual.getObjective(n) > N_[n]) {
                N_[n] = individual.getObjective(n);
            }
        }
    } // updateReference
} // MOEAD

