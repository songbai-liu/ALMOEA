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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.Vector;

public class MOEAD_ANN extends Algorithm {

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
    /**
     * Lambda vectors
     */
    // Vector<Vector<Double>> lambda_ ;
    double[][] lambda_;
    /**
     * T: neighbour size
     */
    int T_;
    /**
     * Neighborhood
     */
    int[][] neighborhood_;
    /**
     * delta: probability that parent solutions are selected from neighbourhood
     */
    double delta_;
    /**
     * nr: maximal number of solutions replaced by each child solution
     */
    int nr_;
    Solution[] indArray_;
    String functionType_;
    int evaluations_;
    /**
     * Operators
     */
    Operator crossover_;
    Operator mutation_;

    String dataDirectory_;

    /**
     * Constructor
     *
     * @param problem
     *            Problem to solve
     */
    public MOEAD_ANN(Problem problem) {
        super(problem);

        functionType_ = "_TCHE1";

    } // MOEAD

    public SolutionSet execute() throws JMException, ClassNotFoundException {
        int maxEvaluations;

        evaluations_ = 0;
        maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations"))
                .intValue();
        populationSize_ = ((Integer) this.getInputParameter("populationSize"))
                .intValue();
        //dataDirectory_ = this.getInputParameter("dataDirectory").toString();
        //System.out.println("POPSIZE: " + populationSize_);

        population = new SolutionSet(populationSize_);
        indArray_ = new Solution[problem_.getNumberOfObjectives()];

        //T_ = ((Integer) this.getInputParameter("T")).intValue();
        //nr_ = ((Integer) this.getInputParameter("nr")).intValue();
        //delta_ = ((Double) this.getInputParameter("delta")).doubleValue();


         T_ = 20;
         delta_ = 0.9;
         nr_ = 2;

        neighborhood_ = new int[populationSize_][T_];

        z_ = new double[problem_.getNumberOfObjectives()];
        // lambda_ = new Vector(problem_.getNumberOfObjectives()) ;
        lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];

        H_= 16;

        // STEP 1. Initialization
        // STEP 1.1. Compute euclidean distances between weight vectors and find
        // T
        initUniformWeight();
        // for (int i = 0; i < 300; i++)
        // System.out.println(lambda_[i][0] + " " + lambda_[i][1]) ;

        initNeighborhood();

        // STEP 1.2. Initialize population
        initPopulation();

        // STEP 1.3. Initialize z_
        initIdealPoint();

        // STEP 2. Update
        do {

            SolutionSet offspringPopulation = new AGSMetaCrossover(population,problem_).doCrossover(0);

            for (int i = 0; i < populationSize_; i++) {
                int type;
                double rnd = PseudoRandom.randDouble();
                if (rnd < delta_) // if (rnd < realb)
                {
                    type = 1; // neighborhood
                } else {
                    type = 2; // whole population
                }

                Solution child = offspringPopulation.get(i);
                evaluations_++;

                // STEP 2.4. Update z_
                updateReference(child);

                // STEP 2.5. Update of solutions
                updateProblem(child, i, type);
            } // for
        } while (evaluations_ < maxEvaluations);

        return population;
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
        }
        else{
            String dataFileName;
            dataFileName = "W" + problem_.getNumberOfObjectives() + "D_"
                    + populationSize_ + ".dat";

            try {
                // Open the file
                FileInputStream fis = new FileInputStream(dataDirectory_ + "/"
                        + dataFileName);
                InputStreamReader isr = new InputStreamReader(fis);
                BufferedReader br = new BufferedReader(isr);

                int i = 0;
                int j = 0;
                String aux = br.readLine();
                while (aux != null) {
                    StringTokenizer st = new StringTokenizer(aux);
                    j = 0;
                    while (st.hasMoreTokens()) {
                        double value = (new Double(st.nextToken())).doubleValue();
                        lambda_[i][j] = value;
                        j++;
                    }
                    aux = br.readLine();
                    i++;
                }
                br.close();
            } catch (Exception e) {
                System.out
                        .println("initUniformWeight: failed when reading for file: "
                                + dataDirectory_ + "/" + dataFileName);
                e.printStackTrace();
            }
        } // else


    } // initUniformWeight

    /**
     *
     */
    public void initNeighborhood() {
        double[] x = new double[populationSize_];
        int[] idx = new int[populationSize_];

        for (int i = 0; i < populationSize_; i++) {
            // calculate the distances based on weight vectors
            for (int j = 0; j < populationSize_; j++) {
                x[j] = Utils.distVector(lambda_[i], lambda_[j]);
                // x[j] = dist_vector(population[i].namda,population[j].namda);
                idx[j] = j;
                // System.out.println("x["+j+"]: "+x[j]+
                // ". idx["+j+"]: "+idx[j]) ;
            } // for

            // find 'niche' nearest neighboring subproblems
            Utils.minFastSort(x, idx, populationSize_, T_);
            // minfastsort(x,idx,population.size(),niche);

            System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
        } // for
    } // initNeighborhood

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
    void initIdealPoint() throws JMException, ClassNotFoundException {
        for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
            z_[i] = 1.0e+30;
            indArray_[i] = new Solution(problem_);
            problem_.evaluate(indArray_[i]);
            evaluations_++;
        } // for

        for (int i = 0; i < populationSize_; i++) {
            updateReference(population.get(i));
        } // for
    } // initIdealPoint

    /**
     *
     * @param individual
     */
    void updateReference(Solution individual) {
        for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
            if (individual.getObjective(n) < z_[n]) {
                z_[n] = individual.getObjective(n);

                indArray_[n] = individual;
            }
        }
    } // updateReference


    void updateProblem(Solution indiv, int id, int type) {
        // indiv: child solution
        // id: the id of current subproblem
        // type: update solutions in - neighborhood (1) or whole population
        // (otherwise)
        int size;
        int time;

        time = 0;

        if (type == 1) {
            size = neighborhood_[id].length;
        } else {
            size = population.size();
        }
        int[] perm = new int[size];

        Utils.randomPermutation(perm, size);

        for (int i = 0; i < size; i++) {
            int k;
            if (type == 1) {
                k = neighborhood_[id][perm[i]];
            } else {
                k = perm[i]; // calculate the values of objective function
                // regarding the current subproblem
            }
            double f1, f2;

            f1 = fitnessFunction(population.get(k), lambda_[k]);
            f2 = fitnessFunction(indiv, lambda_[k]);

            if (f2 < f1) {
                population.replace(k, new Solution(indiv));
                // population[k].indiv = indiv;
                time++;
            }
            // the maximal number of solutions updated is not allowed to exceed
            // 'limit'
            if (time >= nr_) {
                return;
            }
        }
    } // updateProblem

    double fitnessFunction(Solution individual, double[] lambda) {
        double fitness;
        fitness = 0.0;

        if (functionType_.equals("_TCHE1")) {
            double maxFun = -1.0e+30;

            for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
                double diff = Math.abs(individual.getObjective(n) - z_[n]);

                double feval;
                if (lambda[n] == 0) {
                    feval = 0.0001 * diff;
                } else {
                    feval = diff * lambda[n];
                }
                if (feval > maxFun) {
                    maxFun = feval;
                }
            } // for

            fitness = maxFun;
        } // if
        else {
            System.out.println("MOEAD.fitnessFunction: unknown type "
                    + functionType_);
            System.exit(-1);
        }
        return fitness;
    } // fitnessEvaluation
} // MOEAD

