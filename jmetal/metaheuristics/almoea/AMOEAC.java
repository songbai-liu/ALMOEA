package jmetal.metaheuristics.almoea;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.Spea2Fitness;
import jmetal.util.archive.SPEA2DensityArchive;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class AMOEAC extends Algorithm{
	int populationSize;
	int maxEvaluations;
	int evaluations;
	SolutionSet population;
	SolutionSet offspringPopulation;
	SolutionSet union;
	
	Distance distance = new Distance();
	
	int numVars;
	int numObjs;
	
	static double learningRate = 0.1;
	static int K = 10;
	static int numLayers = 3;
	MLPModel cMLP;
	MLPModel dMLP;

	double[] upBounds;
	double[] lowBounds;
	
	double[] zIdeal;
	double[] zNadir;
	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public AMOEAC(Problem problem) {
		super(problem);
	} // LMOEA
	
	public void initialization() throws JMException, ClassNotFoundException {
		// Read the parameters
		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
		
		// Initialize the variables
		population = new SolutionSet(populationSize);
		offspringPopulation = new SolutionSet(populationSize);
		union = new SolutionSet(2*populationSize);
		
		Solution newSolution;
		for (int j = 0; j < populationSize; j++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population.add(newSolution);
		} // for		
		evaluations = 0;
		
		numVars = problem_.getNumberOfVariables();
		numObjs = problem_.getNumberOfObjectives();
		
		upBounds = new double[numVars];
		lowBounds = new double[numVars];

		
		for(int var=0;var<numVars;var++) {
			upBounds[var] = problem_.getUpperLimit(var);
			lowBounds[var] = problem_.getLowerLimit(var);
		}
		
		zIdeal = new double[numObjs];
		zNadir = new double[numObjs];
		for(int m=0;m<numObjs;m++){
			zIdeal[m] = Double.MAX_VALUE;
			zNadir[m] = Double.MIN_VALUE;
		}
		
		//Construct the Auto-Encoders based on neural networks for each task
		cMLP = new MLPModel(numVars, K, numLayers, learningRate);
		dMLP = new MLPModel(numVars, K, numLayers, learningRate);
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		while(evaluations <= 1.0*maxEvaluations){
			AGSBalanceableCrossover1 crossover = new AGSBalanceableCrossover1(population, problem_,cMLP,dMLP);
			//AGSCompetitiveCrossover crossover = new AGSCompetitiveCrossover(population, problem_,cMLP);
            offspringPopulation = crossover.doCrossover();
            //cMLP = crossover.getCurrentMLPModel();
            cMLP = crossover.getConvergenceMLPModel();
            dMLP = crossover.getDiversityMLPModel();
			environmentalSelection();	
			evaluations += populationSize;
		}
		
		while(evaluations <= maxEvaluations){
			offspringPopulation = new CompetitiveDECrossover(population, problem_).doCrossover();
			environmentalSelection();	
			evaluations += populationSize;
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population);
		return final_ranking.getSubfront(0);
	}	
	
	public void environmentalSelection() {
		// Create the solutionSet union of solutionSet and offSpring
		union = ((SolutionSet) population).union(offspringPopulation);
		population.clear();
		offspringPopulation.clear();
		SolutionSet candidateSet = getStSolutionSet(union,(int)(1.1*populationSize));
		estimateIdealNadirPoint(candidateSet);
		normalization(candidateSet);
		List<SolutionSet> list = new <SolutionSet>ArrayList();
		for(int i=0;i<candidateSet.size();i++){
			SolutionSet sols = new SolutionSet();
			sols.add(candidateSet.get(i));
		    list.add(sols);  
		}
		list = new HierarchicalClustering(list).clusteringAnalysis(populationSize);
		if(list.size() != populationSize){
			System.out.println("The number of clusters after hierarchical clustering: ListSize = "+list.size());
			System.exit(0);
		}
		bestSolutionSelection(list);
	}//environmetalSelection
	
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
	
	public void estimateIdealNadirPoint(SolutionSet solSet){
		for(int p=0;p<solSet.size();p++){
			for(int i=0;i<numObjs;i++){
				if(solSet.get(p).getObjective(i) < zIdeal[i])
					zIdeal[i] = solSet.get(p).getObjective(i);
				if(solSet.get(p).getObjective(i) > zNadir[i])
					zNadir[i] = solSet.get(p).getObjective(i);
			}
		}
	}
	
	public void normalization(SolutionSet solSet){
		double value;
		for(int p=0;p<solSet.size();p++){
			double sum = 0.0;
			
			for(int i=0;i<numObjs;i++){
				value = (solSet.get(p).getObjective(i)-zIdeal[i])/(zNadir[i]-zIdeal[i]);
				//value = (solSet.get(p).getObjective(i)-zIdeal[i]);
				solSet.get(p).setNormalizedObjective(i, value);
				sum += value;
			}
			solSet.get(p).setSumValue(sum);
			
			for(int i=0;i<numObjs;i++){
				solSet.get(p).setIthTranslatedObjective(i, solSet.get(p).getNormalizedObjective(i)/sum);
			}
		}
	}
	
	public void bestSolutionSelection(List<SolutionSet> list){
    	
    	for(int k=0; k<problem_.getNumberOfObjectives();k++){
    		double minClustering2Axis = Math.acos(Math.abs(list.get(0).getCentroid().getNormalizedObjective(k)
    						/list.get(0).getCentroid().getDistanceToIdealPoint()));
      		int minClustering2AxisID = 0;
      		for(int i=1;i<list.size();i++){
      			SolutionSet sols = list.get(i);
      			if(sols.size() == 0){
      				System.out.println("ElitSolutionSelection_Diversity_SolsSize = "+sols.size());
      				System.exit(0);
      			}
      			
      			double angle1 = Math.acos(Math.abs(sols.getCentroid().getNormalizedObjective(k)/sols.getCentroid().getDistanceToIdealPoint()));
      			//System.out.println(angle1);
      			if(angle1 < minClustering2Axis){
      				minClustering2Axis = angle1;
      				minClustering2AxisID = i;
      			}//if
      		}//for
      		double minSolution2Axis = Math.acos(list.get(minClustering2AxisID).get(0).getNormalizedObjective(k)
      				/list.get(minClustering2AxisID).get(0).getDistanceToIdealPoint());;
      		int minSolution2AxisID = 0;
      		for(int j=1;j<list.get(minClustering2AxisID).size();j++){
      			Solution sol = list.get(minClustering2AxisID).get(j);
      			double ang = Math.acos(list.get(minClustering2AxisID).get(j).getNormalizedObjective(k)
      					/list.get(minClustering2AxisID).get(j).getDistanceToIdealPoint());
      			if(ang < minSolution2Axis){
      				minSolution2Axis = ang;
      				minSolution2AxisID = j;
      			}
      		}//for
      		double rnd = PseudoRandom.randDouble();
      		if(rnd < 1.1){
      			population.add(list.get(minClustering2AxisID).get(minSolution2AxisID));
              	list.remove(minClustering2AxisID);
      		}
    	}
    	
    	
	   Iterator<SolutionSet> it = list.iterator();
		while(it.hasNext()){
			SolutionSet sols = it.next();
			if(sols.size() == 0){
				System.out.println("In best solution selection, SolsSize2 = "+sols.size());
				System.exit(0);
			}
			double minFitness = 1.0e+30;
			int minFitnessID = -1;
			if(sols.size() == 0){
				System.out.println("size = 0!");
				System.exit(0);
			}
			for(int j=0;j<sols.size();j++){
				Solution sol2 = sols.get(j);
				double fitness = sol2.getSumValue();
				if(minFitness > fitness){
					minFitness = fitness;
					minFitnessID = j;
				}	
			}//for
			population.add(sols.get(minFitnessID));
			it.remove();
		}//while
		if(list.size() != 0){
			System.out.println("In best solution selection, ListSize2 = "+list.size());
			System.exit(0);
		}
   }
}
