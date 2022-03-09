package jmetal.metaheuristics.almoea;

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

public class ALMOEA extends Algorithm{
	int populationSize;
	int maxEvaluations;
	int evaluations;
	SolutionSet population;
	SolutionSet offspringPopulation;
	SolutionSet union;
	
	Operator mutationOperator;
	Operator crossoverOperator;
	Operator selectionOperator;
	
	Distance distance = new Distance();
	
	int numVars;
	int numObjs;
	
	static double learningRate = 0.1;
	static int K = 10;
	static int numLayers = 3;
	MLPModel MLP;

	double[] upBounds;
	double[] lowBounds;
	
	double[] zIdeal;
	double[] zNadir;
	double[][] lamda1;
	double[][] lamda2;
	int div1_,div2_;
	int div3_,div4_;
	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public ALMOEA(Problem problem) {
		super(problem);
	} // LMOEA
	
	public void initialization() throws JMException, ClassNotFoundException {
		// Read the parameters
		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
		div1_ = ((Integer) getInputParameter("div1")).intValue();
		div2_ = ((Integer) getInputParameter("div2")).intValue();
		div3_ = ((Integer) getInputParameter("div3")).intValue();
		div4_ = ((Integer) getInputParameter("div4")).intValue();
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");
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
		// Read the operators
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");
		
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
		
		/* generate two-layer weight vectors */
		VectorGenerator vg;
		vg = new TwoLevelWeightVectorGenerator(div1_, div2_,numObjs);
		lamda1 = vg.getVectors();
		vg = new TwoLevelWeightVectorGenerator(div3_, div4_,numObjs);
		lamda2 = vg.getVectors();
		//Construct the Auto-Encoders based on neural networks for each task
		MLP = new MLPModel(numVars, K, numLayers, learningRate);
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		while(evaluations <= 0.2*maxEvaluations){
			AGSCompetitiveCrossover crossover = new AGSCompetitiveCrossover(population, problem_,MLP);
            offspringPopulation = crossover.doCrossover();
            MLP = crossover.getCurrentMLPModel();
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
	
	public double[] learningViaDAE(XReal xsol) throws JMException, ClassNotFoundException {
		double[] input = new double[numVars];
		for(int i=0;i<numVars;i++) { 
			input[i] = (xsol.getValue(i) - lowBounds[i])/(upBounds[i] - lowBounds[i]); 
		}
		double[] output = MLP.baseMLP.computeOut(input);
		double value;
		for(int var=0;var<numVars;var++) {
			value = output[var]*(upBounds[var] - lowBounds[var]) + lowBounds[var];
			if(value < lowBounds[var]) {
				value = lowBounds[var];
			}
			if(value > upBounds[var]) {
				value = upBounds[var];
			}
			output[var] = value;
		}
		return output;
	}
	
	public void environmentalSelection() {
		// Create the solutionSet union of solutionSet and offSpring
		union = ((SolutionSet) population).union(offspringPopulation);
		population.clear();
		offspringPopulation.clear();
		SolutionSet st = getStSolutionSet(union,populationSize);
		estimateIdealNadirPoint(st);
		normalization(union);
		population = new DecompositionBasedSelection(union, numObjs, lamda2).environmentalSelection();
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
	
	public SolutionSet getLabel(SolutionSet[] input){
		 int size = input[1].size();
	        SolutionSet labelSet = new SolutionSet(size);
	        for(int i=0;i<size;i++) {
	        	int rdInt = PseudoRandom.randInt(0, input[0].size() - 1);
	        	if(PseudoRandom.randDouble(0, 1) < 0.0) {
	        		labelSet.add(input[0].get(rdInt));
	        	}else {
	        		Solution sol = input[1].get(i);
	                double minDis = distance(sol, input[0].get(0));
	                int minIndex = 0;
	                for(int j=1;j<input[0].size();j++) {
	                    double dis = distance(sol, input[0].get(j));
	                    if(dis < minDis) {
	                        minDis = dis;
	                        minIndex = j;
	                    }
	                }
	                labelSet.add(input[0].get(minIndex));
	        	}
	        }
	        return labelSet;
	}
	
	public SolutionSet[] getClassifiedSet(){
		SolutionSet[] st = new SolutionSet[2];
		NondominatedRanking ranking = new NondominatedRanking(population);
		int size = populationSize/2;
		st[0] = new SolutionSet(size);
		st[1] = new SolutionSet(size);
		int remain = size;
		int index = 0;
		SolutionSet front = null;
		// Obtain the next front
		front = ranking.getSubfront(index);

		while ((remain > 0) && (remain >= front.size())) {
			// Assign crowding distance to individuals
			distance.crowdingDistanceAssignment(front,
					problem_.getNumberOfObjectives());
			// Add the individuals of this front
			for (int k = 0; k < front.size(); k++) {
				st[0].add(front.get(k));
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
				st[0].add(front.get(k));
			} // for
			for (int k = remain; k < front.size(); k++) {
				st[1].add(front.get(k));
			} // for
			remain = 0;
			index++;
		} // if
		
		for(int i=index;i<ranking.getNumberOfSubfronts();i++){
			front = ranking.getSubfront(i);
			for (int k = 0; k < front.size(); k++) {
				st[1].add(front.get(k));
			} // for
		}
		return st;
	}
	
	public double distance(Solution s1, Solution s2) {
		double dis = 0.0;
		double sum1 = s1.getSumValue();
		double sum2 = s2.getSumValue();
		double v1,v2;
		for(int i=0;i<numObjs;i++) {
			v1 = s1.getNormalizedObjective(i)/sum1;
			v2 = s2.getNormalizedObjective(i)/sum2;
			dis += (v1 - v2) * (v1 - v2);
		}
		dis = Math.sqrt(dis);
		return dis;
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
		}
	}

	public void doMutation(double probability, Solution solution)
			throws JMException {
		double rnd, delta1, delta2, mut_pow, deltaq;
		double y, yl, yu, val, xy;
		XReal x = new XReal(solution);
		for (int var = 0; var < solution.numberOfVariables(); var++) {
			if (PseudoRandom.randDouble() <= probability) {
				y = x.getValue(var);
				yl = x.getLowerBound(var);
				yu = x.getUpperBound(var);
				delta1 = (y - yl) / (yu - yl);
				delta2 = (yu - y) / (yu - yl);
				rnd = PseudoRandom.randDouble();
				mut_pow = 1.0 / (20 + 1.0);
				if (rnd <= 0.5) {
					xy = 1.0 - delta1;
					val = 2.0 * rnd + (1.0 - 2.0 * rnd)
							* (Math.pow(xy, (20 + 1.0)));
					deltaq = java.lang.Math.pow(val, mut_pow) - 1.0;
				} else {
					xy = 1.0 - delta2;
					val = 2.0
							* (1.0 - rnd)
							+ 2.0
							* (rnd - 0.5)
							* (java.lang.Math.pow(xy,
							(20 + 1.0)));
					deltaq = 1.0 - (java.lang.Math.pow(val, mut_pow));
				}
				y = y + deltaq * (yu - yl);
				if (y < yl)
					y = yl;
				if (y > yu)
					y = yu;
				x.setValue(var, y);
			}
		} // for
	}
}
