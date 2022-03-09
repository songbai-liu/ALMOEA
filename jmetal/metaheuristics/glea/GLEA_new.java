package jmetal.metaheuristics.glea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.moead.Utils;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class GLEA_new extends Algorithm{
	
	private int populationSize_;//population size
	private int numObj_; //number of objectives
	private int numVar_; //number of variables
	
	/**
	 * Stores the population
	 */
	private SolutionSet population_;
	private SolutionSet offspringPopulation_;
	private SolutionSet unionPopulation_;	
	private int generations_;
	private int maxGenerations_;	
	/**
	 * Operators
	 */
	private Operator mutationOperator_;	
	private int groupSize_; //number of groups
	private int[][] group_;
	double[] zideal_;
	double[] nadir_;
	double[][] lamada1_;
	double[][] lamada2_;
	int div1_,div2_;
	int div3_,div4_;
	double[] UB_;
	double[] LB_;
	
	public GLEA_new(Problem problem) {
		super(problem);
		numObj_ = problem.getNumberOfObjectives();
		numVar_ = problem.getNumberOfVariables();
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		while(generations_ <= maxGenerations_){
			if((generations_) % (100) == 0) { randomGrouping(); }
			reproduction();
			environmentalSelection();	
			generations_++;
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population_);
		return final_ranking.getSubfront(0);
	}
	
	
	
	public void initialization() throws JMException, ClassNotFoundException{
		generations_ = 0;
		populationSize_ = ((Integer) getInputParameter("populationSize")).intValue();
		maxGenerations_ = ((Integer) getInputParameter("maxGenerations")).intValue();
		div1_ = ((Integer) getInputParameter("div1")).intValue();
		div2_ = ((Integer) getInputParameter("div2")).intValue();
		div3_ = ((Integer) getInputParameter("div3")).intValue();
		div4_ = ((Integer) getInputParameter("div4")).intValue();
		population_ = new SolutionSet(populationSize_);
		offspringPopulation_ = new SolutionSet(populationSize_);
		unionPopulation_ = new SolutionSet(2*populationSize_);
		//Read the operators
		mutationOperator_ = operators_.get("mutation");
		//Create the initial population
		Solution newSolution;
		for (int i = 0; i < populationSize_; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population_.add(newSolution);
		} // for
		zideal_ = new double[problem_.getNumberOfObjectives()];
		nadir_ = new double[problem_.getNumberOfObjectives()];
		/* generate two-layer weight vectors */
		VectorGenerator vg;
		vg = new TwoLevelWeightVectorGenerator(div1_, div2_,numObj_);
		lamada1_ = vg.getVectors();
		vg = new TwoLevelWeightVectorGenerator(div3_, div4_,numObj_);
		lamada2_ = vg.getVectors();
		UB_ = new double[numVar_];
		LB_ = new double[numVar_];
		for (int i = 0; i < numVar_; i++) {
			UB_[i] = problem_.getUpperLimit(i);
			LB_[i] = problem_.getLowerLimit(i);
		} // for
		//groupSize_ = 5;
		groupSize_ = 4 + numVar_/100;
		randomGrouping();
	}
	
	public void reproduction() throws JMException{
		double learningRate;
		Solution[] parents = new Solution[3];
		for (int i = 0; i < population_.size(); i++) {
			// obtain parents
			matingSelection(parents,i);
			Solution child = new Solution(parents[0]);
			learningRate = PseudoRandom.randDouble(0.0, 1.0);
			//encode
			double[][] encodedParents = encode(parents);
            //search in the compressed space
			double[] newEncode = searchInTransferedSpace(encodedParents,learningRate);
			//decode
			decode(newEncode,child);
			//evaluation of the child
			mutationOperator_.execute(child);
			problem_.evaluate(child);
			offspringPopulation_.add(child);
		}
	}//for reproduction
	
	public void matingSelection(Solution[] parents,int i) {
		parents[0] = population_.get(i);
		int rdInt1 = PseudoRandom.randInt(0, populationSize_-1);
		while(rdInt1 == i) {
			rdInt1 = PseudoRandom.randInt(0, populationSize_-1);
		}
		int rdInt2 = PseudoRandom.randInt(0, populationSize_-1);
		while(rdInt2 == i || rdInt2 == rdInt1) {
			rdInt2 = PseudoRandom.randInt(0, populationSize_-1);
		}
		parents[1] = population_.get(rdInt1);
		parents[2] = population_.get(rdInt2);
	}
	
	public double[][] encode(Solution[] parents) throws JMException {
		int len = parents.length;
		double[][] encodedParents = new double[len][groupSize_];
		for(int l=0;l<len;l++) {
			double[] encodedSolution = new double[numVar_];
			encodedSolution = problem_.normalizeVariableSet(parents[l]);
			for(int g=0;g<groupSize_;g++) {
				encodedParents[l][g] = 0;
				for(int var=0;var<group_[g].length;var++) {
					encodedParents[l][g] += encodedSolution[group_[g][var]];
				}
				encodedParents[l][g] = encodedParents[l][g]/group_[g].length;
			}
		}
		return encodedParents; 
	}
	
	public double[] searchInTransferedSpace(double[][] encodedParents, double learningRate) {
		double[] newEncode = new double[groupSize_];
		for(int j=0;j<groupSize_;j++) {
			newEncode[j] = encodedParents[0][j] + learningRate*(encodedParents[1][j] - encodedParents[2][j]);
			if(newEncode[j] < encodedParents[0][j]*0.3) { newEncode[j] = encodedParents[0][j]*0.3; } 
			if(newEncode[j] > encodedParents[0][j]*3.0) { newEncode[j] = encodedParents[0][j]*3.0; }
		}
		return newEncode;
	}
	
	public void decode(double[] newEncode, Solution child) throws JMException {
		XReal offspring = new XReal(child);
		for(int g=0;g<groupSize_;g++) {
			double sum = 0;
			for(int var=0;var<group_[g].length;var++) {
				int cVar = group_[g][var];
				sum += (offspring.getValue(cVar)-LB_[cVar])/(UB_[cVar] - LB_[cVar]);
			}
			for(int var=0;var<group_[g].length;var++) {
				int cVar = group_[g][var];
				double normalizedValue = (offspring.getValue(cVar)-LB_[cVar])/(UB_[cVar] - LB_[cVar]);
				double decodedValue = (normalizedValue/sum) * newEncode[g] * group_[g].length;
				decodedValue = decodedValue*(UB_[cVar] - LB_[cVar]) + LB_[cVar];
				if(decodedValue < LB_[cVar]) {
					decodedValue = LB_[cVar];
				}
				if(decodedValue > UB_[cVar]) {
					decodedValue = UB_[cVar];
				}
				offspring.setValue(cVar, decodedValue);
			}
		}
	}
	
	public void randomGrouping(){
		int gSize_ = numVar_/groupSize_;
		group_ = new int[groupSize_][];
		for(int g=0;g<groupSize_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVar_-(groupSize_-1)*gSize_;//the variable size of the last group
		group_[groupSize_-1] = new int[lSize_];
		int[] permutation = new int[numVar_];
		Utils.randomPermutation(permutation, numVar_);
		int t = 0;
		for(int g=0;g<groupSize_-1;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = permutation[t];
				t++;
			}
		}
		//assign variable to the last group
		for(int m=0;m<lSize_;m++){
			group_[groupSize_-1][m] = permutation[t];
			t++;
		}
	}
	
	public void linearGrouping(){
		int gSize_ = numVar_/groupSize_;
		group_ = new int[groupSize_][];
		for(int g=0;g<groupSize_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVar_-(groupSize_-1)*gSize_;//the variable size of the last group
		group_[groupSize_-1] = new int[lSize_];
		
		int t = 0;
		for(int g=0;g<groupSize_-1;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = t;
				t++;
			}
		}
		//assign variable to the last group
		for(int m=0;m<lSize_;m++){
			group_[groupSize_-1][m] = t;
			t++;
		}
	}
	/*
	 * Estimate the Ideal Point 
	 */
	public void estimateIdealPoint(SolutionSet solutionSet){
		for(int i=0; i<numObj_;i++){
			zideal_[i] = 1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) < zideal_[i]){
					zideal_[i] = solutionSet.get(j).getObjective(i);
				}//if
			}//for
		}//for
	}
	
	/*
	 * Estimate the Ideal Point 
	 */
	public void estimateNadirPoint(SolutionSet solutionSet){
		for(int i=0; i<numObj_;i++){
			nadir_[i] = -1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) > nadir_[i]){
					nadir_[i] = solutionSet.get(j).getObjective(i);
				}//if
			}//for
		}//for
	}
	
	public void normalizedPopulation(SolutionSet solutionSet){
		double minNomal = Double.MAX_VALUE;
		int minIndex = 0;
		for(int i=0; i<solutionSet.size();i++){
			Solution sol = solutionSet.get(i);
			double sum = 0.0;
			double normal = 0.0;
			for(int j=0; j<problem_.getNumberOfObjectives();j++){
				//double value = (sol.getObjective(j)-zideal_[j])/(nadir_[j]-zideal_[j]);
				double value = (sol.getObjective(j)-zideal_[j]);
				sol.setNormalizedObjective(j, value);
				sum += value;
				normal += value*value;
			}
			normal = Math.sqrt(normal);
			if(normal < minNomal) {
				minNomal = normal;
				minIndex = i;
			}
			sol.setDistanceToIdealPoint(normal);
			sol.setSumValue(sum);
			sol.setFitness(sum);
		}
	}
	
	public void environmentalSelection(){
		// Create the solutionSet union of solutionSet and offSpring
		unionPopulation_ = ((SolutionSet) population_).union(offspringPopulation_);
		population_.clear();
		offspringPopulation_.clear();
		SolutionSet st = getStSolutionSet(unionPopulation_,populationSize_);
		estimateIdealPoint(st);
		estimateNadirPoint(st);
		normalizedPopulation(st);
		population_ = new PopulationClassification(st,numObj_,lamada1_).classification()[0];
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
}