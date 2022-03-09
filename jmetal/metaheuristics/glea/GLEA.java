package jmetal.metaheuristics.glea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.moead.Utils;
import jmetal.util.JMException;
import jmetal.util.Permutation;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.ranking.ThetaRanking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class GLEA extends Algorithm{
	
	private int populationSize_;//population size
	private int numObj_; //number of objectives
	private int numVar_; //number of variables
	private int interval_;//the sliding interval
	
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
	private Operator crossoverOperator_;
	private Operator mutationOperator_;	
	private int groupSize_; //number of groups
	private int[][] group_;
	private int[][] bGroup_;
	private int[] runningGroup1;
	private int[] allVariables;
	double[] zideal_;
	double[] nadir_;
	double[][] lamada1_;
	double[][] lamada2_;
	int div1_,div2_;
	int div3_,div4_;
	double[] UB_;
	double[] LB_;
	public GLEA(Problem problem) {
		super(problem);
		numObj_ = problem.getNumberOfObjectives();
		numVar_ = problem.getNumberOfVariables();
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		int g = 0;
		while(generations_ <= maxGenerations_){
			if((generations_) % (interval_) == 0) { randomGrouping1(); } 
			//reproduction();
			reproduction(generations_);
			environmentalSelection();
			generations_++;
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population_);
		return final_ranking.getSubfront(0);
	}
	
	public int[][] broadeningGrouping(int[][] group){
		int numOfGroup = group.length;
		int[][] groups = new int[numOfGroup][];
		int T_ = 0;
		groups[0] = group[0];
		while(T_ < numOfGroup-1){
			T_++;
			int[] groupA = groups[T_-1];
			int[] groupB = group[T_];
			groups[T_] = new int[groupA.length + groupB.length];
			System.arraycopy(groupA, 0, groups[T_], 0, groupA.length);
			System.arraycopy(groupB, 0, groups[T_], groupA.length, groupB.length);
		}
		return groups;
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
		crossoverOperator_ = operators_.get("crossover");
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
		allVariables = new int[numVar_];
		for (int i = 0; i < numVar_; i++) {
			UB_[i] = problem_.getUpperLimit(i);
			LB_[i] = problem_.getLowerLimit(i);
			allVariables[i] = i;
		} // for
		groupSize_ = 5 + numVar_/100;
		//groupSize_ = 5;
		//randomGrouping();
		//linearGrouping();
		//bGroup_ = broadeningGrouping(group_);
		interval_ = 100;
	}
	
	public void reproduction() throws JMException{
		estimateIdealPoint(population_);
		estimateNadirPoint(population_);
		normalizedPopulation(population_);
		SolutionSet[] solSet = new PopulationClassification(population_,numObj_,lamada2_).classification();
		//Group-based Learning for generating the offsprings
		groupWeightLearning(solSet[0]);
		//groupNormalLearning(solSet[0]);
		groupBroadLearning(solSet[0], solSet[1]);
	}//for reproduction
	
	public void reproduction(int gen) throws JMException{
		if(gen < 0.0*maxGenerations_) {
			groupNormalLearning(population_);
		}else if(gen <= 1.0*maxGenerations_) {
			estimateIdealPoint(population_);
			estimateNadirPoint(population_);
			normalizedPopulation(population_);
			SolutionSet[] solSet = new PopulationClassification(population_,numObj_,lamada2_).classification();
			//Group-based Learning for generating the offsprings
			groupWeightLearning(solSet[0]);
			//groupNormalLearning(solSet[0]);
			groupBroadLearning(solSet[0], solSet[1]);
		}else {
			groupNormalLearning(population_);
		}
		
	}//for reproduction
	
	public void groupBroadLearning(SolutionSet s1, SolutionSet s2) throws JMException {
		double r1;
		for (int i = 0; i < s2.size(); i++) {
			Solution child = new Solution(s2.get(i));
			XReal xChild = new XReal(child);
			//int rd1 = PseudoRandom.randInt(0, s1.size()-1);
			for(int j=0;j<groupSize_;j++){
				double val;
				r1 = PseudoRandom.randDouble(0.0, 1.0);
				int rd1 = PseudoRandom.randInt(0, s1.size()-1);
				XReal xParent = new XReal(s1.get(rd1));
				for(int m=0;m<group_[j].length;m++){
					val = xChild.getValue(group_[j][m]) + r1*(xParent.getValue(group_[j][m]) - xChild.getValue(group_[j][m]));
					if (val < xChild.getLowerBound(group_[j][m]))
						val = xChild.getLowerBound(group_[j][m]);
					if (val > xChild.getUpperBound(group_[j][m]))
						val = xChild.getUpperBound(group_[j][m]);
					xChild.setValue(group_[j][m], val);
				}
			}
			
			mutationOperator_.execute(child);
			problem_.evaluate(child);
			offspringPopulation_.add(child);
		}
	}
	
	public void groupWeightLearning(SolutionSet s1) throws JMException {
		double [][] weight =  new double[3][groupSize_];
		double r1;
		double alph = 0.3;
		double beta = 3.0;
		for (int i = 0; i < s1.size(); i++) {
			Solution child = new Solution(s1.get(i));
			XReal xChild = new XReal(child);
			r1 = PseudoRandom.randDouble(0.0, 1.0);
			int rd1 = PseudoRandom.randInt(0, s1.size()-1);
			while(rd1 == i) {
				rd1 = PseudoRandom.randInt(0, s1.size()-1);
			}
			int rd2 = PseudoRandom.randInt(0, s1.size()-1);
			while(rd2 == i || rd2 == rd1) {
				rd2 = PseudoRandom.randInt(0, s1.size()-1);
			}
			XReal ss1 = new XReal(s1.get(rd1));
			XReal ss2 = new XReal(s1.get(rd2));
			
			for(int j=0;j<groupSize_;j++){
				weight[0][j] = 0;
				weight[1][j] = 0;
				weight[2][j] = 0;
				double sum = 0;
				for(int g=0;g<group_[j].length;g++){//transfer to weighted space
					sum += (xChild.getValue(group_[j][g]) - LB_[group_[j][g]])/(UB_[group_[j][g]] - LB_[group_[j][g]]);
					weight[0][j] += (xChild.getValue(group_[j][g]) - LB_[group_[j][g]])/(UB_[group_[j][g]] - LB_[group_[j][g]])/group_[j].length;
					weight[1][j] += (ss1.getValue(group_[j][g]) - LB_[group_[j][g]])/(UB_[group_[j][g]] - LB_[group_[j][g]])/group_[j].length;
					weight[2][j] += (ss2.getValue(group_[j][g]) - LB_[group_[j][g]])/(UB_[group_[j][g]] - LB_[group_[j][g]])/group_[j].length;
				}
				double w = weight[0][j] + r1*(weight[1][j] - weight[2][j]);//learning in the weighted space, i.e., weighted optimization
				if(w < alph*weight[0][j]) {
					w = alph*weight[0][j]; 
				}else if(w > beta*weight[0][j]){ 
					w = beta*weight[0][j]; 
				}
				double val;
				for(int m=0;m<group_[j].length;m++){//Recover to the original variable space
					val = w*group_[j].length*(xChild.getValue(group_[j][m]) - LB_[group_[j][m]])/(UB_[group_[j][m]] - LB_[group_[j][m]])/sum;
					val = val * (UB_[group_[j][m]] - LB_[group_[j][m]]) + LB_[group_[j][m]];
					if (val < xChild.getLowerBound(group_[j][m]))
						val = xChild.getLowerBound(group_[j][m]);
					if (val > xChild.getUpperBound(group_[j][m]))
						val = xChild.getUpperBound(group_[j][m]);
					xChild.setValue(group_[j][m], val);
				}
			}
			mutationOperator_.execute(child);
			problem_.evaluate(child);
			offspringPopulation_.add(child);
		}
	}
	
	public void groupNormalLearning(SolutionSet s1) throws JMException {
		double r1, r2, r3;
		for (int i = 0; i < s1.size(); i++) {
			Solution child = new Solution(s1.get(i));
			XReal xChild = new XReal(child);
			r1 = PseudoRandom.randDouble(0.0, 1.0); 
			int rd1 = PseudoRandom.randInt(0, s1.size()-1);
			while(rd1 == i) {
				rd1 = PseudoRandom.randInt(0, s1.size()-1);
			}
			int rd2 = PseudoRandom.randInt(0, s1.size()-1);
			while(rd2 == rd1 || rd2 == i) {rd2 = PseudoRandom.randInt(0, s1.size()-1); }
			
			XReal xParent1;
			XReal xParent2;
			xParent1 = new XReal(s1.get(rd1));
			xParent2 = new XReal(s1.get(rd2));
			
			for(int j=0;j<groupSize_;j++){
				double val;
				for(int m=0;m<group_[j].length;m++){
					val = xChild.getValue(group_[j][m]) + r1*(xParent1.getValue(group_[j][m]) - xParent2.getValue(group_[j][m]));
					if (val < xChild.getLowerBound(group_[j][m]))
						val = xChild.getLowerBound(group_[j][m]);
					if (val > xChild.getUpperBound(group_[j][m]))
						val = xChild.getUpperBound(group_[j][m]);
					xChild.setValue(group_[j][m], val);
				}
			}
			mutationOperator_.execute(child);
			problem_.evaluate(child);
			offspringPopulation_.add(child);
		}
	}
	
	public void randomGrouping(){
		int remain = numVar_%groupSize_;
		int divisor = numVar_/groupSize_;
		group_ = new int[groupSize_][];
		int iter = 0;
		while(iter < groupSize_) {
			if(iter < remain) group_[iter] = new int[divisor+1];
			else group_[iter] = new int[divisor];
			iter++;
		}
		
		int[] permutation = new int[numVar_];
		Utils.randomPermutation(permutation, numVar_);
		int t = 0;
		for(int g=0;g<groupSize_;g++){
			for(int m=0;m<divisor;m++){
				group_[g][m] = permutation[t];
				t++;
			}
		}
		
		for(int g=0;g<remain;g++){
			group_[g][divisor] = permutation[t];
			t++;
		}
	}
	
	public void randomGrouping1(){
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