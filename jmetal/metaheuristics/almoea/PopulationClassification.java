package jmetal.metaheuristics.almoea;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.Permutation;

public class PopulationClassification {
	private SolutionSet population_;
	private SolutionSet winner_;
	private SolutionSet loser_;
	private int numbObj_;
	private int size_;
	private double[][] lamda_;
	
	public PopulationClassification(SolutionSet population, int number_obj, double[][] lamda){
		population_ = population;
		numbObj_ = number_obj;
		size_ = lamda.length;
		winner_ = new SolutionSet(size_);
		loser_ = new SolutionSet(population.size() - size_);
		lamda_ = lamda;
	}
	
	public SolutionSet[] classification(){
		SolutionSet[] solSet = new SolutionSet[2];
		SolutionSet[] clusters = new SolutionSet[size_];
		// Create the solutionSet union of solutionSet and offSpring
		for(int c=0; c<size_; c++){
			clusters[c] = new SolutionSet();
		}
		//æ¯?ä¸ªè§£éƒ½è·Ÿä¸Žå®ƒè·?ç¦»æœ€è¿‘çš„å?‚è€ƒå?‘é‡?ç»‘å®š
		associating(population_, clusters);
		int remain = size_;
		SolutionSet candidateSet = new SolutionSet();
		while(remain > 0){
			double minfit;
			int minID;
			double fit;
			candidateSet.clear();
			for(int i=0;i<size_;i++){
				if(clusters[i].size() > 0){
					minfit = computeWS1(clusters[i].get(0), lamda_[i]);
					minID = 0;
					for(int c=1;c<clusters[i].size();c++){
						fit = computeWS1(clusters[i].get(c), lamda_[i]);
						if(fit < minfit){
							minfit = fit;
							minID = c;
						}
					}
					//æ¯?ä¸ªclusteré‡Œä¸Žå?‚è€ƒå?‘é‡?åŠ æ?ƒå’Œæœ€å°?çš„ä½œä¸ºwinner
					candidateSet.add(clusters[i].get(minID));
					clusters[i].remove(minID);
				}
			}
			if(remain >= candidateSet.size()){
				for(int r=0;r<candidateSet.size();r++){
					winner_.add(candidateSet.get(r));
				}
				remain = remain - candidateSet.size();
				if(remain == 0){
					break;
				}
			}else{
				int[] perm = new Permutation().intPermutation(candidateSet.size());
				for (int k = 0; k < remain; k++) {
					winner_.add(candidateSet.get(perm[k]));
				}
				for (int k = remain; k < candidateSet.size(); k++) {
					loser_.add(candidateSet.get(perm[k]));
				}
				break;
			}	
		}//while	
		
		for(int i=0;i<size_;i++){
			if(clusters[i].size() > 0){
				for(int j=0;j<clusters[i].size();j++){
					loser_.add(clusters[i].get(j));
				}
			}
		}
		
		solSet[0] = winner_;
		solSet[1] = loser_;
		return solSet;
	}
	
	public void associating(SolutionSet sols, SolutionSet[] clusters){
		for(int i=0; i<sols.size(); i++){
			Solution sol = sols.get(i);
			double dis_min = computeD2(sol, lamda_[0]);
			int id_min = 0;
			for(int j=1; j<lamda_.length;j++){
				double dis = computeD2(sol, lamda_[j]);
				if(dis < dis_min){
					dis_min = dis;
					id_min = j;
				}
			}
			clusters[id_min].add(sol);
		}
	}
	
	public double computeD2(Solution sol, double[] ref){
		double ip = 0;
		double refLenSQ = 0;
		double norm = 0.0;
		//double distance = 0.0;
		double[] d = new double[2];
		for (int j = 0; j < numbObj_; j++) {
			ip += sol.getNormalizedObjective(j) * ref[j];
			refLenSQ += (ref[j] * ref[j]);
			norm += sol.getNormalizedObjective(j) * sol.getNormalizedObjective(j);
		}
		refLenSQ = Math.sqrt(refLenSQ);
        norm = Math.sqrt(norm);
		d[0] = Math.abs(ip) / refLenSQ;

		d[1] = 0;
	    d[1] = norm*norm - d[0]*d[0];
		d[1] = Math.sqrt(d[1]);
		return d[1];
	}
	
	public double computeWS1(Solution sol, double[] lamada){
		double fitness = 0;
		for(int m=0;m<numbObj_;m++){
			fitness += sol.getNormalizedObjective(m)*lamada[m];
		}
		return fitness;
	}

}
