package jmetal.util.ranking;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PairData;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.comparators.RankComparator;

public class NormEnsembleFitnessRanking implements Ranking {
	private SolutionSet solutionSet_;

	private List<SolutionSet> ranking_;

	double[][] lambda_;

	int K_;

	int obj_;

	SolutionSet[] sets_;

	public NormEnsembleFitnessRanking(SolutionSet solutionSet,
			double[][] lambda, int K) {
		this.solutionSet_ = solutionSet;
		this.lambda_ = lambda;

		this.K_ = K;

		this.obj_ = solutionSet.get(0).getNumberOfObjectives();

		sets_ = new SolutionSet[lambda_.length];

		for (int i = 0; i < sets_.length; i++)
			sets_[i] = new SolutionSet();

		doMaximumRanking();
	}

	public void doMaximumRanking() {

		for (int i = 0; i < solutionSet_.size(); i++) {
			Solution sol = solutionSet_.get(i);
			sol.setRank(Integer.MAX_VALUE);

			int indexes[] = findWeightVectors(sol);

			for (int j = 0; j < indexes.length; j++) {
				int k = indexes[j];

				sets_[k].add(sol);
			}
		}

		for (int i = 0; i < sets_.length; i++) {

			for (int j = 0; j < sets_[i].size(); j++)
				setFitness(sets_[i].get(j), lambda_[i]);

			sets_[i].sort(new FitnessComparator());

			for (int j = 0; j < sets_[i].size(); j++) {
				if (j < sets_[i].get(j).getRank())
					sets_[i].get(j).setRank(j);
			}
		}

		solutionSet_.sort(new RankComparator());

		ranking_ = new ArrayList<SolutionSet>();

		int curRank = solutionSet_.get(0).getRank();
		SolutionSet set = new SolutionSet();

		for (int i = 0; i < solutionSet_.size(); i++) {
			Solution sol = solutionSet_.get(i);

			if (sol.getRank() != curRank) {
				ranking_.add(set);
				set = new SolutionSet();
				curRank = sol.getRank();
				set.add(sol);
			} else
				set.add(sol);

		}

		ranking_.add(set);

	}

	void setFitness(Solution individual, double[] lambda) {
		double fitness;
		fitness = 0.0;

		double maxFun = -1.0e+30;

		for (int n = 0; n < obj_; n++) {
			double diff = Math.abs(individual.getNormalizedObjective(n));

			double feval;
			if (lambda[n] == 0) {
				feval = diff / 0.000001;
			} else {
				feval = diff / lambda[n];
			}
			if (feval > maxFun) {
				maxFun = feval;
			}
			
			//fitness += individual.getNormalizedObjective(n)*lambda[n];
		} // for

		fitness = maxFun;

		individual.setFitness(fitness);

	}

	int[] findWeightVectors(Solution sol) {
		int indexes[] = new int[K_];

		PairData[] pairs = new PairData[lambda_.length];
		for (int i = 0; i < lambda_.length; i++) {
			double dist = getPerpendicularDistance(sol, lambda_[i]);
			pairs[i] = new PairData(i, dist);
		}

		Arrays.sort(pairs);

		for (int i = 0; i < K_; i++)
			indexes[i] = pairs[i].getID();

		return indexes;
	}

	double getPerpendicularDistance(Solution sol, double[] ref) {
		double ip = 0;
		double refLenSQ = 0;

		double[] d = new double[2];
		for (int j = 0; j < obj_; j++) {

			ip += sol.getNormalizedObjective(j) * ref[j];
			refLenSQ += (ref[j] * ref[j]);
		}
		refLenSQ = Math.sqrt(refLenSQ);

		d[0] = Math.abs(ip) / refLenSQ;

		d[1] = 0;

		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			d[1] += (sol.getNormalizedObjective(i) - d[0] * (ref[i] / refLenSQ))
					* (sol.getNormalizedObjective(i) - d[0]
							* (ref[i] / refLenSQ));
		}
		d[1] = Math.sqrt(d[1]);

		return (d[0]+10.0*d[1]);
	}

	public SolutionSet getSubfront(int layer) {
		return ranking_.get(layer);
	}

	public int getNumberOfSubfronts() {
		return ranking_.size();
	}

}
