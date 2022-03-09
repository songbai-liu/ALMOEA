package jmetal.metaheuristics.moead;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.Permutation;

public class DecompositionBasedSelection {
    private SolutionSet population_;
    private SolutionSet winner_;
    private SolutionSet loser_;
    private int numbObj_;
    private int size_;
    private double[][] lamda_;

    public DecompositionBasedSelection(SolutionSet population, int number_obj, double[][] lamda){
        population_ = population;
        numbObj_ = number_obj;
        size_ = lamda.length;
        winner_ = new SolutionSet(size_);
        loser_ = new SolutionSet(population.size() - size_);
        lamda_ = lamda;
    }

    public SolutionSet environmentalSelection(){
        SolutionSet[] solSet = new SolutionSet[2];
        SolutionSet[] clusters = new SolutionSet[size_];
        // Create the solutionSet union of solutionSet and offSpring
        for(int c=0; c<size_; c++){
            clusters[c] = new SolutionSet();
        }
        //每个解都跟与它距离最近的参考向量绑定
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
                    minfit = computeTCH(clusters[i].get(0), lamda_[i]);
                    minID = 0;
                    for(int c=1;c<clusters[i].size();c++){
                        fit = computeTCH(clusters[i].get(c), lamda_[i]);
                        if(fit < minfit){
                            minfit = fit;
                            minID = c;
                        }
                    }
                    //每个cluster里与参考向量加权和最小的作为winner
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
        return solSet[0];
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

    public double computeTCH(Solution individual, double[] lambda){
        double fitness;
        fitness = 0.0;

        double maxFun = -1.0e+30;

        for (int n = 0; n < numbObj_; n++) {
            double diff = Math.abs(individual.getNormalizedObjective(n));

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
        return fitness;
    }

    public double computeWS1(Solution sol, double[] lamada){
        double fitness = 0;
        for(int m=0;m<numbObj_;m++){
            fitness += sol.getNormalizedObjective(m)*lamada[m];
        }
        return fitness;
    }

}
