package jmetal.problems.LMF;

import java.util.ArrayList;
import java.util.List;

public class GroupingStrategy {
	
	/*
	 * niv indicates the number of independent variables of each group
	 * nov indicates the number of overlap variables between two neighbor groups
	 * nsv indicates the number of shared variables among all groups
	 * nvg indicates the number of variable groups
	 * startpos indicates the start position for variable grouping
	*/
	public static int[][] overlapGrouping(int[] niv,int[] nov,int nsv,int nvg,int startPos){
		int[][] group = new int[nvg][];
		int pointer = startPos;
		int[] len = new int[nvg];
		for(int i=0; i<nvg; i++){//assign the independent variables to the last part of each group
			len[i] = niv[i] + nov[i] + nsv;
			group[i] = new int[len[i]];
			for(int j=nov[i]+nsv; j<len[i]; j++){
				group[i][j] = pointer;
				pointer = pointer+1;
			}
		}
		
		for(int i=0; i<nvg; i++){//assign the overlap variables to the first part of each group
			for(int j=0; j<nov[i]; j++){
				if(i == 0)
					group[i][j] = group[nvg-1][len[nvg-1]-1-j];
		        else
		        	group[i][j] = group[i-1][len[i-1]-j-1];
			}
			
			for(int j=0; j<nsv; j++){
				
				
			}
		}
		
		for(int j=0; j<nsv; j++){//assign the shared variables to the middle part of each group
			for(int i=0; i<nvg; i++){
				group[i][j+nov[i]] = pointer;
			}
			pointer = pointer+1; 
		}
		return group;
	}
	
	public static int[][] uniformDeepGrouping(int[] group, int num){
	    int size = group.length;
	    int remain = size % num;
	    int divisor = size/num;
	    int[][] deepGroup = new int[num][];
	    int iter = 0;
	    while(iter < num) {
	    	if(iter < remain) {
	    		deepGroup[iter] = new int[divisor+1];
	    	}else {
	    		deepGroup[iter] = new int[divisor];
	    	}
	    	iter++;
	    }
	    int t = 0;
	    for(int i=0;i<num;i++){
	    	for(int j=0;j<divisor;j++){
	    		deepGroup[i][j] = group[t];
	    		t++;
	    	}
	    }
	    
	    for(int i=0;i<remain;i++) {
	    	deepGroup[i][divisor] = group[t];
	    	t++;
	    }
	    return deepGroup;    
	}
	
	//The first term of the arithmetic sequence is a, and the common difference is d
	public static int[][] deepGrouping(int[] group,int a,int d){
		int span = a;
	    int remain = group.length;
	    int t = 0;
	    List<Integer[]> dg = new ArrayList<Integer[]>();
	    while(remain > span + a){
	    	Integer[] g = new Integer[span];
	    	for(int i=0; i<span; i++){
	    		g[i] = group[t];
		        t = t + 1;
	    	}
	    	dg.add(g);
		    remain = remain - span;
		    span = span + d;
	    }
	        
	    if(remain > 0){
	    	Integer[] g = new Integer[remain];
	    	for(int i=0;i<remain;i++){
	    		g[i] = group[t];
		        t = t + 1;
	    	}
	    	dg.add(g);
	    }
	    
	    int size = dg.size();
	    int[][] deepGroup = new int[size][];
	    for(int i=0;i<size;i++){
	    	deepGroup[i] = new int[dg.get(i).length];
	    	for(int j=0;j<dg.get(i).length;j++){
	    		deepGroup[i][j] = dg.get(i)[j];
	    	}
	    }
	    
	    return deepGroup;    
	}
	
	public static double[] getWeightsLogisticMap(double c1, double r, int num) {
		double[] w = new double[num];
		w[0] = c1;
		double sum = w[0];
		for(int i = 0; i<num-1; i++) {
			w[i+1] = r*w[i]*(1-w[i]);
			sum += w[i+1];
		}
		
		for(int i = 0; i<num; i++) {
			w[i] = w[i]/sum;
		}
		return w;
	}
}
