package jmetal.metaheuristics.almoea.BaseMLP;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;

public class IHierarchicalClustering {
List<SolutionSet> list = new <SolutionSet>ArrayList();
	
	public IHierarchicalClustering(List list){
		this.list = list;
	}
	public List<SolutionSet> clusteringAnalysis(int clusteringSize){
		int size = list.size();
		double minDistance = computeDistanceOnUHP(list.get(0).getCentroid(),list.get(1).getCentroid());
		
		double distance = 0.0;
		double[] minDistances = new double[size];
		int[] minIndexs = new int[size];
		int[] index = new int[4];
		index[2] = 0;
		index[3] = 1;
		for(int i=0; i<size; i++){
			double min;
			int rd = PseudoRandom.randInt(0, size-1);
			while(rd == i) {
				rd = PseudoRandom.randInt(0, size-1);
			}
			min = computeDistanceOnUHP(list.get(i).getCentroid(),list.get(rd).getCentroid());
			index[0] = i;
			index[1] = rd;
			
			for(int j=0;j<size;j++){
				if(i != j){
					distance = computeDistanceOnUHP(list.get(i).getCentroid(),list.get(j).getCentroid());
					if(min > distance){
						min = distance;
						index[0] = i;
						index[1] = j;
					}
				}
			}
			minDistances[i] =  min;
			minIndexs[i] = index[1];
			if(minDistance > min){
				minDistance = min;
				index[2] = index[0];
				index[3] = index[1];
			}
		}
		if(index[2] == index[3]) {
			index[2] = 0;
			index[3] = 1;
			
		}
		while(size > clusteringSize){
			SolutionSet sols = (list.get(index[2]).union(list.get(index[3])));
			list.get(index[2]).setRemove(true);	
			list.remove(index[3]);
			list.add(index[3], sols);
			
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[2] && !list.get(i).isRemove()){
					double min;
		    		int sb = PseudoRandom.randInt(0, size-1);
		    		while(list.get(sb).isRemove() || sb == i) {
		    			sb = PseudoRandom.randInt(0, size-1);
		    		}
		    		min = computeDistanceOnUHP(list.get(i).getCentroid(),list.get(sb).getCentroid());
		    		
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j){
		    				ss = computeDistanceOnUHP(list.get(i).getCentroid(),list.get(j).getCentroid());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minDistances[i] = min;
			    	minIndexs[i] = sb;
				}//if
			}//for
			
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[3] && !list.get(i).isRemove()){
					double min;
		    		int sb = PseudoRandom.randInt(0, size-1);
		    		while(list.get(sb).isRemove() || sb == i) {
		    			sb = PseudoRandom.randInt(0, size-1);
		    		}
		    		min = computeDistanceOnUHP(list.get(i).getCentroid(),list.get(sb).getCentroid());
		    		
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j){
		    				ss = computeDistanceOnUHP(list.get(i).getCentroid(),list.get(j).getCentroid());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minDistances[i] = min;
			    	minIndexs[i] = sb;
				}//if
			}//for
			
			int srd = PseudoRandom.randInt(0, size-1);
			while(list.get(srd).isRemove()) {
				srd = PseudoRandom.randInt(0, size-1);
			}
			double sAngle = minDistances[srd];
			index[2] = srd;
			index[3] = minIndexs[srd];
			
			for(int k=0;k<list.size();k++){
				if(!list.get(k).isRemove()){
					if(sAngle > minDistances[k]){
						sAngle = minDistances[k];
						index[2] = k;
						index[3] = minIndexs[k]; 	
					}
				}
			}
	
			size--;
		}//while
	
		Iterator<SolutionSet> iterator = list.iterator();
		while(iterator.hasNext()){
			if(iterator.next().isRemove()){
				iterator.remove();
			}
		}
		return this.list;
	}
	
	public double computeDistanceOnUHP(Solution so1, Solution so2){
		double distance = 0.0;
		double innerProduc = 0.0;
		double v1, v2;
		int numObj = list.get(0).get(0).getNumberOfObjectives();
		for(int i=0; i<numObj; i++){
			v1 = (so1.getNormalizedObjective(i)-1)/(so1.getSumValue()-numObj);
			v2 = (so2.getNormalizedObjective(i)-1)/(so2.getSumValue()-numObj);
			innerProduc += Math.pow((v1-v2), 2);
		}
		distance = Math.sqrt(innerProduc);
		return distance;
	}
}
