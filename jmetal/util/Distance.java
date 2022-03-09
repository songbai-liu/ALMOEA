//  Distance.java
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

package jmetal.util;

import java.util.Arrays;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;
import jmetal.util.comparators.NormalizedObjectiveComparator;
import jmetal.util.comparators.ObjectiveComparator;

/**
 * This class implements some utilities for calculating distances
 */
public class Distance {      
    
  /** 
  * Constructor.
  */
  public Distance() {
    //do nothing.
  } // Distance
        
    
  /** 
  * Returns a matrix with distances between solutions in a 
  * <code>SolutionSet</code>.
  * @param solutionSet The <code>SolutionSet</code>.
  * @return a matrix with distances.
  */
  public double [][] distanceMatrix(SolutionSet solutionSet) {
    Solution solutionI, solutionJ;

    //The matrix of distances
    double [][] distance = new double [solutionSet.size()][solutionSet.size()];        
    //-> Calculate the distances
    for (int i = 0; i < solutionSet.size(); i++){
      distance[i][i] = 0.0;
      solutionI = solutionSet.get(i);
      for (int j = i + 1; j < solutionSet.size(); j++){
        solutionJ = solutionSet.get(j);
        distance[i][j] = this.distanceBetweenObjectives(solutionI,solutionJ);                
        distance[j][i] = distance[i][j];            
      } // for
    } // for        
    
    //->Return the matrix of distances
    return distance;
  } // distanceMatrix
    
 /** Returns the minimum distance from a <code>Solution</code> to a 
  * <code>SolutionSet according to the objective values</code>.
  * @param solution The <code>Solution</code>.
  * @param solutionSet The <code>SolutionSet</code>.
  * @return The minimum distance between solution and the set.
 * @throws JMException 
  */  
  public double distanceToSolutionSetInObjectiveSpace(Solution    solution, 
		                                  SolutionSet solutionSet) throws JMException{
    //At start point the distance is the max
    double distance = Double.MAX_VALUE;    
        
    // found the min distance respect to population
    for (int i = 0; i < solutionSet.size();i++){            
      double aux = this.distanceBetweenObjectives(solution,solutionSet.get(i));
      if (aux < distance)
        distance = aux;
    } // for
    
    //->Return the best distance
    return distance;
  } // distanceToSolutionSetinObjectiveSpace
    
  /** Returns the minimum distance from a <code>Solution</code> to a 
   * <code>SolutionSet according to the encodings.variable values</code>.
   * @param solution The <code>Solution</code>.
   * @param solutionSet The <code>SolutionSet</code>.
   * @return The minimum distance between solution and the set.
  * @throws JMException 
   */  
   public double distanceToSolutionSetInSolutionSpace(Solution    solution, 
 		                                  SolutionSet solutionSet) throws JMException{
     //At start point the distance is the max
     double distance = Double.MAX_VALUE;    
         
     // found the min distance respect to population
     for (int i = 0; i < solutionSet.size();i++){            
       double aux = this.distanceBetweenSolutions(solution,solutionSet.get(i));
       if (aux < distance)
         distance = aux;
     } // for
     
     //->Return the best distance
     return distance;
   } // distanceToSolutionSetInSolutionSpace
    
 /** Returns the distance between two solutions in the search space.
  *  @param solutionI The first <code>Solution</code>. 
  *  @param solutionJ The second <code>Solution</code>.
  *  @return the distance between solutions.
 * @throws JMException 
  */
  public double distanceBetweenSolutions(Solution solutionI, Solution solutionJ) 
  throws JMException{                
    double distance = 0.0;
    if ((solutionI.getDecisionVariables() != null) && 
    		(solutionJ.getDecisionVariables() != null)) {
      Variable[] decisionVariableI = solutionI.getDecisionVariables();
      Variable[] decisionVariableJ = solutionJ.getDecisionVariables();    
    
      double diff;    //Auxiliar var
      //-> Calculate the Euclidean distance
      for (int i = 0; i < decisionVariableI.length; i++){
        diff = decisionVariableI[i].getValue() -
               decisionVariableJ[i].getValue();
        distance += Math.pow(diff,2.0);
      } // for    
    }    
    //-> Return the euclidean distance
    return Math.sqrt(distance);
  } // distanceBetweenSolutions
    
 /** Returns the distance between two solutions in objective space.
  *  @param solutionI The first <code>Solution</code>.
  *  @param solutionJ The second <code>Solution</code>.
  *  @return the distance between solutions in objective space.
  */
  public double distanceBetweenObjectives(Solution solutionI, Solution solutionJ){                
    double diff;    //Auxiliar var
    double distance = 0.0;
    //-> Calculate the euclidean distance
    for (int nObj = 0; nObj < solutionI.numberOfObjectives();nObj++){
      diff = solutionI.getObjective(nObj) - solutionJ.getObjective(nObj);
      distance += Math.pow(diff,2.0);           
    } // for   
        
    //Return the euclidean distance
    return Math.sqrt(distance);
  } // distanceBetweenObjectives.

 /** Assigns crowding distances to all solutions in a <code>SolutionSet</code>.
  * @param solutionSet The <code>SolutionSet</code>.
  * @param nObjs Number of objectives.
  */
  public void crowdingDistanceAssignment(SolutionSet solutionSet, int nObjs) {
    int size = solutionSet.size();        
                
    if (size == 0)
      return;
    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double distance;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objetiveMinn = front.get(0).getObjective(i);      
      objetiveMaxn = front.get(front.size()-1).getObjective(i);      
      
      //Set de crowding distance            
      front.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      front.get(size-1).setCrowdingDistance(Double.POSITIVE_INFINITY);                                      
           
      
      for (int j = 1; j < size-1; j++) {
        distance = front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i);                    
        distance = distance / (objetiveMaxn - objetiveMinn);        
        distance += front.get(j).getCrowdingDistance();                
        front.get(j).setCrowdingDistance(distance);   
      } // for
    } // for        
  } // crowdingDistanceAssing            

public void sbcrowdingDistanceAssignment(SolutionSet solutionSet, int nObjs) {
    int size = solutionSet.size();        
                
    if (size == 0)
      return;
    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1.0);
      solutionSet.get(1).setCrowdingDistance(1.0);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    //double objetiveMaxn;
    //double objetiveMinn;
    //double distance;
    double []distance = new double [size];
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
    
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);      
      objectiveMax[i] = front.get(front.size()-1).getObjective(i);      
      
      //Set de crowding distance            
     // front.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
     // front.get(size-1).setCrowdingDistance(Double.POSITIVE_INFINITY);                                      
           
      
      for (int j = 1; j < size-1; j++) {
        distance[j] = front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i);                    
        distance[j] = distance[j] / (objectiveMax[i] - objectiveMin[i]);        
        distance[j] += front.get(j).getCrowdingDistance();                
        front.get(j).setCrowdingDistance(distance[j]); 
      } // for
      
      distance[0] =  2*(front.get(1).getObjective(i) - front.get(0).getObjective(i));
      distance[0] = distance[0] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[0] += front.get(0).getCrowdingDistance(); 
      front.get(0).setCrowdingDistance(distance[0]);   
      distance[0] = 0;
      distance[size-1] = 2*(front.get(size-1).getObjective(i) - front.get( size-2).getObjective(i));
      distance[size-1] = distance[size-1] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[size-1] += front.get(size-1).getCrowdingDistance();  
      front.get(size-1).setCrowdingDistance(distance[size-1]); 
      distance[size-1] = 0;	      
      
    } // for      
    
    for(int n = 0; n<size; n++){
 	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
 		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
 	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
 		   minCrowdistance = front.get(n).getCrowdingDistance(); 
    } 

 //对拥挤距离进行归一化处理
    double sumCrowdingDistance = 0.0;
    double aveCrowdingDistance = 0.0;

    for(int n = 0 ; n < size; n ++){
 	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
// 	   front.get(n).setCrowdingDistance(front.get(n).getCrowdingDistance()/10);    效果不好？所有的拥挤距离的值较小 无法估计种群的多样性。
 	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
// 	   System.out.println("maxCrowdistance="+maxCrowdistance);
// 	   System.out.println(n+"---------"+front.get(n).getCrowdingDistance());
    }
  } // crowdingDistanceAssing  

public void crowdingDistanceAssignment1(SolutionSet solutionSet, int nObjs) {
    int size = solutionSet.size();        
                
    if (size == 0)
      return;
    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double distance = 0.0;
    
    for(int i=0;i<size;i++){
    	for(int j=0;j<i;j++){
    		for(int p=0;p<nObjs;p++){
    			distance += Math.pow(Math.abs(front.get(i).getObjective(p)-front.get(j).getObjective(p)), nObjs);
    		}
    	}
    	for(int k=i+1;k<size;k++){
    		for(int q=0;q<nObjs;q++){
    			distance += Math.pow(Math.abs(front.get(i).getObjective(q)-front.get(k).getObjective(q)), nObjs);
    			}
    	}
    	distance = Math.pow(distance, 1/nObjs);
    	distance = distance/size;
    	front.get(i).setCrowdingDistance(distance);
    }
                     
  } // crowdingDistanceAssing 

public void crowdingDistanceAssignment2(SolutionSet solutionSet, int nObjs) {
    int size = solutionSet.size();        
                
    if (size == 0)
      return;
    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double[][] distance = new double[size][size-1];
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
    
    for(int i=0;i<size;i++){
    	for(int j=0;j<i;j++){
    		for(int p=0;p<nObjs;p++){
    			distance[i][j] += Math.pow(Math.abs(front.get(i).getObjective(p)-front.get(j).getObjective(p)), 2);
    			distance[i][j] = Math.sqrt(distance[i][j]);
    		}
    		distance[i][j] = 1/distance[i][j];
    	}
    	for(int k=i+1;k<size;k++){
    		for(int q=0;q<nObjs;q++){
    			distance[i][k-1] += Math.pow(Math.abs(front.get(i).getObjective(q)-front.get(k).getObjective(q)), 2);
    			distance[i][k-1] = Math.sqrt(distance[i][k-1]);
    			}
    		distance[i][k-1] = 1/distance[i][k-1];
    	}
    	double sumDistance = 0.0;
    	for(int t=0;t<distance[i].length;t++){
    		sumDistance += distance[i][t];
    	}
    	front.get(i).setCrowdingDistance(distance[i].length/sumDistance);
    }
    
    for(int n = 0; n<size; n++){
  	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
  		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
  	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
  		   minCrowdistance = front.get(n).getCrowdingDistance(); 
     }
    
    for(int n = 0 ; n < size; n ++){
  	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
  	}
                     
  } // crowdingDistanceAssing

public void crowdingDistanceAssignment5(SolutionSet solutionSet, int nObjs) {
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] objective_SDE  = new double [size];
    double [] objective_SDE1  = new double [size];
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值
//      objetiveMinn = front.get(0).getObjective(i);      
//      objetiveMaxn = front.get(front.size()-1).getObjective(i); 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
      
      double maxDistance = -1.0e+30;
	   double minDistance = 1.0e+30;  
      
      //Set de crowding distance   这里将边界点的拥挤距离设置为无穷大是不合理的，因为要考虑到比边界点可能特别的差                                                 
      for (int j = 1; j < size-1; j++) {
          distance[j] =front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i) ;                    
          distance[j] = distance[j] /(objectiveMax[i] -  objectiveMin[i]) ;  			   
          distance[j] += front.get(j).getCrowdingDistance();                
          front.get(j).setCrowdingDistance(distance[j]);   
        } // for
      distance[0] =  2*(front.get(1).getObjective(i) - front.get(0).getObjective(i));
      distance[0] = distance[0] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[0] += front.get(0).getCrowdingDistance(); 
      front.get(0).setCrowdingDistance(distance[0]);   
      distance[0] = 0;
      distance[size-1] = 2*(front.get(size-1).getObjective(i) - front.get( size-2).getObjective(i));
      distance[size-1] = distance[size-1] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[size-1] += front.get(size-1).getCrowdingDistance();  
      front.get(size-1).setCrowdingDistance(distance[size-1]); 
      distance[size-1] = 0;	      
    } // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
//	   front.get(n).setCrowdingDistance(front.get(n).getCrowdingDistance()/10);    效果不好？所有的拥挤距离的值较小 无法估计种群的多样性。
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
//	   System.out.println("maxCrowdistance="+maxCrowdistance);
//	   System.out.println(n+"---------"+front.get(n).getCrowdingDistance());
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;
//   for(int k = 0; k < front.size(); k++){
//	   if(front.get(k).getCrowdingDistance() >= (1 + 0.5) *aveCrowdingDistance){
//		   front.get(k).setCrowdingDistance((front.get(k).getCrowdingDistance())/5 );
//	   } else{
//		   front.get(k).setCrowdingDistance( front.get(k).getCrowdingDistance() );
//	   }   
//   }

   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;
    double minDistanceObj = 1.0e+30 ;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((objectiveMax[m] - objectiveMin[m]),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	for(int m = 0 ; m < nObjs; m++){
//    		objective_SDE[k] += front.get(k).getObjective(m);
    		objectiveValue[k] += Math.pow((front.get(k).getObjective(m) - objectiveMin[m]) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }	
    int rnd = PseudoRandom.randInt(0, size - 1);
    
    for(int k = 0; k < front.size(); k++){
//    	for(int j = 0; j < size; j++){
    		for(int m = 0; m <nObjs; m++){
    			if(front.get(k).getObjective(m) < front.get(rnd).getObjective(m))
    			{
    				objective_SDE[k] += Math.pow((front.get(rnd).getObjective(m) - front.get(k).getObjective(m))/(objectiveMax[m] - objectiveMin[m]) , 2);
    			}	    		
    		}
	   objective_SDE[k] = Math.sqrt(objective_SDE[k]);
    		objective_SDE[k] = objective_SDE[k]/nObjs;


//    	System.out.println(k+"---objective_SDE[k] ="+objective_SDE[k] );
    }
  
    avgobjectiveValue =Math.floor((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;
//System.out.println("avgobjectiveValue="+avgobjectiveValue+"---aveCrowdingDistance="+aveCrowdingDistance );
    double a = 0;
    a = 0.5;// 0.8 - 0.2 * (double)iter/maxiter;
    double b = 1 - a;
    //double sb =  PseudoRandom.randDouble(1.0, 3.0);
    //int sb =  PseudoRandom.randInt(-2, 5);
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
//    	if (objectiveValue[k] <= avgobjectiveValue ){
//    		if( front.get(k).getCrowdingDistance() > aveCrowdingDistance){	
//    			objectiveValue[k] = 1 -  objectiveValue[k] / maxdistance;
//    			front.get(k).setCrowdingDistance( objectiveValue[k] + front.get(k).getCrowdingDistance()); 
//    			}else{
//	    			objectiveValue[k] = 1 -  objectiveValue[k] / maxdistance;
//	    			front.get(k).setCrowdingDistance( objectiveValue[k] + front.get(k).getCrowdingDistance());	
//    			}
//    		
//    	}else {	    		
//    		if( front.get(k).getCrowdingDistance() > 4 * aveCrowdingDistance){
////    			System.out.println(k+"---objectiveValue[k] ="+objectiveValue[k] + "--maxdistance ="+maxdistance);
//    			objectiveValue[k] = 1 -  2 * objectiveValue[k] / maxdistance;
////    			System.out.println(k+"---objectiveValue[k] ="+objectiveValue[k] );
//    			front.get(k).setCrowdingDistance( objectiveValue[k]/5 + (0.1 * front.get(k).getCrowdingDistance())); 
//    		}else{
//    			objectiveValue[k] = 1 -  2 * objectiveValue[k] / maxdistance;
//    			front.get(k).setCrowdingDistance( objectiveValue[k]/5 + front.get(k).getCrowdingDistance()); 
//    		}
//    	}
    	//front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance() + objective_SDE[k]); 
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance;
    	    front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance() + objective_SDE[k]); 
    	    //front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance());
//    	    System.out.println(k+"----CrowdingDistance() = "+front.get(k).getCrowdingDistance());
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
            objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;	
            front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance()+ objective_SDE[k]); 
            //front.get(k).setCrowdingDistance( objectiveValue[k] +  (1+sb/5)*front.get(k).getCrowdingDistance()+ objective_SDE[k]);
            //front.get(k).setCrowdingDistance( objectiveValue[k] +  (1+sb/5)*front.get(k).getCrowdingDistance());
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
	        objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
	        front.get(k).setCrowdingDistance( objectiveValue[k] / 5 +  front.get(k).getCrowdingDistance()/ 5+ objective_SDE[k]); 
	        //front.get(k).setCrowdingDistance( objectiveValue[k]  +  front.get(k).getCrowdingDistance() + objective_SDE[k]);
	        //front.get(k).setCrowdingDistance( (1+sb/5)*objectiveValue[k]  +  front.get(k).getCrowdingDistance() + objective_SDE[k]);
	        //front.get(k).setCrowdingDistance( (1+sb/10)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
	        objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
	        front.get(k).setCrowdingDistance( objectiveValue[k] / 5 + front.get(k).getCrowdingDistance() + objective_SDE[k]); 
	        //front.get(k).setCrowdingDistance( objectiveValue[k]  + front.get(k).getCrowdingDistance() + objective_SDE[k]);
	        //front.get(k).setCrowdingDistance( (1+sb/5)*objectiveValue[k]  + (1+sb/5)*front.get(k).getCrowdingDistance() + objective_SDE[k]);
	        //front.get(k).setCrowdingDistance( (1+sb/10)*objectiveValue[k]  + (1+sb/10)*front.get(k).getCrowdingDistance());
//	        System.out.println(k+"----CrowdingDistance() = "+front.get(k).getCrowdingDistance());
    	}
    }
//    for(int k = 0; k < front.size(); k++){
//    	if( objectiveValue[k] <= avgobjectiveValue)
//    	{ 
//    	    objectiveValue[k] = (1 - objectiveValue[k]/maxdistance) ;//maxdistance
//    	}else{
//            objectiveValue[k] =  (1 - objectiveValue[k]/maxdistance ) /5; // 
//    	}
//    }
//    
//    double a = 0;
//    for( int k = 0 ; k < front.size(); k++)	{	    	
//    	    a = (objectiveValue[k]);	
//    		front.get(k).setCrowdingDistance(front.get(k).getCrowdingDistance()+(a)); 	    			
//    }
    
 //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离end----------------//
    
    for( int i = 0 ; i < size; i++){
//    System.out.println(i+"----CrowdingDistance() = "+front.get(i).getCrowdingDistance());
    }
  } // crowdingDistanceAssing5

public void crowdingDistanceAssignment44(SolutionSet solutionSet, int nObjs) {//拥挤度距离+第一套数据
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] objective_SDE  = new double [size];
    double [] objective_SDE1  = new double [size];
    double [] d1 = new double[size];
    double [] d2 = new double[size];
    double aveD1 = 0.0;
    double aveD2 = 0.0;
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
      
      double maxDistance = -1.0e+30;
	   double minDistance = 1.0e+30;  
      
      //Set de crowding distance   这里将边界点的拥挤距离设置为无穷大是不合理的，因为要考虑到比边界点可能特别的差                                                 
      for (int j = 1; j < size-1; j++) {
          distance[j] =front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i) ;                    
          distance[j] = distance[j] /(objectiveMax[i] -  objectiveMin[i]) ;  			   
          distance[j] += front.get(j).getCrowdingDistance();                
          front.get(j).setCrowdingDistance(distance[j]);   
        } // for
      distance[0] =  2*(front.get(1).getObjective(i) - front.get(0).getObjective(i));
      distance[0] = distance[0] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[0] += front.get(0).getCrowdingDistance(); 
      front.get(0).setCrowdingDistance(distance[0]);   
      distance[0] = 0;
      distance[size-1] = 2*(front.get(size-1).getObjective(i) - front.get( size-2).getObjective(i));
      distance[size-1] = distance[size-1] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[size-1] += front.get(size-1).getCrowdingDistance();  
      front.get(size-1).setCrowdingDistance(distance[size-1]); 
      distance[size-1] = 0;	      
    } // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;


   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;
    double minDistanceObj = 1.0e+30 ;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((objectiveMax[m] - objectiveMin[m]),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	for(int m = 0 ; m < nObjs; m++){
//    		objective_SDE[k] += front.get(k).getObjective(m);
    		objectiveValue[k] += Math.pow((front.get(k).getObjective(m) - objectiveMin[m]) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }	
    
    for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*objectiveMax[m];
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }
    /*double maxD1 = Double.MIN_VALUE;
    double minD1 = Double.MAX_VALUE;
    double maxD2 = Double.MIN_VALUE;
    double minD2 = Double.MAX_VALUE;
    //对d1,d2进行归一化
    for(int k=0; k<front.size(); k++){
    	if(maxD1 < d1[k]) maxD1 = d1[k];
    	if(minD1 > d1[k]) minD1 = d1[k];
    	if(maxD2 < d2[k]) maxD1 = d2[k];
    	if(minD1 > d2[k]) minD1 = d2[k];
    }
    for(int k=0; k<front.size(); k++){
    	d1[k] = (d1[k]-minD1)/(maxD1-minD1);
    	d2[k] = (d2[k]-minD2)/(maxD2-minD2);
    }*/
    
    
    for(int k=0; k<front.size(); k++){
    	aveD1 += d1[k];
    	aveD2 += d2[k];
    }
    aveD1 = Math.floor(aveD1/size);
    //aveD1 = aveD1/size;
    aveD2 = Math.floor(aveD2/size);
    //aveD2 = aveD2/size;
    
    int rnd = PseudoRandom.randInt(0, size - 1);
    for(int k = 0; k < front.size(); k++){
//    	for(int j = 0; j < size; j++){
    		for(int m = 0; m <nObjs; m++){
    			if(front.get(k).getObjective(m) < front.get(rnd).getObjective(m))
    			{
    				objective_SDE[k] += Math.pow((front.get(rnd).getObjective(m) - front.get(k).getObjective(m))/(objectiveMax[m] - objectiveMin[m]) , 2);
    			}	    		
    		}
	   objective_SDE[k] = Math.sqrt(objective_SDE[k]);
    		objective_SDE[k] = objective_SDE[k]/nObjs;
    }
  
    avgobjectiveValue =Math.floor((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
    //avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;

   // System.out.println((avgobjectiveValue)/maxdistance +" "+aveCrowdingDistance+" "+aveD1+" "+aveD2);
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance;
    	    front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance()); 
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
    		if(d2[k]>0.9*aveD2){
    			int sb4 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( objectiveValue[k] +  (3.0+((sb4-5)/3))*front.get(k).getCrowdingDistance());
    		}else{
    			int sb1 =  PseudoRandom.randInt(0, 10);
                objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;	 
                front.get(k).setCrowdingDistance( objectiveValue[k] +  (2.0+(sb1-5)/3)*front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
    		if(d1[k] < 1.4*aveD1&d2[k]>1.1*aveD2){
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (5.5)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else{
    			int sb2 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1+(sb2-5)/8)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
    		if(d1[k] < 1.4*aveD1&d2[k]>1.1*aveD2){
    			int sb3 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (1+(sb3-5)/8)*objectiveValue[k]  + (1+(sb3-5)/8)*front.get(k).getCrowdingDistance());
    		}else{
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1/5)*objectiveValue[k]  + (1/5)*front.get(k).getCrowdingDistance());
    		}
    	}
    }
  } // crowdingDistanceAssing4 


public void crowdingDistanceAssignment7(SolutionSet solutionSet, int nObjs) {//SDE
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] d1 = new double[size];
    double [] d2 = new double[size];
    double aveD1 = 0.0;
    double aveD2 = 0.0;
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
    
    DistanceShifted SDEdistance_  = new DistanceShifted();
	for (int i = 0; i < front.size(); i++) {
		Solution sss = front.get(i);
		double disSDE = Double.MAX_VALUE; 
		for(int j =0; j <  front.size() ; j ++){ 
			Solution so2 = front.get(j);
			double tempDis = SDEdistance_.distanceBetweenObjectivesShifted(sss, so2);
			if(tempDis  < disSDE){
				disSDE = tempDis;
			} 
		}   
		sss.setCrowdingDistance(disSDE);                
	} // for 
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
    } // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;


   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;
    double minDistanceObj = 1.0e+30 ;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((objectiveMax[m] - objectiveMin[m]),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	for(int m = 0 ; m < nObjs; m++){
//    		objective_SDE[k] += front.get(k).getObjective(m);
    		objectiveValue[k] += Math.pow((front.get(k).getObjective(m) - objectiveMin[m]) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }	
    
    for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*objectiveMax[m];
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }
    
    for(int k=0; k<front.size(); k++){
    	aveD1 += d1[k];
    	aveD2 += d2[k];
    }
    aveD1 = (aveD1/size);
    aveD2 = (aveD2/size);
 
    //avgobjectiveValue = ((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
    avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;
    //double sb =  PseudoRandom.randDouble(1.0, 3.0);
    //int sb1 =  PseudoRandom.randInt(-4, 5);
    //int sb2 =  PseudoRandom.randInt(-4, 5);
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance;
    	    front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance()); 
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
    		if(d2[k]>1.1*aveD2){
    			int sb4 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( objectiveValue[k] +  (3.0+((sb4-5)/6))*front.get(k).getCrowdingDistance());
    		}else{
    			int sb1 =  PseudoRandom.randInt(0, 10);
                objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;	 
                front.get(k).setCrowdingDistance( objectiveValue[k] +  (2.0+(sb1-5)/4)*front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
    		if(d1[k] < 1.1*aveD1&d2[k]>1.2*aveD2){
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (5.0)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else{
    			int sb2 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1+(sb2-5)/6)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
    		if(d1[k] < 1.2*aveD1&d2[k]>1.2*aveD2){
    			int sb3 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (1+(sb3-5)/8)*objectiveValue[k]  + (1+(sb3-5)/8)*front.get(k).getCrowdingDistance());
    		}else{
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1/5)*objectiveValue[k]  + (1/5)*front.get(k).getCrowdingDistance());
    		}
    	}
    }
  } // crowdingDistanceAssing4

public void crowdingDistanceAssignment4(SolutionSet solutionSet, int nObjs) {//SDE+第二套数据
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] d1 = new double[size];
    double [] d2 = new double[size];
    double aveD1 = 0.0;
    double aveD2 = 0.0;
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
    
    DistanceShifted SDEdistance_  = new DistanceShifted();
    double [][] SDEdistance  =  SDEdistance_.distanceMatrixShifted(front);
    int t = 1 ;
	for (int i = 0; i < SDEdistance.length; i++) {
		Arrays.sort(SDEdistance[i]);
		double tDistance =  SDEdistance[i][t] ; 
		front.get(i).setCrowdingDistance(tDistance);                
	} // for 
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
    } // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;


   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;
    double minDistanceObj = 1.0e+30 ;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((objectiveMax[m] - objectiveMin[m]),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	for(int m = 0 ; m < nObjs; m++){
//    		objective_SDE[k] += front.get(k).getObjective(m);
    		objectiveValue[k] += Math.pow((front.get(k).getObjective(m) - objectiveMin[m]) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }	
    
    for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*objectiveMax[m];
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }
    
    for(int k=0; k<front.size(); k++){
    	aveD1 += d1[k];
    	aveD2 += d2[k];
    }
    aveD1 = Math.floor(aveD1/size);
   // aveD1 = (aveD1/size);
    aveD2 = Math.floor(aveD2/size);
   // aveD2 = (aveD2/size);
 
    avgobjectiveValue =Math.floor((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
      //avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;
   
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance;
    	    front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance()); 
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
    		if(d2[k]>1.1*aveD2){
    			int sb4 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( objectiveValue[k] +  (3.0+((sb4-5)/6))*front.get(k).getCrowdingDistance());
    		}else{
    			int sb1 =  PseudoRandom.randInt(0, 10);
                objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;	 
                front.get(k).setCrowdingDistance( objectiveValue[k] +  (2.0+(sb1-5)/4)*front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
    		if(d1[k] < 1.1*aveD1&d2[k]>1.2*aveD2){
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (5)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else{
    			int sb2 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1+(sb2-5)/6)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
    		if(d1[k] < 1.2*aveD1&d2[k]>1.2*aveD2){
    			int sb3 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (1+(sb3-5)/8)*objectiveValue[k]  + (1+(sb3-5)/8)*front.get(k).getCrowdingDistance());
    		}else{
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1/5)*objectiveValue[k]  + (1/5)*front.get(k).getCrowdingDistance());
    		}
    	}
    }
  } // crowdingDistanceAssing4

public void crowdingDistanceAssignment45(SolutionSet solutionSet, int nObjs) {//拥挤度距离+第二套数据
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] objective_SDE  = new double [size];
    double [] objective_SDE1  = new double [size];
    double [] d1 = new double[size];
    double [] d2 = new double[size];
    double aveD1 = 0.0;
    double aveD2 = 0.0;
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
      
      double maxDistance = -1.0e+30;
	   double minDistance = 1.0e+30;  
      
      //Set de crowding distance   这里将边界点的拥挤距离设置为无穷大是不合理的，因为要考虑到比边界点可能特别的差                                                 
      for (int j = 1; j < size-1; j++) {
          distance[j] =front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i) ;                    
          distance[j] = distance[j] /(objectiveMax[i] -  objectiveMin[i]) ;  			   
          distance[j] += front.get(j).getCrowdingDistance();                
          front.get(j).setCrowdingDistance(distance[j]);   
        } // for
      distance[0] =  2*(front.get(1).getObjective(i) - front.get(0).getObjective(i));
      distance[0] = distance[0] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[0] += front.get(0).getCrowdingDistance(); 
      front.get(0).setCrowdingDistance(distance[0]);   
      distance[0] = 0;
      distance[size-1] = 2*(front.get(size-1).getObjective(i) - front.get( size-2).getObjective(i));
      distance[size-1] = distance[size-1] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[size-1] += front.get(size-1).getCrowdingDistance();  
      front.get(size-1).setCrowdingDistance(distance[size-1]); 
      distance[size-1] = 0;	      
    } // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;


   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;
    double minDistanceObj = 1.0e+30 ;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((objectiveMax[m] - objectiveMin[m]),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	for(int m = 0 ; m < nObjs; m++){
//    		objective_SDE[k] += front.get(k).getObjective(m);
    		objectiveValue[k] += Math.pow((front.get(k).getObjective(m) - objectiveMin[m]) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }	
    
    for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*objectiveMax[m];
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }
    /*double maxD1 = Double.MIN_VALUE;
    double minD1 = Double.MAX_VALUE;
    double maxD2 = Double.MIN_VALUE;
    double minD2 = Double.MAX_VALUE;
    //对d1,d2进行归一化
    for(int k=0; k<front.size(); k++){
    	if(maxD1 < d1[k]) maxD1 = d1[k];
    	if(minD1 > d1[k]) minD1 = d1[k];
    	if(maxD2 < d2[k]) maxD1 = d2[k];
    	if(minD1 > d2[k]) minD1 = d2[k];
    }
    for(int k=0; k<front.size(); k++){
    	d1[k] = (d1[k]-minD1)/(maxD1-minD1);
    	d2[k] = (d2[k]-minD2)/(maxD2-minD2);
    }*/
    
    
    for(int k=0; k<front.size(); k++){
    	aveD1 += d1[k];
    	aveD2 += d2[k];
    }
    //aveD1 = Math.floor(aveD1/size);
    aveD1 = aveD1/size;
    //aveD2 = Math.floor(aveD2/size);
    aveD2 = aveD2/size;
  
    avgobjectiveValue =Math.floor((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
    //avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;

   // System.out.println((avgobjectiveValue)/maxdistance +" "+aveCrowdingDistance+" "+aveD1+" "+aveD2);
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance;
    	    front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance()); 
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
    		if(d2[k]>1.1*aveD2){
    			int sb4 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( objectiveValue[k] +  (3.0+((sb4-5)/6))*front.get(k).getCrowdingDistance());
    		}else{
    			int sb1 =  PseudoRandom.randInt(0, 10);
                objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;	 
                front.get(k).setCrowdingDistance( objectiveValue[k] +  (2.0+(sb1-5)/4)*front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
    		if(d1[k] < 1.1*aveD1&d2[k]>1.2*aveD2){
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (5)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else{
    			int sb2 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1+(sb2-5)/6)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
    		if(d1[k] < 1.2*aveD1&d2[k]>1.2*aveD2){
    			int sb3 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (1+(sb3-5)/8)*objectiveValue[k]  + (1+(sb3-5)/8)*front.get(k).getCrowdingDistance());
    		}else{
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1/5)*objectiveValue[k]  + (1/5)*front.get(k).getCrowdingDistance());
    		}
    	}
    }
  } // crowdingDistanceAssing4 

public void crowdingDistanceAssignment46(SolutionSet solutionSet, int nObjs) {//拥挤度距离+第二套数据
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] objective_SDE  = new double [size];
    double [] objective_SDE1  = new double [size];
    double [] d1 = new double[size];
    double [] d2 = new double[size];
    double aveD1 = 0.0;
    double aveD2 = 0.0;
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
      
      double maxDistance = -1.0e+30;
	   double minDistance = 1.0e+30;  
      
      //Set de crowding distance   这里将边界点的拥挤距离设置为无穷大是不合理的，因为要考虑到比边界点可能特别的差                                                 
      for (int j = 1; j < size-1; j++) {
          distance[j] =front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i) ;                    
          distance[j] = distance[j] /(objectiveMax[i] -  objectiveMin[i]) ;  			   
          distance[j] += front.get(j).getCrowdingDistance();                
          front.get(j).setCrowdingDistance(distance[j]);   
        } // for
      distance[0] =  2*(front.get(1).getObjective(i) - front.get(0).getObjective(i));
      distance[0] = distance[0] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[0] += front.get(0).getCrowdingDistance(); 
      front.get(0).setCrowdingDistance(distance[0]);   
      distance[0] = 0;
      distance[size-1] = 2*(front.get(size-1).getObjective(i) - front.get( size-2).getObjective(i));
      distance[size-1] = distance[size-1] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[size-1] += front.get(size-1).getCrowdingDistance();  
      front.get(size-1).setCrowdingDistance(distance[size-1]); 
      distance[size-1] = 0;	      
    } // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;


   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;
    double minDistanceObj = 1.0e+30 ;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((objectiveMax[m] - objectiveMin[m]),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	for(int m = 0 ; m < nObjs; m++){
//    		objective_SDE[k] += front.get(k).getObjective(m);
    		objectiveValue[k] += Math.pow((front.get(k).getObjective(m) - objectiveMin[m]) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }	
    
    for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*objectiveMax[m];
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }
    /*double maxD1 = Double.MIN_VALUE;
    double minD1 = Double.MAX_VALUE;
    double maxD2 = Double.MIN_VALUE;
    double minD2 = Double.MAX_VALUE;
    //对d1,d2进行归一化
    for(int k=0; k<front.size(); k++){
    	if(maxD1 < d1[k]) maxD1 = d1[k];
    	if(minD1 > d1[k]) minD1 = d1[k];
    	if(maxD2 < d2[k]) maxD1 = d2[k];
    	if(minD1 > d2[k]) minD1 = d2[k];
    }
    for(int k=0; k<front.size(); k++){
    	d1[k] = (d1[k]-minD1)/(maxD1-minD1);
    	d2[k] = (d2[k]-minD2)/(maxD2-minD2);
    }*/
    
    
    for(int k=0; k<front.size(); k++){
    	aveD1 += d1[k];
    	aveD2 += d2[k];
    }
    //aveD1 = Math.floor(aveD1/size);
    aveD1 = aveD1/size;
    //aveD2 = Math.floor(aveD2/size);
    aveD2 = aveD2/size;
  
    avgobjectiveValue =Math.floor((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
   // avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;

   // System.out.println((avgobjectiveValue)/maxdistance +" "+aveCrowdingDistance+" "+aveD1+" "+aveD2);
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance;
    	    front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance()); 
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
    		if(d2[k]>1.1*aveD2){
    			int sb4 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( objectiveValue[k] +  (3.0+((sb4-5)/6))*front.get(k).getCrowdingDistance());
    		}else{
    			int sb1 =  PseudoRandom.randInt(0, 10);
                objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;	 
                front.get(k).setCrowdingDistance( objectiveValue[k] +  (2.0+(sb1-5)/4)*front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
    		if(d1[k] < 1.1*aveD1&d2[k]>1.2*aveD2){
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (5)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else{
    			int sb2 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1+(sb2-5)/6)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
    		if(d1[k] < 1.2*aveD1&d2[k]>1.2*aveD2){
    			int sb3 =  PseudoRandom.randInt(0, 10);
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    			front.get(k).setCrowdingDistance( (1+(sb3-5)/8)*objectiveValue[k]  + (1+(sb3-5)/8)*front.get(k).getCrowdingDistance());
    		}else{
    			objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;	
    	        front.get(k).setCrowdingDistance( (1/5)*objectiveValue[k]  + (1/5)*front.get(k).getCrowdingDistance());
    		}
    	}
    }
  } // crowdingDistanceAssing4 
public void crowdingDistanceAssignment54(SolutionSet solutionSet, int nObjs) {//拥挤度距离+第二套数据
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] objective_SDE  = new double [size];
    double [] objective_SDE1  = new double [size];
    double [] d1 = new double[size];
    double [] d2 = new double[size];
    double aveD1 = 0.0;
    double aveD2 = 0.0;
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
      
      double maxDistance = -1.0e+30;
	   double minDistance = 1.0e+30;  
      
      //Set de crowding distance   这里将边界点的拥挤距离设置为无穷大是不合理的，因为要考虑到比边界点可能特别的差                                                 
      for (int j = 1; j < size-1; j++) {
          distance[j] =front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i) ;                    
          distance[j] = distance[j] /(objectiveMax[i] -  objectiveMin[i]) ;  			   
          distance[j] += front.get(j).getCrowdingDistance();                
          front.get(j).setCrowdingDistance(distance[j]);   
        } // for
      distance[0] =  2*(front.get(1).getObjective(i) - front.get(0).getObjective(i));
      distance[0] = distance[0] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[0] += front.get(0).getCrowdingDistance(); 
      front.get(0).setCrowdingDistance(distance[0]);   
      distance[0] = 0;
      distance[size-1] = 2*(front.get(size-1).getObjective(i) - front.get( size-2).getObjective(i));
      distance[size-1] = distance[size-1] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[size-1] += front.get(size-1).getCrowdingDistance();  
      front.get(size-1).setCrowdingDistance(distance[size-1]); 
      distance[size-1] = 0;	      
    } // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;


   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;
    double minDistanceObj = 1.0e+30 ;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((objectiveMax[m] - objectiveMin[m]),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	for(int m = 0 ; m < nObjs; m++){
//    		objective_SDE[k] += front.get(k).getObjective(m);
    		objectiveValue[k] += Math.pow((front.get(k).getObjective(m) - objectiveMin[m]) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }	
    
   /* for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*objectiveMax[m];
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }*/
    for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*objectiveMax[m];
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }
    /*double maxD1 = Double.MIN_VALUE;
    double minD1 = Double.MAX_VALUE;
    double maxD2 = Double.MIN_VALUE;
    double minD2 = Double.MAX_VALUE;
    //对d1,d2进行归一化
    for(int k=0; k<front.size(); k++){
    	if(maxD1 < d1[k]) maxD1 = d1[k];
    	if(minD1 > d1[k]) minD1 = d1[k];
    	if(maxD2 < d2[k]) maxD1 = d2[k];
    	if(minD1 > d2[k]) minD1 = d2[k];
    }
    for(int k=0; k<front.size(); k++){
    	d1[k] = (d1[k]-minD1)/(maxD1-minD1);
    	d2[k] = (d2[k]-minD2)/(maxD2-minD2);
    }*/
    
    
    for(int k=0; k<front.size(); k++){
    	aveD1 += d1[k];
    	aveD2 += d2[k];
    }
    //aveD1 = Math.floor(aveD1/size);
    aveD1 = aveD1/size;
    //aveD2 = Math.floor(aveD2/size);
    aveD2 = aveD2/size;
  
    //avgobjectiveValue =Math.floor((sumobjectiveValue) /(front.size())) ;
//    avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
    avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;

   // System.out.println((avgobjectiveValue)/maxdistance +" "+aveCrowdingDistance+" "+aveD1+" "+aveD2);
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance;
    		if(d1[k]<aveD1){
    			front.get(k).setCrowdingDistance(1.1*objectiveValue[k] +  1.1*front.get(k).getCrowdingDistance());
    		}else{
    			front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance());
    		}   
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
    		objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;
    		if(d1[k]<aveD1){
    			int sb1 =  PseudoRandom.randInt(0, 10);  			
    			front.get(k).setCrowdingDistance( objectiveValue[k] + (2+(sb1-5)/4)*front.get(k).getCrowdingDistance());
    		}else{	 
                front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
    		objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    		if(d1[k] < 2.0*aveD1&d2[k]>1.0*aveD2){
    			front.get(k).setCrowdingDistance( (10.0)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else{
    			double sb2 =  PseudoRandom.randDouble(0.2, 1.0);	
    	        front.get(k).setCrowdingDistance( (sb2)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
    		objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    		if(d1[k] < 2.0*aveD1&d2[k]>1.0*aveD2){
    			int sb3 =  PseudoRandom.randInt(0, 10);
    			front.get(k).setCrowdingDistance( (4+(sb3-5)/4)*objectiveValue[k]  + (4+(sb3-5)/4)*front.get(k).getCrowdingDistance());
    		}else{	
    	        front.get(k).setCrowdingDistance( (1/5)*objectiveValue[k]  + (1/5)*front.get(k).getCrowdingDistance());
    		}
    	}
    }
  } // crowdingDistanceAssing4

public void crowdingDistanceAssignment55(SolutionSet solutionSet, int nObjs) {//拥挤度距离+第二套数据
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double [] objective_SDE  = new double [size];
    double [] objective_SDE1  = new double [size];
    double [] d1 = new double[size];
    double [] d2 = new double[size];
    double aveD1 = 0.0;
    double aveD2 = 0.0;
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
    }
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((1 - 0),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	double ss = 0.0;
    	for(int m = 0 ; m < nObjs; m++){
    		ss = (front.get(k).getObjective(m) - objectiveMin[m])/(objectiveMax[m] - objectiveMin[m]);
            front.get(k).setNormalizedObjective(m, ss);
    		objectiveValue[k] += Math.pow((front.get(k).getNormalizedObjective(m)) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }
  //avgobjectiveValue =Math.floor((sumobjectiveValue) /(front.size())) ;
  //avgobjectiveValue =Math.ceil((sumobjectiveValue) /(front.size())) ;
  avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;
     
  /*for (int i = 0; i<nObjs; i++) { 
      //Set de crowding distance   这里将边界点的拥挤距离设置为无穷大是不合理的，因为要考虑到比边界点可能特别的差     
	  front.sort(new NormalizedObjectiveComparator(i));
	  objectiveMin[i] = front.get(0).getNormalizedObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getNormalizedObjective(i); //最大的值
      for (int j = 1; j < size-1; j++) {
          distance[j] =front.get(j+1).getNormalizedObjective(i) - front.get(j-1).getNormalizedObjective(i) ;                    
          distance[j] = distance[j] /(objectiveMax[i] -  objectiveMin[i]) ;  			   
          distance[j] += front.get(j).getCrowdingDistance();                
          front.get(j).setCrowdingDistance(distance[j]);   
        } // for
      distance[0] =  2*(front.get(1).getNormalizedObjective(i) - front.get(0).getNormalizedObjective(i));
      distance[0] = distance[0] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[0] += front.get(0).getCrowdingDistance(); 
      front.get(0).setCrowdingDistance(distance[0]);   
      distance[0] = 0;
      distance[size-1] = 2*(front.get(size-1).getNormalizedObjective(i) - front.get( size-2).getNormalizedObjective(i));
      distance[size-1] = distance[size-1] / (objectiveMax[i] -  objectiveMin[i]); 
      distance[size-1] += front.get(size-1).getCrowdingDistance();  
      front.get(size-1).setCrowdingDistance(distance[size-1]); 
      distance[size-1] = 0;	      
    } // for*/
  
  DistanceShifted SDEdistance_  = new DistanceShifted();
  double [][] SDEdistance  =  SDEdistance_.normalDistanceMatrixShifted(front);
  int t = 1 ;
	for (int i = 0; i < SDEdistance.length; i++) {
		Arrays.sort(SDEdistance[i]);
		double tDistance =  SDEdistance[i][t] ; 
		front.get(i).setCrowdingDistance(tDistance);                
	} // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;


   //------------------------对边界的点不再设置为无穷大，而是设置为平均的拥挤距离start----------------//
    
   /* for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getObjective(m) - objectiveMin[m])*(objectiveMax[m]);
    		nl += (objectiveMax[m]*objectiveMax[m]);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl))*((front.get(k).getObjective(n)-objectiveMin[n]) - d1[k]*(objectiveMax[n]/nl));
    	}
    	d2[k] = Math.sqrt(d2[k]);
    }*/
    for(int k=0; k<front.size(); k++){
    	double nl = 0.0 ;
    	d1[k] = 0.0;
    	d2[k] = 0.0;
    	for(int m=0; m<nObjs;m++){
    		d1[k] += (front.get(k).getNormalizedObjective(m) - 0)*1;
    		nl += (1*1);
    	}
    	nl = Math.sqrt(nl);
    	d1[k] = Math.abs(d1[k]/nl);
    	for(int n=0; n<nObjs;n++){
    		d2[k] +=((front.get(k).getNormalizedObjective(n)-0) - d1[k]*(1/nl))*((front.get(k).getNormalizedObjective(n)-0) - d1[k]*(1/nl));
    		//d2[k] +=((front.get(k).getNormalizedObjective(n)-0) )*((front.get(k).getNormalizedObjective(n)-0));
    	}
    	//d2[k] = d2[k] - d1[k]*d1[k];
    	d2[k] = Math.sqrt(d2[k]);
    }
    
    for(int k=0; k<front.size(); k++){
    	aveD1 += d1[k];
    	aveD2 += d2[k];
    }
    //aveD1 = Math.floor(aveD1/size);
    aveD1 = aveD1/size;
    //aveD2 = Math.floor(aveD2/size);
    aveD2 = aveD2/size;
  
 

   // System.out.println((avgobjectiveValue)/maxdistance +" "+aveCrowdingDistance+" "+aveD1+" "+aveD2);
    for(int k = 0; k < front.size(); k++){    
    	/**
    	 * 这里主要针对收敛性的两种情况的考虑，收敛性好时的多样性的考虑，以及收敛性较差时的多样性的考虑
    	 * 当收敛性好，则多样性值大、大于平均的多样性值，表示多样性较好，否则表示多样性不好
    	 * 当收敛不好的情况下，多样性的值较大时，则表示是属于独立较远的点，当多样性小于平均值，则表示单纯的是、多样性不够好
    	 * */
    	if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() >  aveCrowdingDistance){ 
    		objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance ;
    		if(d1[k]>aveD1){
    			front.get(k).setCrowdingDistance(0.9*objectiveValue[k] +  0.9*front.get(k).getCrowdingDistance());
    		}else{
    			front.get(k).setCrowdingDistance( objectiveValue[k] +  front.get(k).getCrowdingDistance());
    		}   
    	}else if( objectiveValue[k] <= avgobjectiveValue && front.get(k).getCrowdingDistance() <=  aveCrowdingDistance){
    		objectiveValue[k] =  1 -  1.0 * objectiveValue[k] / maxdistance;
    		if(d1[k]>aveD1){  			
    			front.get(k).setCrowdingDistance( objectiveValue[k] + (0.6)*front.get(k).getCrowdingDistance());
    		}else{
    			int sb1 =  PseudoRandom.randInt(0, 5);
                front.get(k).setCrowdingDistance( objectiveValue[k] +  (1+(sb1-3)/8)*front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() > 1 * aveCrowdingDistance){
    		objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    		if(d2[k]<aveD2){
    			front.get(k).setCrowdingDistance( (0.2)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else if(d1[k]>1.0*avgobjectiveValue){
    			front.get(k).setCrowdingDistance( (0.2)*objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}else{
    			//double sb2 =  PseudoRandom.randDouble(0.2, 1.0);	
    	        front.get(k).setCrowdingDistance( objectiveValue[k]  +  front.get(k).getCrowdingDistance());
    		}
    	}else if(objectiveValue[k] > avgobjectiveValue && front.get(k).getCrowdingDistance() <= 1 * aveCrowdingDistance){
    		objectiveValue[k] =  1 - 1.0 * objectiveValue[k] / maxdistance;
    		if(d2[k]<aveD2){
    			front.get(k).setCrowdingDistance( (0.2)*objectiveValue[k]  + (0.2)*front.get(k).getCrowdingDistance());
    		}else if(d1[k]>1.0*avgobjectiveValue){
    			front.get(k).setCrowdingDistance( (0.2)*objectiveValue[k]  +  0.2*front.get(k).getCrowdingDistance());
    		}else{
    			int sb3 =  PseudoRandom.randInt(0, 5);
    	        front.get(k).setCrowdingDistance( (1+(sb3-3)/8)*objectiveValue[k]  + (1+(sb3-3)/8)*front.get(k).getCrowdingDistance());
    		}
    	}
    }
  } // crowdingDistanceAssing4
public void crowdingDistanceAssignment77(SolutionSet solutionSet, int nObjs) {//拥挤度距离+第二套数据
    int size = solutionSet.size();        
                
    if (size == 0)
      return;    
    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if
        
    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(1);
      solutionSet.get(1).setCrowdingDistance(1);
      return;
    } // if       
        
    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }
        
    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);        
        
    double objetiveMaxn;
    double objetiveMinn;
    double [] objectiveMin = new double[nObjs];
    double [] objectiveMax = new double[nObjs];
    double []distance = new double [size];
    double [] objectiveValue = new double [size];
    double maxCrowdistance  = -1.0e+30 ;
    double minCrowdistance  =  1.0e+30 ;
                
    for (int i = 0; i<nObjs; i++) {          
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      objectiveMin[i] = front.get(0).getObjective(i);   //最小的值 
      objectiveMax[i] = front.get(front.size()-1).getObjective(i); //最大的值
    }
    double maxdistance = 0.0;
    double sumobjectiveValue =0.0;
    double avgobjectiveValue =0.0;

    for(int m = 0 ; m < nObjs; m++){
    	maxdistance += Math.pow((1 - 0),2);
    }
        maxdistance = Math.sqrt(maxdistance);
   
    for(int k = 0; k < front.size(); k++){
    	double ss = 0.0;
    	for(int m = 0 ; m < nObjs; m++){
    		ss = (front.get(k).getObjective(m) - objectiveMin[m])/(objectiveMax[m] - objectiveMin[m]);
            front.get(k).setNormalizedObjective(m, ss);
    		objectiveValue[k] += Math.pow((front.get(k).getNormalizedObjective(m)) , 2); //objectiveMin[m]	    		
    	}
    	   objectiveValue[k] = Math.sqrt(objectiveValue[k]);
    	   sumobjectiveValue += objectiveValue[k];
    }
  avgobjectiveValue =((sumobjectiveValue) /(front.size())) ;
  
  DistanceShifted SDEdistance_  = new DistanceShifted();
  double [][] SDEdistance  =  SDEdistance_.normalDistanceMatrixShifted(front);
  int t = 1 ;
	for (int i = 0; i < SDEdistance.length; i++) {
		Arrays.sort(SDEdistance[i]);
		double tDistance =  SDEdistance[i][t] ; 
		front.get(i).setCrowdingDistance(tDistance);                
	} // for 
    
   for(int n = 0; n<size; n++){
	   if (front.get(n).getCrowdingDistance() > maxCrowdistance)
		   maxCrowdistance = front.get(n).getCrowdingDistance(); 
	   if(front.get(n).getCrowdingDistance() < minCrowdistance)
		   minCrowdistance = front.get(n).getCrowdingDistance(); 
   } 

//对拥挤距离进行归一化处理
   double sumCrowdingDistance = 0.0;
   double aveCrowdingDistance = 0.0;

   for(int n = 0 ; n < size; n ++){
	   front.get(n).setCrowdingDistance((front.get(n).getCrowdingDistance()-minCrowdistance )/(maxCrowdistance - minCrowdistance));
	   sumCrowdingDistance +=  front.get(n).getCrowdingDistance();  
   }
   
   aveCrowdingDistance = sumCrowdingDistance/size;

    for(int k = 0; k < front.size(); k++){    
    	objectiveValue[k] = 1 -  1.0 * objectiveValue[k] / maxdistance ;
    	front.get(k).setCrowdingDistance(0.1*objectiveValue[k] +  0.9*front.get(k).getCrowdingDistance());			
    }
  } // crowdingDistanceAssing4

} // Distance

