//  Distance.java
//
//  Author:
//       lbd

package jmetal.util;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;
import jmetal.util.comparators.ObjectiveComparator;

/**
 * This class implements some utilities for calculating distances
 */
public class DistanceShifted {      
    
  /** 
  * Constructor.
  */
  public DistanceShifted() {
    //do nothing.
  } // Distance
        
    
  /** 
  * Returns a matrix with distances between solutions in a 
  * <code>SolutionSet</code>.
  * @param solutionSet The <code>SolutionSet</code>.
  * @return a matrix with distances.
  */
  public double [][] distanceMatrixShifted(SolutionSet solutionSet) {
    Solution solutionI, solutionJ;

    //The matrix of distances
    double [][] distance = new double [solutionSet.size()][solutionSet.size()];        
    //-> Calculate the distances
    for (int i = 0; i < solutionSet.size(); i++){
      distance[i][i] = 0.0;
      solutionI = solutionSet.get(i);
      for (int j = i + 1; j < solutionSet.size(); j++){
        solutionJ = solutionSet.get(j);
        distance[i][j] = this.distanceBetweenObjectivesShifted(solutionI,solutionJ);                
        distance[j][i] = this.distanceBetweenObjectivesShifted(solutionJ,solutionI);             
      } // for
    } // for 
    //->Return the matrix of distances
    return distance;
  } // distanceMatrix
  
 public double [][] normalDistanceMatrixShifted(SolutionSet solutionSet) {
        Solution solutionI, solutionJ;

        //The matrix of distances
        double [][] distance = new double [solutionSet.size()][solutionSet.size()];        
        //-> Calculate the distances
        for (int i = 0; i < solutionSet.size(); i++){
          distance[i][i] = 0.0;
          solutionI = solutionSet.get(i);
          for (int j = i + 1; j < solutionSet.size(); j++){
            solutionJ = solutionSet.get(j);
            distance[i][j] = this.distanceBetweenNormalizedObjectivesShifted(solutionI,solutionJ);                
            distance[j][i] = this.distanceBetweenNormalizedObjectivesShifted(solutionJ,solutionI);             
          } // for
        } // for        
    
    //->Return the matrix of distances
    return distance;
  } // distanceMatrix
 
 public double [][] translateDistanceMatrixShifted(SolutionSet solutionSet) {
     Solution solutionI, solutionJ;

     //The matrix of distances
     double [][] distance = new double [solutionSet.size()][solutionSet.size()];        
     //-> Calculate the distances
     for (int i = 0; i < solutionSet.size(); i++){
       distance[i][i] = 0.0;
       solutionI = solutionSet.get(i);
       for (int j = i + 1; j < solutionSet.size(); j++){
         solutionJ = solutionSet.get(j);
         distance[i][j] = this.distanceBetweenObjectivesShiftedTranslated(solutionI,solutionJ);                
         distance[j][i] = this.distanceBetweenObjectivesShiftedTranslated(solutionJ,solutionI);             
       } // for
     } // for        
 
 //->Return the matrix of distances
 return distance;
} // distanceMatrix
  
  public double [][] distanceMatrixShifted1(SolutionSet solutionSet) {
	    Solution solutionI, solutionJ;

	    //The matrix of distances
	    double [][] distance = new double [solutionSet.size()][solutionSet.size()];        
	    //-> Calculate the distances
	    for (int i = 0; i < solutionSet.size(); i++){
	      distance[i][i] = 0.0;
	      solutionI = solutionSet.get(i);
	      for (int j = 0; j < solutionSet.size(); j++){
	        solutionJ = solutionSet.get(j);
	        if(i == j){
	        	distance[i][j] = Double.POSITIVE_INFINITY;
	        }else{
	        	 distance[i][j] = this.distanceBetweenObjectivesShifted(solutionI,solutionJ); 
	        }
	        //distance[i][j] = this.distanceBetweenObjectivesShifted(solutionI,solutionJ);                
	        //distance[j][i] = this.distanceBetweenObjectivesShifted(solutionJ,solutionI);             
	      } // for
	    } // for        
	    
	    //->Return the matrix of distances
	    return distance;
	  } // distanceMatrix
  
  /** 
   * Returns a matrix with distances between solutions in a 
   * <code>SolutionSet</code>.
   * @param solutionSet The <code>SolutionSet</code>.
   * @return a matrix with distances.
   */
   public double [][] distanceMatrixShiftedalpha(SolutionSet solutionSet,
		   double alpha) {
     Solution solutionI, solutionJ;

     //The matrix of distances
     double [][] distance = new double [solutionSet.size()][solutionSet.size()];        
     //-> Calculate the distances
     for (int i = 0; i < solutionSet.size(); i++){
       distance[i][i] = 0.0;
       solutionI = solutionSet.get(i);
       for (int j = i + 1; j < solutionSet.size(); j++){
         solutionJ = solutionSet.get(j);
         distance[i][j] = this.distanceBetweenObjectivesShiftedalpha(solutionI,solutionJ, alpha);                
         distance[j][i] = this.distanceBetweenObjectivesShiftedalpha(solutionJ,solutionI, alpha);             
       } // for
     } // for        
     
     //->Return the matrix of distances
     return distance;
   } // distanceMatrix
   public double kNNdistanceToSolutionSetInObjectiveSpaceShifted(int k,
		  Solution    solution, 
          SolutionSet solutionSet) throws JMException{
		//At start point the distance is the max
		double[] minkDistances = new double[k];
		for(int i = 0; i< k; i ++){
			minkDistances[i] = Double.MAX_VALUE;    
		}

		// found the min distance respect to population
		for (int i = 0; i < solutionSet.size();i++){            
			double aux = this.distanceBetweenObjectivesShifted(solution,solutionSet.get(i));
			int location = -1;
			for(int j = 0; j < minkDistances.length; j ++){ // = 0; k <)
				if (aux < minkDistances[j]){
					location = j;
					break;
				}
			}
			if(location == -1){
				// bigger than all the k min values, ignored
			}
			else if(location >=0 && location <minkDistances.length){
			//	System.out.println("location == " + location);
			//	System.out.println("length ==" + minkDistances.length);
				
				for(int tt = minkDistances.length -2; tt >= location&& tt >=0 && tt+1< minkDistances.length; tt--) {
				//	System.out.println("tt == " + tt);
					minkDistances[tt+1] = minkDistances[tt];
				}
				minkDistances[location] = aux;
			}
 		} // for
		double sumd = 0;
		for(int i = 0; i < minkDistances.length; i ++) {
			sumd += minkDistances[i];
		}
		double dist = sumd/(double)minkDistances.length;
		//->Return the best distance
		return dist;
} // distanceToSolutionSetinObjectiveSpace
 
 /** Returns the minimum distance from a <code>Solution</code> to a 
  * <code>SolutionSet according to the objective values</code>.
  * @param solution The <code>Solution</code>.
  * @param solutionSet The <code>SolutionSet</code>.
  * @return The minimum distance between solution and the set.
 * @throws JMException 
  */  
  public double distanceToSolutionSetInObjectiveSpaceShifted(Solution    solution, 
		                                  SolutionSet solutionSet) throws JMException{
    //At start point the distance is the max
    double distance = Double.MAX_VALUE;    
        
    // found the min distance respect to population
    for (int i = 0; i < solutionSet.size();i++){            
      double aux = this.distanceBetweenObjectivesShifted(solution,solutionSet.get(i));
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
   public double distanceToSolutionSetInSolutionSpaceShifted(Solution    solution, 
 		                                  SolutionSet solutionSet) throws JMException{
     //At start point the distance is the max
     double distance = Double.MAX_VALUE;    
         
     // found the min distance respect to population
     for (int i = 0; i < solutionSet.size();i++){            
       double aux = this.distanceBetweenSolutionsShifted(solution,solutionSet.get(i));
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
  public double distanceBetweenSolutionsShifted(Solution solutionI, Solution solutionJ) 
  throws JMException{                
    double distance = 0.0;
    if ((solutionI.getDecisionVariables() != null) && 
    		(solutionJ.getDecisionVariables() != null)) {
      Variable[] decisionVariableI = solutionI.getDecisionVariables();
      Variable[] decisionVariableJ = solutionJ.getDecisionVariables();    
    
      double diff = 0;    //Auxiliar var
      //-> Calculate the Euclidean distance
      for (int i = 0; i < decisionVariableI.length; i++){
    	  if(decisionVariableI[i].getValue() < decisionVariableJ[i].getValue()){
    		  //only consider the distance where solution I performs better than J
    		  // the other worse dimensions are just ignored
    		  // in other words, dis == 0 is used to measure the distance on these worse dimensions
    		  diff = decisionVariableI[i].getValue() -  decisionVariableJ[i].getValue();
    	  }
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
  public double distanceBetweenObjectivesShifted(Solution solutionI, Solution solutionJ){                
    double diff = 0;    //Auxiliar var
    double distance = 0.0;
    //-> Calculate the euclidean distance
    for (int nObj = 0; nObj < solutionI.numberOfObjectives();nObj++){
    	diff = 0;
    	if(solutionI.getObjective(nObj) < solutionJ.getObjective(nObj)){
    		diff = Math.abs(solutionI.getObjective(nObj) - solutionJ.getObjective(nObj));
    	}
    	distance += diff*diff;           
    } // for   
        
    //Return the euclidean distance
    return Math.sqrt(distance);
  } // distanceBetweenObjectives.LpNorm
  
  public double distanceBetweenNormalizedObjectivesShifted(Solution solutionI, Solution solutionJ){                
	    double diff = 0;    //Auxiliar var
	    double distance = 0.0;
	    //-> Calculate the euclidean distance
	    for (int nObj = 0; nObj < solutionI.numberOfObjectives();nObj++){
	    	diff = 0;
	    	if(solutionI.getNormalizedObjective(nObj) < solutionJ.getNormalizedObjective(nObj)){
	    		diff = Math.abs(solutionI.getNormalizedObjective(nObj) - solutionJ.getNormalizedObjective(nObj));
	    	}
	    	distance += diff*diff;           
	    } // for   
	        
	    //Return the euclidean distance
	    return Math.sqrt(distance);
	  } // distanceBetweenObjectives.LpNorm
  
  public double distanceBetweenObjectivesShiftedTranslated(Solution solutionI, Solution solutionJ){                
	    double diff = 0;    //Auxiliar var
	    double distance = 0.0;
	    //-> Calculate the euclidean distance
	    for (int nObj = 0; nObj < solutionI.numberOfObjectives();nObj++){
	    	diff = 0;
	    	if(solutionI.getIthTranslatedObjective(nObj) < solutionJ.getIthTranslatedObjective(nObj)){
	    		diff = solutionI.getIthTranslatedObjective(nObj) - solutionJ.getIthTranslatedObjective(nObj);
	    	}
	    	distance += diff*diff;                
	    } // for   
	        
	    //Return the euclidean distance
	    return Math.sqrt(distance);
	  } // distanceBetweenObjectives.LpNorm
	  
  public double distanceBetweenObjectivesShiftedReverse(Solution solutionI, Solution solutionJ){                
	    double diff = 0;    //Auxiliar var
	    double distance = 0.0;
	    //-> Calculate the euclidean distance
	    for (int nObj = 0; nObj < solutionI.numberOfObjectives();nObj++){
	    	if(solutionI.getObjective(nObj) > solutionJ.getObjective(nObj)){
	    		diff = solutionI.getObjective(nObj) - solutionJ.getObjective(nObj);
	    	}
	    	distance += Math.pow(diff,2.0);           
	    } // for   
	        
	    //Return the euclidean distance
	    return Math.sqrt(distance);
	  } // distanceBetweenObjectives.LpNorm
  public double distanceBetweenObjectivesShiftedLpNorm(Solution solutionI, Solution solutionJ,double lp){                
	    double diff = 0;    //Auxiliar var
	    double distance = 0.0;
	    //-> Calculate the euclidean distance
	    for (int nObj = 0; nObj < solutionI.numberOfObjectives();nObj++){
	    	if(solutionI.getObjective(nObj) < solutionJ.getObjective(nObj)){
	    		diff = solutionJ.getObjective(nObj) - solutionI.getObjective(nObj);
	    	}
	    	distance += Math.pow(diff,lp);           
	    } // for   
	    distance = Math.pow(distance,1.0/lp);  
	    //Return the euclidean distance
	    return distance;
	  } // distanceBetweenObjectives.       
  
  public double distanceBetweenObjectivesShiftedalpha(Solution solutionI, Solution solutionJ, double alpha){                
	    double diff = 0;    //Auxiliar var
	    double distance = 0.0;
	    //-> Calculate the euclidean distance
	    for (int nObj = 0; nObj < solutionI.numberOfObjectives();nObj++){
	    	if(solutionI.getObjective(nObj) < solutionJ.getObjective(nObj)){
	    		diff =  Math.abs(solutionI.getObjective(nObj) - solutionJ.getObjective(nObj));
	    	}
	    	else {
	    		diff = alpha * Math.abs(solutionI.getObjective(nObj) - solutionJ.getObjective(nObj));
	    	}
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


public double[][] distanceMatrixShiftedLpNorm(SolutionSet solutionSet,
		double lpnorm) {
	 Solution solutionI, solutionJ;

	    //The matrix of distances
	    double [][] distance = new double [solutionSet.size()][solutionSet.size()];        
	    //-> Calculate the distances
	    for (int i = 0; i < solutionSet.size(); i++){
	      distance[i][i] = 0.0;
	      solutionI = solutionSet.get(i);
	      for (int j = i + 1; j < solutionSet.size(); j++){
	        solutionJ = solutionSet.get(j);
	        distance[i][j] = this.distanceBetweenObjectivesShiftedLpNorm(solutionI,solutionJ, lpnorm);                
	        distance[j][i] = this.distanceBetweenObjectivesShiftedLpNorm(solutionJ,solutionI, lpnorm);             
	      } // for
	    } // for        
	    
	    //->Return the matrix of distances
	    return distance;
}
} // Distance

