//  Front.java
//
//  Authors:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2013 Antonio J. Nebro, Juan J. Durillo
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
//

//  CREDIT
//  This class is based on the code of the WFG group (http://www.wfg.csse.uwa.edu.au/hypervolume/)
//  Copyright (C) 2010 Lyndon While, Lucas Bradstreet.

package jmetal.qualityIndicator.fastHypervolume.wfg;

import jmetal.core.SolutionSet;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: antelverde
 * Date: 25/07/13
 * Time: 16:26
 * To change this template use File | Settings | File Templates.
 */
public class Front {
  private int numberOfPoints_ ;
  private int dimension_ ;
  public Point[] points_ ;
  private boolean maximizing_ ;
  public int nPoints_ ;

  private Comparator pointComparator ;

  public Front() {
    maximizing_ = true ;
    pointComparator = new PointComparator(maximizing_) ;

  }

  public Front (int numberOfPoints, int dimension, SolutionSet solutionSet) {
    maximizing_ = true ;
    pointComparator = new PointComparator(maximizing_) ;
    numberOfPoints_ = numberOfPoints ;
    dimension_ = dimension  ;
    nPoints_ = numberOfPoints_ ;

    points_ = new Point[numberOfPoints_] ;
    for (int i = 0 ; i < numberOfPoints_; i++) {
      double [] p = new double[dimension] ;
      for (int j = 0 ; j < dimension; j++) {
         p[j] = solutionSet.get(i).getObjective(j) ;
      }
      points_[i] = new Point(p) ;
    }
  }

  public Front (int numberOfPoints, int dimension) {
    maximizing_ = true ;
    pointComparator = new PointComparator(maximizing_) ;
    numberOfPoints_ = numberOfPoints ;
    dimension_ = dimension  ;
    nPoints_ = numberOfPoints_ ;

    points_ = new Point[numberOfPoints_] ;
    for (int i = 0 ; i < numberOfPoints_; i++) {
      double [] p = new double[dimension] ;
      for (int j = 0 ; j < dimension; j++) {
        p[j] = 0.0 ;
      }
      points_[i] = new Point(p) ;
    }
  }

  public Front(int numberOfPoints, int dimension, List<double[]> listOfPoints) {
    maximizing_ = true ;
    pointComparator = new PointComparator(maximizing_) ;
    numberOfPoints_ = numberOfPoints ;
    dimension_ = dimension  ;

    points_ = new Point[numberOfPoints_] ;
    for (int i = 0 ; i < numberOfPoints_; i++) {
      points_[i] = new Point(listOfPoints.get(i)) ;
    }
  }


  public void readFront(String fileName, Point points) throws IOException {
    FileInputStream fis   = new FileInputStream(fileName)  ;
    InputStreamReader isr = new InputStreamReader(fis)    ;
    BufferedReader br      = new BufferedReader(isr)      ;

    List<double []> list = new ArrayList<double []>();
    int numberOfObjectives = 0;
    String aux = br.readLine();
    while (aux!= null) {
      StringTokenizer st = new StringTokenizer(aux);
      int i = 0;
      numberOfObjectives = st.countTokens();
      if (numberOfObjectives != points.getNumberOfObjectives()){
    	  System.out.println("objectives number not met");
    	  return;
      }
      double [] vector = new double[st.countTokens()];
      boolean remove_flag = false;
      while (st.hasMoreTokens()) {
        double value = new Double(st.nextToken());
        if (value > points.objectives_[i]){
        	remove_flag=true;break;
        }
        vector[i] = value;
        i++;
      }
      if(remove_flag == false)
    	  list.add(vector);
      aux = br.readLine();
    }
    numberOfPoints_ = list.size() ;
    dimension_ = numberOfObjectives ;
    points_ = new Point[numberOfPoints_] ;
    nPoints_ = numberOfPoints_ ;
    for (int i = 0 ; i < numberOfPoints_; i++) {
      points_[i] = new Point(list.get(i)) ;
    }
  }
  
  
  public void loadFront(SolutionSet solutionSet, int notLoadingIndex) {          
      
      if (notLoadingIndex >= 0 && notLoadingIndex < solutionSet.size()) {
          numberOfPoints_ = solutionSet.size()-1;
      } else {
          numberOfPoints_ = solutionSet.size();
      }
          
      nPoints_ = numberOfPoints_;
      dimension_ = solutionSet.get(0).numberOfObjectives();
      
      points_ = new Point[numberOfPoints_];
      
      int index = 0;
      for (int i = 0; i < solutionSet.size(); i++) {
        if (i != notLoadingIndex) {          
          double [] vector = new double[dimension_];
          for (int j = 0; j < dimension_; j++) {
              vector[j] = solutionSet.get(i).getObjective(j);
          }                        
          points_[index++] = new Point(vector);            
        }
      }
  }
public void loadFront1(double solutionSet [][]) {          
      
      numberOfPoints_ = solutionSet.length;
          
      nPoints_ = numberOfPoints_;
      dimension_ = solutionSet[0].length;
      
      points_ = new Point[numberOfPoints_];
      
      int index = 0;
      for (int i = 0; i < nPoints_; i++) {
          double [] vector = new double[dimension_];
          vector = solutionSet[i];
          points_[index++] = new Point(vector);   
      }
  }
  

  public void printFront() {
    System.out.println("Objectives:       " + dimension_) ;
    System.out.println("Number of points: " + numberOfPoints_) ;

    for (Point point : points_) {
      System.out.println(point) ;
    }
  }

  public int getNumberOfObjectives() {
    return dimension_ ;
  }

  public int getNumberOfPoints() {
    return numberOfPoints_;
  }

  public Point getPoint(int index) {
    return points_[index] ;
  }

  public Point[] getPoints() {
    return points_ ;
  }

  public void setToMazimize() {
    maximizing_ = true ;
    pointComparator = new PointComparator(maximizing_) ;
  }

  public void setToMinimize() {
    maximizing_ = false ;
    pointComparator = new PointComparator(maximizing_) ;
  }

  public void sort() {
    Arrays.sort(points_, pointComparator);
  }

  public Point getReferencePoint(int fun) {
    Point referencePoint = new Point(dimension_) ;

	for (int j=0;j<dimension_;j++){
		if(fun == 6){
			referencePoint.setObjective(j,0.5);
		}else if(fun>6&fun<=8){
			referencePoint.setObjective(j,1.0);
		}else if(fun>8&fun<=11){
			if(j!=dimension_-1){
				referencePoint.setObjective(dimension_-1-j,Math.pow(Math.sqrt(2)/2, j));
			}else{
				referencePoint.setObjective(0,referencePoint.getObjective(1));
			}
		}else if(fun>12&fun<=21){
			referencePoint.setObjective(j,2.0*(j+1));
		}
	}

    return referencePoint ;
  }
  public Point getMinReferencePoint() {
	    Point referencePoint = new Point(dimension_) ;

	    double [] minObjectives = new double[dimension_] ;
	    for (int i = 0; i < dimension_; i++)
	      minObjectives[i] = Double.MAX_VALUE ;

	    for (int i = 0; i < points_.length; i++)
	      for (int j = 0 ; j < dimension_; j++)
	        if (minObjectives[j] > points_[i].objectives_[j])
	          minObjectives[j] = points_[i].objectives_[j] ;

	    //for (int i = 0; i < solution.getNumberOfObjectives(); i++) {
	    //  if (maxObjectives[i] < referencePoint_.objectives_[i])
	    //    referencePoint_.objectives_[i] = maxObjectives[i] ;
	    //
	    // }
	    for (int i = 0; i < dimension_; i++)
	      referencePoint.objectives_[i] = minObjectives[i] ;

	    return referencePoint ;
	  }
}
