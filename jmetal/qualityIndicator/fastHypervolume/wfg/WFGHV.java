//  WFGHV.java
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

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Problem;
import jmetal.util.Configuration;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Created with IntelliJ IDEA.
 * User: Antonio J. Nebro
 * Date: 25/07/13
 * Time: 17:50
 * To change this template use File | Settings | File Templates.
 */
public class WFGHV {
  Front [] fs_ ;
  Point referencePoint_ ;
  boolean maximizing_  ;
  int currentDeep_  ;
  int currentDimension_ ;
  int maxNumberOfPoints_ ;
  int maxNumberOfObjectives_ ;
  final int OPT = 2 ;
  Comparator pointComparator_;
  public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
	      BufferedWriter bw      = new BufferedWriter(osw)        ;               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");
	        bw.newLine();        
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD

  public WFGHV(int dimension, int maxNumberOfPoints) {
    referencePoint_ = null ;
    maximizing_ = false ;
    currentDeep_ = 0 ;
    currentDimension_ = dimension ;
    maxNumberOfPoints_ = maxNumberOfPoints ;
    maxNumberOfObjectives_ = dimension ;
    pointComparator_ = new PointComparator(true) ;

    int maxd = maxNumberOfPoints_ - (OPT /2 + 1) ;
    fs_ = new Front[maxd] ;
    for (int i = 0; i < maxd; i++) {
      fs_[i] = new Front(maxNumberOfPoints, dimension) ;
    }
  }

  public WFGHV(int dimension, int maxNumberOfPoints, Solution referencePoint) {
    referencePoint_ = new Point(referencePoint) ;
    maximizing_ = false ;
    currentDeep_ = 0 ;
    currentDimension_ = dimension ;
    maxNumberOfPoints_ = maxNumberOfPoints ;
    maxNumberOfObjectives_ = dimension ;
    pointComparator_ = new PointComparator(true) ;

    int maxd = maxNumberOfPoints_ - (OPT /2 + 1) ;
    fs_ = new Front[maxd] ;
    for (int i = 0; i < maxd; i++) {
      fs_[i] = new Front(maxNumberOfPoints, dimension) ;
    }
  }

  public WFGHV(int dimension, int maxNumberOfPoints, Point referencePoint) {
    referencePoint_ = referencePoint ;
    maximizing_ = false ;
    currentDeep_ = 0 ;
    currentDimension_ = dimension ;
    maxNumberOfPoints_ = maxNumberOfPoints ;
    maxNumberOfObjectives_ = dimension ;
    pointComparator_ = new PointComparator(true) ;

    int maxd = maxNumberOfPoints_ - (OPT /2 + 1) ;
    fs_ = new Front[maxd] ;
    for (int i = 0; i < maxd; i++) {
      fs_[i] = new Front(maxNumberOfPoints, dimension) ;
    }
  }

  public int getLessContributorHV(SolutionSet set) {

    Front wholeFront   = new Front();

    wholeFront.loadFront(set, -1);

    int index= 0;
    double contribution = Double.POSITIVE_INFINITY;

    for (int i = 0; i < set.size(); i++) {
      double [] v = new double[set.get(i).numberOfObjectives()];
      for (int j = 0; j < v.length; j++){
        v[j] = set.get(i).getObjective(j);
      }
      Point p = new Point(v);
      double aux = this.getExclusiveHV(wholeFront, i);
      if ((aux) < contribution) {
        index = i;
        contribution = aux;
      }
      set.get(i).setCrowdingDistance(aux);
    }

    return index;
  }

  public double getHV(Front front, Solution referencePoint) {
    referencePoint_ = new Point(referencePoint) ;
    double volume = 0.0 ;
    sort(front) ;

    if (currentDimension_ == 2)
      volume = get2DHV(front) ;
    else {
      volume = 0.0 ;

      currentDimension_ -- ;
      for (int i = front.nPoints_-1; i >= 0; i--) {
        volume += Math.abs(front.getPoint(i).objectives_[currentDimension_] -
                referencePoint_.objectives_[currentDimension_])*
                this.getExclusiveHV(front, i) ;
      }
      currentDimension_ ++ ;
    }

    return volume ;
  }

  public double getHV(Front front) {
    double volume = 0.0 ;
    if (front.nPoints_ == 0){
    	volume = 0.0;
	}else if(front.nPoints_ == 1){
		volume = hv2point(front.getPoint(0), referencePoint_);
	}else if (front.nPoints_ == 2){
		double [] mid = new double [front.getNumberOfObjectives()];
		for (int j=0;j<front.getNumberOfObjectives();j++){
			mid[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(1).objectives_[j]);
		}
		Point midp = new Point(mid);
		volume = hv2point(front.getPoint(0),referencePoint_)+hv2point(front.getPoint(1),referencePoint_)-hv2point(midp,referencePoint_);
		
	}else if (front.nPoints_==3){
		double [] w01 = new double [front.getNumberOfObjectives()];
		double [] w02 = new double [front.getNumberOfObjectives()];
		double [] w12 = new double [front.getNumberOfObjectives()];
		double [] w012 = new double [front.getNumberOfObjectives()];
		for (int j=0;j<front.getNumberOfObjectives();j++){
			w01[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(1).objectives_[j]);
			w02[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(2).objectives_[j]);
			w12[j] = Math.max(front.getPoint(1).objectives_[j], front.getPoint(2).objectives_[j]);
		}
		for (int j=0;j<front.getNumberOfObjectives();j++){
			w012[j] = Math.max(w02[j], front.getPoint(1).objectives_[j]);
		}
		Point p01 = new Point(w01);Point p02 = new Point(w02);Point p12 = new Point(w12);Point p012 = new Point(w012);
		volume = hv2point(front.getPoint(0),referencePoint_)+hv2point(front.getPoint(1),referencePoint_)+hv2point(front.getPoint(2),referencePoint_)
				-hv2point(p01,referencePoint_)-hv2point(p02,referencePoint_)-hv2point(p12,referencePoint_)+hv2point(p012,referencePoint_);
	}else{
	    sort(front) ;
	
	    if (currentDimension_ == 2)
	      volume = get2DHV(front) ;
	    else {
	      volume = 0.0 ;
	
	      currentDimension_ -- ;
	      for (int i = front.nPoints_-1; i >= 0; i--) {
	        volume += Math.abs(front.getPoint(i).objectives_[currentDimension_] -
	                referencePoint_.objectives_[currentDimension_])*
	                this.getExclusiveHV(front, i) ;
	      }
	      currentDimension_ ++ ;
	    }
	}
    return volume ;
  }

  public double get2DHV(Front front) {
    double hv = 0.0 ;

    hv = Math.abs((front.getPoint(0).getObjectives()[0] - referencePoint_.objectives_[0]) *
            (front.getPoint(0).getObjectives()[1] - referencePoint_.objectives_[1])) ;

    for (int i = 1; i < front.nPoints_; i++) {
      hv += Math.abs((front.getPoint(i).getObjectives()[0] - referencePoint_.objectives_[0]) *
              (front.getPoint(i).getObjectives()[1] - front.getPoint(i-1).getObjectives()[1])) ;

    }

    return hv ;

  }

  public double getInclusiveHV(Point p) {
    double volume = 1 ;
    for (int i = 0; i < currentDimension_; i++) {
      volume *= Math.abs(p.objectives_[i] - referencePoint_.objectives_[i]) ;
    }

    return volume ;
  }

  public double getExclusiveHV(Front front, int point) {
    double volume ;

    volume = getInclusiveHV(front.getPoint(point)) ;
    if (front.nPoints_ > point + 1) {
      makeDominatedBit(front, point);
      double v = getHV(fs_[currentDeep_-1]) ;
      volume -= v ;
      currentDeep_ -- ;
    }

    return volume ;
  }

  public void makeDominatedBit(Front front, int p) {
    int z = front.nPoints_ - 1 - p ;

    for (int i = 0 ; i < z ; i++)
      for (int j = 0 ; j < currentDimension_; j++) {
        fs_[currentDeep_].getPoint(i).objectives_[j] = worse(front.points_[p].objectives_[j], front.points_[p+1+i].objectives_[j], false) ;
      }


    Point t ;
    fs_[currentDeep_].nPoints_ = 1 ;

    for (int i = 1; i < z; i++) {
      int j = 0 ;
      boolean keep = true ;
      while (j < fs_[currentDeep_].nPoints_ && keep) {
        switch (dominates2way(fs_[currentDeep_].points_[i], fs_[currentDeep_].points_[j])) {
          case -1:
            t = fs_[currentDeep_].points_[j] ;
            fs_[currentDeep_].nPoints_--;
            fs_[currentDeep_].points_[j] = fs_[currentDeep_].points_[fs_[currentDeep_].nPoints_];
            fs_[currentDeep_].points_[fs_[currentDeep_].nPoints_] = t;
            break;
          case  0: j++; break;
          // case  2: printf("Identical points!\n");
          default: keep = false;
        }
      }
      if (keep) {t = fs_[currentDeep_].points_[fs_[currentDeep_].nPoints_];
        fs_[currentDeep_].points_[fs_[currentDeep_].nPoints_] = fs_[currentDeep_].points_[i];
        fs_[currentDeep_].points_[i] = t;
        fs_[currentDeep_].nPoints_++;
      }
    }

    currentDeep_++ ;
  }

  private double worse (double x, double y, boolean maximizing) {
    double result ;
    if (maximizing) {
      if (x > y)
        result = y ;
      else
        result = x ;
    }
    else {
      if (x > y)
        result = x ;
      else
        result = y ;
    }
    return result ;
  }

  int dominates2way(Point p, Point q)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise
  // ASSUMING MINIMIZATION
  {
    // domination could be checked in either order

    for (int i = currentDimension_ - 1; i >= 0; i--)
      if (p.objectives_[i] < q.objectives_[i]){
        for (int j = i - 1; j >= 0; j--)
          if (q.objectives_[j] < p.objectives_[j]) return 0;
        return -1;
      }
      else
      if (q.objectives_[i] < p.objectives_[i]){
        for (int j = i - 1; j >= 0; j--)
          if (p.objectives_[j] < q.objectives_[j]) return 0;
        return  1;
      }
    return 2;
  }

  public void sort(Front front) {
    Arrays.sort(front.points_, 0, front.nPoints_, pointComparator_);
  }
  public static double hv2point(Point Point1, Point ref){
  double x=ref.objectives_[0]-Point1.objectives_[0];
	for (int j=1;j<Point1.getNumberOfObjectives();j++){
		 x = x*(ref.objectives_[j]-Point1.objectives_[j]);
	}
	return x;
}
  public static void main(String args[]) throws IOException {
    Front front = new Front() ;
    
    Point referencePoint =null;
    String problem_name ="";
    int m =2;
    for (int i = 13; i <=13; i++){
    	//get problem name bounds and objective number m
    	if (i<5){ //ZDT1-ZDT4
    		problem_name = String.format("ZDT%d",i);
    		referencePoint = new Point(new double [] {2.0, 2.0});
    	}else if (i==5){
    		problem_name = "ZDT6";
    		referencePoint = new Point(new double [] {2.0, 2.0});
    	}else if (i <= 12){ //DTLZ
    		m=10;
    		problem_name = String.format("DTLZ%d", i-5);
    		if (i==6){//DTLZ1
    			double [] tempp = new double [m];
    			for (int j=0;j<m;j++){
    				tempp[j] = 1.0;
    			}
    			referencePoint = new Point(tempp);
    		}else if (i<12){
    			double [] tempp = new double [m];
    			for (int j=0;j<m;j++){
    				tempp[j] = 2.0;
    			}
    			referencePoint = new Point(tempp);
    		}else{
    			double [] tempp = new double [m];
    			for (int j=0;j<m-1;j++){
    				tempp[j] = 2.0;
    			}
    			tempp[m] = 7.0;
    			referencePoint = new Point(tempp);
    		}
    	}else if (i <= 21){
    		m=10;
    		problem_name = String.format("WFG%d", i-12);
			double [] tempp = new double [m];
			for (int j=0;j<m;j++){
				tempp[j] = 2.0*(j+1)+1;
			}
			referencePoint = new Point(tempp);
    	}
    }

    double hv [] = new double [30];
    for (int i=0;i<30;i++){
		front.readFront(String.format("NSGA-II_SBX\\NSGAII_%s_%d_T%d",problem_name,m,i+1), referencePoint);
	
//		//NORMALIZATION
//		for (int j=0;j<front.nPoints_;j++)
//			for(int k=0;k<front.getNumberOfObjectives();k++)
//				front.getPoint(j).objectives_[k] = front.getPoint(j).objectives_[k]/referencePoint.objectives_[k];
//		for (int j=0;j<front.getNumberOfObjectives();j++)
//			referencePoint.objectives_[j] = 1.0;
		if (front.nPoints_ == 0){
			hv[i] = 0.0;
		}else if(front.nPoints_ == 1){
			hv[i] = hv2point(front.getPoint(0), referencePoint);
		}else if (front.nPoints_ == 2){
			double [] mid = new double [front.getNumberOfObjectives()];
			for (int j=0;j<front.getNumberOfObjectives();j++){
				mid[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(1).objectives_[j]);
			}
			Point midp = new Point(mid);
			hv[i] = hv2point(front.getPoint(0),referencePoint)+hv2point(front.getPoint(1),referencePoint)-hv2point(midp,referencePoint);
			
		}else if (front.nPoints_==3){
			double [] w01 = new double [front.getNumberOfObjectives()];
			double [] w02 = new double [front.getNumberOfObjectives()];
			double [] w12 = new double [front.getNumberOfObjectives()];
			double [] w012 = new double [front.getNumberOfObjectives()];
			for (int j=0;j<front.getNumberOfObjectives();j++){
				w01[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(1).objectives_[j]);
				w02[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(2).objectives_[j]);
				w12[j] = Math.max(front.getPoint(1).objectives_[j], front.getPoint(2).objectives_[j]);
			}
			for (int j=0;j<front.getNumberOfObjectives();j++){
				w012[j] = Math.max(w02[j], front.getPoint(1).objectives_[j]);
			}
			Point p01 = new Point(w01);Point p02 = new Point(w02);Point p12 = new Point(w12);Point p012 = new Point(w012);
			hv[i] = hv2point(front.getPoint(0),referencePoint)+hv2point(front.getPoint(1),referencePoint)+hv2point(front.getPoint(2),referencePoint)
					-hv2point(p01,referencePoint)-hv2point(p02,referencePoint)-hv2point(p12,referencePoint)+hv2point(p012,referencePoint);
		}
		else{
			WFGHV wfghv = new WFGHV(referencePoint.getNumberOfObjectives(), front.getNumberOfPoints(), referencePoint) ;
		    hv[i] = wfghv.getHV(front);
		}
    }
    printGD(String.format("NSGA-II_SBX\\NSGAII_%s_%d_HV",problem_name,m),hv);
  }

}
