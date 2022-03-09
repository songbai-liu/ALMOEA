package jmetal.qualityIndicator.fastHypervolume.wfg;

import java.io.IOException;

import jmetal.core.Problem;
import jmetal.core.*;

public class wfghvCalculator2 {
	public jmetal.qualityIndicator.util.MetricsUtil utils_;
	//SolutionSet pf_ = null;
	double [][] pfMatrix_ = null;
	int fun;
	public wfghvCalculator2(double[][] population) {
	    pfMatrix_ = population;  
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor 
	
	public static double hv2point(double[] Point1, Solution ref){
		  double x=ref.getObjective(0)-Point1[0];
			for (int j=1;j<Point1.length;j++){
				 x = x*(ref.getObjective(j)-Point1[j]);
			}
			return x;
		}
	public double calculatewfghv() throws IOException{

	    double hv;
	    int number = pfMatrix_[0].length;
	    
	    Solution referencePoint1 = new Solution( number);
		
		for (int j=0;j<number;j++){
			referencePoint1.setObjective(j,1.5);
		}
		
		//NORMALIZATION
		//for (int j=0;j<pfMatrix_.length;j++)
			//for(int k=0;k<number;k++)
				//pfMatrix_[j][k] = pfMatrix_[j][k]/(1.1*referencePoint1.getObjective(k));
		//invertedFront = utils_.invertedFront(pf_,number);
		for (int j=0;j<pfMatrix_.length;j++)
			for(int k=0;k<number;k++)
				if(pfMatrix_[j][k]>1.5){
					//pf_.remove(j);
					//j--;
					for(int s=0;s<number;s++){
						pfMatrix_[j][s] = 1.5;
					}
					break;
				}
		
		if (pfMatrix_.length == 0){
			hv = 0.0;
		}else if(pfMatrix_.length == 1){
			hv = hv2point(pfMatrix_[0], referencePoint1);
		}else if (pfMatrix_.length == 2){
			double [] mid = new double [number];
			for (int j=0;j<number;j++){
				mid[j] = Math.max(pfMatrix_[0][j], pfMatrix_[1][j]);
			}
			double[] midp = new double[number];
			for(int i=0;i<number;i++){
				midp[i] = mid[i];
			}
			hv = hv2point(pfMatrix_[0],referencePoint1)+hv2point(pfMatrix_[1],referencePoint1)-hv2point(midp,referencePoint1);
			
		}else if (pfMatrix_.length==3){
			double [] w01 = new double [number];
			double [] w02 = new double [number];
			double [] w12 = new double [number];
			double [] w012 = new double [number];
			for (int j=0;j<number;j++){
				w01[j] = Math.max(pfMatrix_[0][j], pfMatrix_[1][j]);
				w02[j] = Math.max(pfMatrix_[0][j], pfMatrix_[2][j]);
				w12[j] = Math.max(pfMatrix_[1][j], pfMatrix_[2][j]);
			}
			for (int j=0;j<number;j++){
				w012[j] = Math.max(w02[j], pfMatrix_[1][j]);
			}
			double[] p01 = new double[number];double[] p02 = new double[number];double[] p12 = new double[number];double[] p012 = new double[number];
			for(int i=0;i<number;i++){
				p01[i] = w01[i];
				p02[i] = w02[i];
				p12[i] = w12[i];
				p012[i] = w012[i];
			}
			hv = hv2point(pfMatrix_[0],referencePoint1)+hv2point(pfMatrix_[1],referencePoint1)+hv2point(pfMatrix_[2],referencePoint1)
					-hv2point(p01,referencePoint1)-hv2point(p02,referencePoint1)-hv2point(p12,referencePoint1)+hv2point(p012,referencePoint1);
		}else{
			WFGHV1 wfghv = new WFGHV1(number, pfMatrix_.length) ;
			Front1 front = new Front1(pfMatrix_.length, number, pfMatrix_);
		    hv = wfghv.getHV(front,referencePoint1);
		}
		return hv;
  } // CalculateHypervolume
}
