package jmetal.qualityIndicator.fastHypervolume.wfg;

import java.io.IOException;

import jmetal.core.Problem;
import jmetal.core.*;

public class wfghvCalculator1 {
	public jmetal.qualityIndicator.util.MetricsUtil utils_;
	SolutionSet pf_ = null;
	double [][] pfMatrix_ = null;
	int fun;
	public wfghvCalculator1(SolutionSet paretoFront) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;  
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor 
	
	public wfghvCalculator1(SolutionSet paretoFront,int fun) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;
	    this.fun = fun;
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor
	public static double hv2point(Solution Point1, Solution ref){
		  double x=ref.getObjective(0)-Point1.getObjective(0);
			for (int j=1;j<Point1.numberOfObjectives();j++){
				 x = x*(ref.getObjective(j)-Point1.getObjective(j));
			}
			return x;
		}
	public double calculatewfghv() throws IOException{
		SolutionSet sb = new SolutionSet(pf_.size());
		for(int ss=0;ss<pf_.size();ss++){
			Solution sss = new Solution(pf_.get(ss));
			sb.add(sss);
		}

	    double hv;
	    int number = pf_.get(0).numberOfObjectives();
	    
	    Solution referencePoint1 = new Solution( number);
		
		for (int j=0;j<number;j++){
			if(fun == 6){//DTLZ1
				referencePoint1.setObjective(j,0.5);
			}else if(fun<=9){//DTLZ2-DTLZ4
				referencePoint1.setObjective(j,1.0);
			}else if(fun>9&&fun<=11){//DTLZ5 and DTLZ6
				if(j!=number-1){
					referencePoint1.setObjective(number-1-j,Math.pow(Math.sqrt(2)/2, j));
				}else{
					referencePoint1.setObjective(0,referencePoint1.getObjective(1));
				}
			}else if(fun == 12){//MaF7: DTLZ7
				if(j!=number-1){
					referencePoint1.setObjective(j,1.0);
				}else{
					referencePoint1.setObjective(j,2.0*(j+1));
				}
			}else if(fun>12&&fun<=21){//WFG1-WFG9
				referencePoint1.setObjective(j,2.0*(j+1));
			}else if(fun == 22){//MaF1: Inverted DTLZ1
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 23){//MaF2: Concave DTLZ2-BZ
				if(j!=number-1){
					referencePoint1.setObjective(number-1-j,Math.pow(Math.cos(Math.PI/8), j+1));
				}else{
					referencePoint1.setObjective(0,referencePoint1.getObjective(1));
				}
			}else if(fun == 24){//MaF3: Convex DTLZ3
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 25){//MaF4: Inverted badly-scaled DTLZ3
				referencePoint1.setObjective(j,Math.pow(2.0, j+1));
			}else if(fun == 26){//MaF5: Concave-badly scaled DTLZ4
				referencePoint1.setObjective(number-1-j,Math.pow(2.0, j+1));
			}else if(fun == 27){//MaF6: DTLZ5
				if(j!=number-1){
					referencePoint1.setObjective(number-1-j,Math.pow(Math.cos(Math.PI/4), j));
				}else{
					referencePoint1.setObjective(0,referencePoint1.getObjective(1));
				}
			}else if(fun == 28){//MaF7: DTLZ7
				if(j!=number-1){
					referencePoint1.setObjective(j,1.0);
				}else{
					referencePoint1.setObjective(j,2.0*(j+1));
				}
			}else if(fun == 29){//MaF8: Multi-Point Distance Minimization Problem
				double[][] point = new double[number][2];
				point[0][0] = 0.0;
				point[0][1] = 1.0;	
				double arc = 2*Math.PI/number;
				for (int i = 1; i < number; i++){
					point[i][0] = point[0][0] - Math.sin(arc*i);
					point[i][1] = point[0][1] - 1.0 + Math.cos(arc*i);
				}
				double maxValue = Double.MIN_VALUE;
				for(int i=0;i<number;i++){
					for(int s=i+1;s<number;s++){
						double value = Math.pow(point[i][0]-point[j][0], 2) + Math.pow(point[i][1]-point[j][1], 2);
						if(value > maxValue){
							maxValue = value;
						}
					}	
				}
				referencePoint1.setObjective(j,Math.sqrt(maxValue));
			}else if(fun == 30){//MaF9: Multi-Line Distance Minimization Problem
				double[][] point = new double[number][2];
				point[0][0] = 0.0;
				point[0][1] = 1.0;	
				double arc = 2*Math.PI/number;
				for (int i = 1; i < number; i++){
					point[i][0] = point[0][0] - Math.sin(arc*i);
					point[i][1] = point[0][1] - 1.0 + Math.cos(arc*i);
				}
				
				double[] k = new double[number];
				double[] f = new double[number];
				for (int i = 0; i < number-1; i++){
					k[i] = (point[i+1][1]-point[i][1])/(point[i+1][0]-point[i][0]);
					f[i] = Math.abs(point[0][1]-k[i]*point[0][0]+k[i]*point[i][0]-point[i][1])/Math.sqrt(1+Math.pow(k[i], 2));
				}
				k[number-1]=(point[number-1][1]-point[0][1])/(point[number-1][0]-point[0][0]);
				f[number-1] = Math.abs(point[0][1]-k[number-1]*point[0][0]+k[number-1]*point[0][0]-point[0][1])
						/Math.sqrt(1+Math.pow(k[number-1], 2));
				
				double maxValue = Double.MIN_VALUE;
				for(int i=0;i<number;i++){
					if(f[i] > maxValue){
						maxValue = f[i];
					}	
				}
				referencePoint1.setObjective(j,Math.sqrt(maxValue));
			}else if(fun>=31&&fun<=33){//MaF10-MaF12: WFG1, WFG2, and WFG9
				referencePoint1.setObjective(j,2.0*(j+1));
			}else if(fun == 34){//MaF13: PF7
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 35){//mDTLZ1
			    referencePoint1.setObjective(j,0.5);
			}else if(fun>=36&&fun<=38){//mDTLZ2-mDTLZ4
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 39){//
				
			}
		}
		
		//NORMALIZATION
		for (int j=0;j<sb.size();j++)
			for(int k=0;k<number;k++)
				sb.get(j).setObjective(k, sb.get(j).getObjective(k)/(1.1*referencePoint1.getObjective(k)) );
		        //sb.get(j).setObjective(k, sb.get(j).getObjective(k)/(referencePoint1.getObjective(k)) );
		//SolutionSet invertedFront;
		//invertedFront = utils_.invertedFront(pf_,number);
		for (int j=0;j<sb.size();j++)
			for(int k=0;k<number;k++)
				if(sb.get(j).getObjective(k)>1.0){
					/*sb.remove(j);
					j--;*/
					for(int s=0;s<number;s++){
					   sb.get(j).setObjective(s, 1.0);
					}
					break;
					//sb.get(j).setObjective(k, 1.0);
				}
		for (int j=0;j<number;j++)
			//referencePoint1.setObjective(j,2.0*j+10);
		    referencePoint1.setObjective(j,1.0);
		    //referencePoint1.setObjective(j,1.1);
		if (sb.size() == 0){
			hv = 0.0;
		}else if(sb.size() == 1){
			hv = hv2point(sb.get(0), referencePoint1);
		}else if (sb.size() == 2){
			double [] mid = new double [number];
			for (int j=0;j<number;j++){
				mid[j] = Math.max(sb.get(0).getObjective(j), sb.get(1).getObjective(j));
			}
			Solution midp = new Solution(number);
			for(int i=0;i<number;i++){
				midp.setObjective(i, mid[i]);
			}
			hv = hv2point(sb.get(0),referencePoint1)+hv2point(sb.get(1),referencePoint1)-hv2point(midp,referencePoint1);
			
		}else if (sb.size()==3){
			double [] w01 = new double [number];
			double [] w02 = new double [number];
			double [] w12 = new double [number];
			double [] w012 = new double [number];
			for (int j=0;j<number;j++){
				w01[j] = Math.max(sb.get(0).getObjective(j), sb.get(1).getObjective(j));
				w02[j] = Math.max(sb.get(0).getObjective(j), sb.get(2).getObjective(j));
				w12[j] = Math.max(sb.get(1).getObjective(j), sb.get(2).getObjective(j));
			}
			for (int j=0;j<number;j++){
				w012[j] = Math.max(w02[j], sb.get(1).getObjective(j));
			}
			Solution p01 = new Solution(number);Solution p02 = new Solution(number);Solution p12 = new Solution(number);Solution p012 = new Solution(number);
			for(int i=0;i<number;i++){
				p01.setObjective(i, w01[i]);
				p02.setObjective(i, w02[i]);
				p12.setObjective(i, w12[i]);
				p012.setObjective(i, w012[i]);
			}
			hv = hv2point(sb.get(0),referencePoint1)+hv2point(sb.get(1),referencePoint1)+hv2point(sb.get(2),referencePoint1)
					-hv2point(p01,referencePoint1)-hv2point(p02,referencePoint1)-hv2point(p12,referencePoint1)+hv2point(p012,referencePoint1);
		}else{
			WFGHV1 wfghv = new WFGHV1(number, sb.size()) ;
			Front1 front = new Front1(sb.size(), number, sb);
		    hv = wfghv.getHV(front,referencePoint1);
		}
		return hv;
  } // CalculateHypervolume
}
