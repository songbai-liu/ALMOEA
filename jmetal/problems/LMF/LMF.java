package jmetal.problems.LMF;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

public abstract class LMF extends Problem {
	String gType_;
	int type_fm;//type for formulation model, setting 0 for addition, 1 for multiplication, while 2 for mixed model
	int type_lk;//type of variable linkage, 0 for linear linkage, 1 for nonlinear linkage, and 2 for mixed/hybrid
	int type_dg;//type for deep grouping distance-related variables, setting 0 for even grouping while 1 for nonuniform grouping
	int type_cv;//type for contribution of variables, setting 0 for balanced contribution while 1 for unbalanced contribution
	
	int[] nip;//Number of independent variables in each position-related variable group
    int[] nop;//Number of overlapping variables in each position-related variable group
    int nsp;//Number of shared variables in each position-related variable group
    int[] nid;//Number of independent variables in each distance-related variable group
    int[] nod;//Number of overlapping variables in each distance-related variable group
    int nsd;// Number of shared variables in each distance-related variable group
    int K;//Number of position-related variables
    int L;//Number of distance-related variables
    int[][] gp;
    int[][] gd;
    int[][][] dgd;

	public LMF(String solutionType, int numberOfVariables,int numberOfObjectives) {
		numberOfObjectives_ = numberOfObjectives;
		numberOfVariables_ = numberOfVariables;
		gType_ = "LMF1";
		hType_ = "concave";
		type_fm = 2;
        type_lk = 2;
        type_dg = 1;
        type_cv = 1;
		
		nip = new int[numberOfObjectives_-1];
		nop = new int[numberOfObjectives_-1];
		int sumNIP = 0;
		for(int i=0; i<numberOfObjectives_-1; i++){
			nip[i] = 3;
			if(numberOfObjectives_ == 2){
				nop[i] = 0;
			}else{
				nop[i] = 1;
			}
			sumNIP += nip[i];
		}
		nsp = 1;
		K = nsp + sumNIP;
		L = numberOfVariables_ - K;
		
		double[] c = new double[numberOfObjectives_];
		c[0] = 3.8*0.1*(1.0-0.1);
		double sumC = c[0];
		for(int i = 1; i < numberOfObjectives_; i++){
			c[i] = 3.8*c[i-1]*(1.0-c[i-1]);
			sumC += c[i];
		}
		nid = new int[numberOfObjectives_];
		int sumNID = 0;
		int nsd_leas = 2;//make sure the number of nsd is at least 2
		int proportion = 5;//control the number of nod as a percentage of the number of nid
		for(int i = 0; i < numberOfObjectives_; i++){
			nid[i] = (int)Math.floor((c[i]/sumC)*(L - nsd_leas));
			sumNID += nid[i];
		}
		nod = new int[numberOfObjectives_];
		for(int i = 0; i < numberOfObjectives_; i++){
			if(i == 0){
				nod[i] = nid[numberOfObjectives_-1]/proportion;
			}else{
				nod[i] = nid[i-1]/proportion;
			}
		}
		nsd = L - sumNID;
		
		gp = GroupingStrategy.overlapGrouping(nip, nop, nsp, numberOfObjectives_-1, 0);
		gd = GroupingStrategy.overlapGrouping(nid, nod, nsd, numberOfObjectives_, K);
		dgd = new int[numberOfObjectives_][][];
		for(int i=0;i<numberOfObjectives_;i++){
			if(i%2 == 0){ 
				dgd[i] = GroupingStrategy.deepGrouping(gd[i], 5, 1); 
			}else{
				dgd[i] = GroupingStrategy.deepGrouping(gd[i], 5, 2); 
			}
			//dgd[i] = GroupingStrategy.uniformDeepGrouping(gd[i], 5);
		}
		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < K; var++) {
			lowerLimit_[var] = -1.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = K; var < numberOfVariables; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 10.0;
		}
		
		if (solutionType.compareTo("BinaryReal") == 0)
			solutionType_ = new BinaryRealSolutionType(this);
		else if (solutionType.compareTo("Real") == 0)
			solutionType_ = new RealSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType
					+ " invalid");
			System.exit(-1);
		}
	}
	
    //type1: type of formulation model, 0 for addition, 1 for multiplication, and 2 for mixed/hybrid
    //type2: type of variable linkage, 0 for linear linkage, 1 for nonlinear linkage, and 2 for mixed/hybrid
	//type3: type of deep grouping on variables, 0 for even grouping, 1 for nonuniform grouping
	//type4: type of contribution on variables, 0 for balanced contribution, 1 for unbalanced contribution
	public LMF(String solutionType, int numberOfVariables,int numberOfObjectives, int type1, int type2, int type3, int type4) {
		numberOfObjectives_ = numberOfObjectives;
		numberOfVariables_ = numberOfVariables;
		gType_ = "LMF1";
		hType_ = "concave";
		type_fm = type1;
        type_lk = type2;
        type_dg = type3;
        type_cv = type4;
		
		nip = new int[numberOfObjectives_-1];
		nop = new int[numberOfObjectives_-1];
		int sumNIP = 0;
		for(int i=0; i<numberOfObjectives_-1; i++){
			nip[i] = 3;
			if(numberOfObjectives_ == 2){
				nop[i] = 0;
			}else{
				nop[i] = 1;
			}
			sumNIP += nip[i];
		}
		nsp = 1;
		K = nsp + sumNIP;
		L = numberOfVariables_ - K;
		
		double[] c = new double[numberOfObjectives_];
		c[0] = 3.8*0.1*(1.0-0.1);
		double sumC = c[0];
		for(int i = 1; i < numberOfObjectives_; i++){
			c[i] = 3.8*c[i-1]*(1.0-c[i-1]);
			sumC += c[i];
		}
		nid = new int[numberOfObjectives_];
		int sumNID = 0;
		int nsd_leas = 2;//make sure the number of nsd is at least 2
		int proportion = 5;//control the number of nod as a percentage of the number of nid
		for(int i = 0; i < numberOfObjectives_; i++){
			nid[i] = (int)Math.floor((c[i]/sumC)*(L - nsd_leas));
			sumNID += nid[i];
		}
		nod = new int[numberOfObjectives_];
		for(int i = 0; i < numberOfObjectives_; i++){
			if(i == 0){
				nod[i] = nid[numberOfObjectives_-1]/proportion;
			}else{
				nod[i] = nid[i-1]/proportion;
			}
		}
		nsd = L - sumNID;
		
		gp = GroupingStrategy.overlapGrouping(nip, nop, nsp, numberOfObjectives_-1, 0);
		gd = GroupingStrategy.overlapGrouping(nid, nod, nsd, numberOfObjectives_, K);
		dgd = new int[numberOfObjectives_][][];
		for(int i=0;i<numberOfObjectives_;i++){
			if(type_dg == 1) {//%nonuniformly deep grouping without knowing the number of groups
				if(i%2 == 0){ 
					dgd[i] = GroupingStrategy.deepGrouping(gd[i], 5, 1); 
				}else{
					dgd[i] = GroupingStrategy.deepGrouping(gd[i], 5, 2); 
				}
			}else {
				dgd[i] = GroupingStrategy.uniformDeepGrouping(gd[i], 5);
			}
		}
		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < K; var++) {
			lowerLimit_[var] = -1.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = K; var < numberOfVariables; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 10.0;
		}
		
		if (solutionType.compareTo("BinaryReal") == 0)
			solutionType_ = new BinaryRealSolutionType(this);
		else if (solutionType.compareTo("Real") == 0)
			solutionType_ = new RealSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType
					+ " invalid");
			System.exit(-1);
		}
	}
	
	double[] evalG(double xI, double[][][] xIII) throws JMException {
		double[] g = null;
		if (gType_.equalsIgnoreCase("LMF1"))
			g = GFunctions.getGLMF1(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF2"))
			g = GFunctions.getGLMF2(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF3"))
			g = GFunctions.getGLMF3(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF4"))
			g = GFunctions.getGLMF4(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF5"))
			g = GFunctions.getGLMF5(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF6"))
			g = GFunctions.getGLMF6(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF7"))
			g = GFunctions.getGLMF7(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF8"))
			g = GFunctions.getGLMF8(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF9"))
			g = GFunctions.getGLMF9(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF10"))
			g = GFunctions.getGLMF10(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF11"))
			g = GFunctions.getGLMF11(xI, xIII,dgd,type_cv);
		else if(gType_.equalsIgnoreCase("LMF12"))
			g = GFunctions.getGLMF12(xI, xIII,dgd,type_cv);
		else {
			System.out.println("Error: g function type " + gType_ + " invalid");
			System.exit(0);
		}
		return g;
	}

	double[] evalH(double[] xI) {
		double[] h = new double[numberOfObjectives_];
		if (hType_.equalsIgnoreCase("linear")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= xI[j];
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= (1 - xI[aux]);
				} // if
			} // for
		}else if(hType_.equalsIgnoreCase("concave")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
			} // for
		}else if(hType_.equalsIgnoreCase("convex")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
				if(i != numberOfObjectives_-1)
					h[i] = Math.pow(h[i], 4);
				else
					h[i] = Math.pow(h[i], 2);
			} // for
		}if (hType_.equalsIgnoreCase("inverted_linear")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= xI[j];
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= (1 - xI[aux]);
				} // if
				h[i] = 1.0 - h[i];
			} // for
		}else if(hType_.equalsIgnoreCase("inverted_concave")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
				h[i] = 1.0 - h[i];
			} // for
		}
		return h;
	}

	protected double sum_avg_abs(double[] x){
	    double avg = 0;
		for(int i=0; i<x.length; i++){
			avg += x[i];
		}
		return Math.abs(avg)/x.length;
	}
	
}
