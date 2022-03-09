package jmetal.problems.LMF;

public class GFunctions {
	
	/*LMF1: GFunction used in LMF1*/
	public static double[] getGLMF1(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}
				}
			}
		}
		
		return g;
	}
	
	/*LMF2: GFunction used in LMF2*/
	public static double[] getGLMF2(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}
					
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getAckley(xII[i][j]);
					}
				}
			}
		}
		
		return g;
	}
	
	/*LMF3: GFunction used in LMF3*/
	public static double[] getGLMF3(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getAckley(xII[i][j]);
					}
					
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}
				}
			}
		}
		
		return g;
	}
	
	/*LMF4: GFunction used in LMF4*/
	public static double[] getGLMF4(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}
					
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRosenbrock(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	/*LMF5: GFunction used in LMF5*/
	public static double[] getGLMF5(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w1[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}
					
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w2[j] * BaseFunctions.getSphere(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w2[j] * BaseFunctions.getAckley(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	/*LMF6: GFunction used in LMF6*/
	public static double[] getGLMF6(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w1[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getAckley(xII[i][j]);
					}
					
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w2[j] * BaseFunctions.getSphere(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w2[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	/*LMF7: GFunction used in LMF7*/
	public static double[] getGLMF7(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w1[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getRosenbrock(xII[i][j]);
					}
					
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w2[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w2[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRosenbrock(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	
	/*LMF8: GFunction used in LMF8*/
	public static double[] getGLMF8(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w1[j] * BaseFunctions.getAckley(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getRosenbrock(xII[i][j]);
					}
					
				}
			}else{
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 3 == 0) {
						g[i] += w2[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}else if(j % 3 == 1) {
						g[i] += w2[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRosenbrock(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	/*LMF9: GFunction used in LMF9*/
	public static double[] getGLMF9(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}
					
				}
			}else if(i%3 == 1){
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}
				}
			}else {
				double[] w3 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w3 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w3[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}else {
						g[i] += w3[j] * BaseFunctions.getAckley(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	/*LMF10: GFunction used in LMF10*/
	public static double[] getGLMF10(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getAckley(xII[i][j]);
					}
					
				}
			}else if(i%3 == 1){
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}
				}
			}else {
				double[] w3 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w3 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w3[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w3[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	/*LMF11: GFunction used in LMF11*/
	public static double[] getGLMF11(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getAckley(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}
					
				}
			}else if(i%3 == 1){
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getAckley(xII[i][j]);
					}
				}
			}else {
				double[] w3 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w3 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w3[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w3[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}
				}
			}
		}
		return g;
	}
	
	/*LMF12: GFunction used in LMF12*/
	public static double[] getGLMF12(double x1, double xII[][][],int[][][] index, int type) {
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				double[] w1 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w1 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w1[j] * BaseFunctions.getSphere(xII[i][j]);
					}else {
						g[i] += w1[j] * BaseFunctions.getRosenbrock(xII[i][j]);
					}
					
				}
			}else if(i%3 == 1){
				double[] w2 = getWeightsLogisticMap(0.23, 3.75, xII[i].length);
				if(type == 0) {
					w2 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w2[j] * BaseFunctions.getSchwefel2(xII[i][j]);
					}else {
						g[i] += w2[j] * BaseFunctions.getRastrigin(xII[i][j]);
					}
				}
			}else {
				double[] w3 = getWeightsLogisticMap(0.23, 3.7, xII[i].length);
				if(type == 0) {
					w3 = getUniformWeights(xII[i].length);
				}
				for(int j=0;j<xII[i].length;j++){
					if(j % 2 == 0) {
						g[i] += w3[j] * BaseFunctions.getSchwefel1(xII[i][j]);
					}else {
						g[i] += w3[j] * BaseFunctions.getAckley(xII[i][j]);
					}
				}
			}
		}
		return g;
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
	
	public static double[] getUniformWeights(int num) {
		double[] w = new double[num];
		for(int i = 0; i<num; i++) {
			w[i] = 1.0/num;
		}
		return w;
	}

}
