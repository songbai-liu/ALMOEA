package jmetal.metaheuristics.almoea;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.almoea.BaseMLP.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.wrapper.XReal;

public class MLPModel {
	
	int features;
	int feaAfterEncode;
	int layerNum;
	double learningRate;
	BPMLP baseMLP;
	
	
	public MLPModel (int numVars, int reD, int layerNum, double rate) {
		features = numVars;
		feaAfterEncode = reD;
		this.layerNum = layerNum;
		learningRate = rate;
		baseMLP = new BPMLP(features,feaAfterEncode,layerNum,learningRate);
	}
	
	public void getTrainingModel(SolutionSet trainingSet, SolutionSet targetSet, double[] ub, double[] lb) throws JMException {
		int size = trainingSet.size();
		int tSize = targetSet.size();
		if(size != tSize){
			System.out.println("the size of the input is not equal to the size of labels");
			System.exit(0);
		}
		double[][] inputs = new double[size][features];
		double[][] targets = new double[size][features];
		for(int i=0;i<size;i++) {
			//double fit1 = trainingSet.get(i).getSumValue();
			double fit1 = trainingSet.get(i).getDistanceToIdealPoint();
			//double fit2 = targetSet.get(i).getSumValue();
			double fit2 = targetSet.get(i).getDistanceToIdealPoint();
			XReal sol1;
			XReal sol2;
			if(fit1 > fit2) {
				sol1 = new XReal(trainingSet.get(i));
				sol2 = new XReal(targetSet.get(i));
			}else {
				sol2 = new XReal(trainingSet.get(i));
				sol1 = new XReal(targetSet.get(i));
			}
			
			for(int j=0;j<features;j++) {
				inputs[i][j] = (sol1.getValue(j)-lb[j])/(ub[j]-lb[j]);
				targets[i][j] = (sol2.getValue(j)-lb[j])/(ub[j]-lb[j]);
			}
		}
		baseMLP.trainModel(inputs, targets, 1);
	}
	
	public void getTrainingModel(SolutionSet trainingSet, SolutionSet targetSet, int epochs, int[] group) throws JMException {
		int size = trainingSet.size();
		int tSize = targetSet.size();
		double[][] inputs = new double[size][features];
		double[][] targets = new double[size][features];
		for(int i=0;i<size;i++) {
			XReal sol = new XReal(trainingSet.get(i));
			XReal tar;
			if(tSize < size) {
				tar = new XReal(targetSet.get(PseudoRandom.randInt(0, tSize-1)));
			}else {
				tar = new XReal(targetSet.get(i));
			}
			for(int j=0;j<features;j++) {
				inputs[i][j] = sol.getValue(group[j]);
				targets[i][j] = tar.getValue(group[j]);
			}
		}
		baseMLP.trainModel(inputs, targets, epochs);
	}
	
	public double[] encode(Solution sol, int dim) throws JMException {
		double[] encodedSolution = new double[feaAfterEncode];
		XReal xsol = new XReal(sol);
		if(dim != features) {
			System.out.println("The variable-related dimensions of the input solution do not match the model at encode a single solution!");
			System.out.println("features = " + features + ", and the dimensions of input is: " + dim);
			System.exit(0);
		}else {
			double[] input = new double[features];
			for(int i=0;i<features;i++) { 
				input[i] = xsol.getValue(i); 
			}
			baseMLP.computeEncodeData(input, encodedSolution);
		}
		return encodedSolution;
	}
	
	public double[][] encode(SolutionSet solSet,int dim) throws JMException{
		int size = solSet.size();
		double[][] encodedSet = new double[size][feaAfterEncode];
		for(int p=0;p<size;p++) {
			XReal xsol = new XReal(solSet.get(p));
			if(dim != features) {
				System.out.println("The variable-related dimensions of the input solution do not match the model at encode a solution set!");
				System.out.println("features = " + features + ", and the dimensions of input is: " + dim);
				System.exit(0);
			}else {
				double[] input = new double[features];
				for(int i=0;i<features;i++) { 
					input[i] = xsol.getValue(i); 
				}
				baseMLP.computeEncodeData(input, encodedSet[p]);
			}
		}
		return encodedSet;
	}
	
	public double[] encode(double[] input) {
		double[] encodedSolution = new double[feaAfterEncode];
		if(input.length != features) {
			System.out.println("The variable-related dimensions of the input solution do not match the model at encode a array of input!");
			System.exit(0);
		}else {
			baseMLP.computeEncodeData(input, encodedSolution);
		}
		return encodedSolution;
	}
	
	public double[] decode(double[] encodeData) {
		double[] decodedSolution = new double[features];
		if(encodeData.length != feaAfterEncode) {
			System.out.println("The dimensions of the input encoded data do not match the model at decode a arry of input!");
			System.exit(0);
		}else {
			baseMLP.computeDecodeData(encodeData, decodedSolution);
		}
		return decodedSolution;
	}

}
