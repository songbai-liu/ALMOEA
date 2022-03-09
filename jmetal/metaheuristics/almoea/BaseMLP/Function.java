package jmetal.metaheuristics.almoea.BaseMLP;

public class Function {
	public static double sigmoid(double z) {
		double f = 0;
		f = 1.0/(1+Math.exp(-z));
		return f;
	}
	
	public static double inverseSigmoid(double z) {
		double f = 0;
		f = Math.log(z/(1-z));
		return f;
	}
	
	public static double tanh(double z) {
		double f = 0;
		f = (Math.exp(z) - Math.exp(-z))/(Math.exp(z) + Math.exp(-z));
		return f;
	}
}
