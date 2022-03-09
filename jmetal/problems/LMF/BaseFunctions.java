package jmetal.problems.LMF;

public class BaseFunctions {
	/*bf1: Sphere Function*/
	public static double getSphere(double x[]) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++)
			sum += x[i]*x[i];		
		return sum/x.length;
	}//Unimodal and Separable
	
	/*bf2: Schwefel 1.2 Function*/
	public static double getSchwefel1(double x[]) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			double mid = 0.0;
			for(int j = 0; j < i; j++) {
				mid += x[j];
			}
			sum += mid*mid;
		}			
		return sum/x.length;
	}//Unimodal and Non-Separable
	
	/*bf3: Schwefel 2.21 Function*/
	public static double getSchwefel2(double x[]) {
		double max = Math.abs(x[0]);
		for (int i = 1; i < x.length; i++) {
			if(max < Math.abs(x[i])) {
				max = Math.abs(x[i]);
			}
		}	
		return max/x.length;
	}//Unimodal and Non-Separable
	
	/*bf4: Ackley Function*/
	public static double getAckley(double x[]) {
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 0; i < x.length; i++) {
			sum1 += ((x[i] * x[i]) / x.length);
			sum2 += (Math.cos(2 * Math.PI * x[i]) / x.length);
		}
		return (-20 * Math.exp(-0.2 * Math.sqrt(sum1)) - Math.exp(sum2) + 20 + Math.E)/x.length;
	}//Multimodal and Separable	
	
	/*bf5: Rastrigin Function*/
	public static double getRastrigin(double x[]) {
		double result = 0.0;
		double a = 10.0;
		double w = 2 * Math.PI;
		for (int i = 0; i < x.length; i++) {
			result += x[i] * x[i] - a * Math.cos(w * x[i]);
		}
		result += a * x.length;
		return result/x.length;
	}//Multimodal and Separable	
	
	/*bf6: Rosenbrock Function*/
	public static double getRosenbrock(double x[]) {
		double sum = 0;
		for (int i = 0; i < x.length - 1; i++) {
			double t = 100 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) + (1 - x[i]) * (1 - x[i]);
			sum += t;
		}
		return sum/x.length;
	}//Multimodal and Non-Separable	

}
