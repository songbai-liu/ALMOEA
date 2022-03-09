package jmetal.metaheuristics.almoea;

import jmetal.util.PseudoRandom;

/**
 * Utilities methods to used by MOEA/D
 */
public class Utils {

	public static double distVector(double[] vector1, double[] vector2) {
		int dim = vector1.length;
		double sum = 0;
		for (int n = 0; n < dim; n++) {
			sum += (vector1[n] - vector2[n]) * (vector1[n] - vector2[n]);
		}
		return Math.sqrt(sum);
	} // distVector

	public static void minFastSort(double x[], int idx[], int n, int m) {
		for (int i = 0; i < m; i++) {
			for (int j = i + 1; j < n; j++) {
				if (x[i] > x[j]) {
					double temp = x[i];
					x[i] = x[j];
					x[j] = temp;
					int id = idx[i];
					idx[i] = idx[j];
					idx[j] = id;
				} // if
			}
		} // for

	} // minFastSort
	
	public static void maxFastSort(double x[], int idx[], int n, int m) {
		for (int i = 0; i < m; i++) {
			for (int j = i + 1; j < n; j++) {
				if (x[i] < x[j]) {
					double temp = x[i];
					x[i] = x[j];
					x[j] = temp;
					int id = idx[i];
					idx[i] = idx[j];
					idx[j] = id;
				} // if
			}
		} // for

	} // minFastSort

	public static void randomPermutation(int[] perm, int size) {
		int[] index = new int[size];
		boolean[] flag = new boolean[size];

		for (int n = 0; n < size; n++) {
			index[n] = n;
			flag[n] = true;
		}

		int num = 0;
		while (num < size) {
			int start = PseudoRandom.randInt(0, size - 1);
			// int start = int(size*nd_uni(&rnd_uni_init));
			while (true) {
				if (flag[start]) {
					perm[num] = index[start];
					flag[start] = false;
					num++;
					break;
				}
				if (start == (size - 1)) {
					start = 0;
				} else {
					start++;
				}
			}
		} // while
	} // randomPermutation
	static void QuickSort(double[] array, int[] idx, int from, int to) {
		if (from < to) {
			double temp = array[to];
			int tempIdx = idx[to];
			int i = from - 1;
			for (int j = from; j < to; j++) {
				if (array[j] <= temp) {
					i++;
					double tempValue = array[j];
					array[j] = array[i];
					array[i] = tempValue;
					int tempIndex = idx[j];
					idx[j] = idx[i];
					idx[i] = tempIndex;
				}
			}
			array[to] = array[i + 1];
			array[i + 1] = temp;
			idx[to] = idx[i + 1];
			idx[i + 1] = tempIdx;
			QuickSort(array, idx, from, i);
			QuickSort(array, idx, i + 1, to);
		}
	}
	
	static void QuickSort(double[] array, int from, int to) {
		if (from < to) {
			double temp = array[to];
			int i = from - 1;
			for (int j = from; j < to; j++) {
				if (array[j] <= temp) {
					i++;
					double tempValue = array[j];
					array[j] = array[i];
					array[i] = tempValue;
				}
			}
			array[to] = array[i + 1];
			array[i + 1] = temp;
			QuickSort(array, from, i);
			QuickSort(array, i + 1, to);
		}
	}
}
