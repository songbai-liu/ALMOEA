package jmetal.util.archive;

import java.util.Comparator;

import jmetal.core.Solution;
import jmetal.util.Spea2Fitness;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.EqualSolutions;
import jmetal.util.comparators.FitnessComparator;

public class SPEA2DensityArchive extends Archive{
	/**
	 * Stores the maximum size of the archive.
	 */
	private int maxSize_;

	/**
	 * Stores a <code>Comparator</code> for dominance checking.
	 */
	private Comparator dominance_;

	/**
	 * Stores a <code>Comparator</code> for fitness checking.
	 */
	private Comparator fitnessComparator_;

	/**
	 * Stores a <code>Comparator</code> for equality checking (in the objective
	 * space).
	 */
	private Comparator equals_;

	/**
	 * Constructor.
	 * 
	 * @param maxSize
	 *            The maximum size of the archive.
	 */
	public SPEA2DensityArchive(int maxSize) {
		super(maxSize);
		maxSize_ = maxSize;
		dominance_ = new DominanceComparator();
		equals_ = new EqualSolutions();
		fitnessComparator_ = new FitnessComparator();
	} // StrengthRawFitnessArchive

	/**
	 * Adds a <code>Solution</code> to the archive. If the <code>Solution</code>
	 * is dominated by any member of the archive then it is discarded. If the
	 * <code>Solution</code> dominates some members of the archive, these are
	 * removed. If the archive is full and the <code>Solution</code> has to be
	 * inserted, all the solutions are ordered by his strengthRawFitness value
	 * and the one having the worst value is removed.
	 * 
	 * @param solution
	 *            The <code>Solution</code>
	 * @return true if the <code>Solution</code> has been inserted, false
	 *         otherwise.
	 */
	public boolean add(Solution solution) {
		int flag = 0;
		int i = 0;
		Solution aux;
		while (i < solutionsList_.size()) {
			aux = solutionsList_.get(i);
			flag = dominance_.compare(solution, aux);
			if (flag == 1) { // The solution to add is dominated
				return false; // Discard the new solution
			} else if (flag == -1) { // A solution in the archive is dominated
				solutionsList_.remove(i); // Remove the dominated solution
			} else {
				if (equals_.compare(aux, solution) == 0) {
					return false;
				}
				i++;
			}
		}
		// Insert the solution in the archive
		solutionsList_.add(solution);

		if (size() > maxSize_) { // The archive is full
			(new Spea2Fitness(this)).fitnessAssign();
			// Remove the last
			remove(indexWorst(fitnessComparator_));
		}
		return true;
	} // add
}
