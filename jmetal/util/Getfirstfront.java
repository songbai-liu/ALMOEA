//Ranking.java
//
//Author:
//   Antonio J. Nebro <antonio@lcc.uma.es>
//   Juan J. Durillo <durillo@lcc.uma.es>
//
//Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU Lesser General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU Lesser General Public License for more details.
//
//You should have received a copy of the GNU Lesser General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.util;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.OverallConstraintViolationComparator;

import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
* This class implements some facilities for ranking solutions. Given a
* <code>SolutionSet</code> object, their solutions are ranked according to
* scheme proposed in NSGA-II; as a result, a set of subsets are obtained. The
* subsets are numbered starting from 0 (in NSGA-II, the numbering starts from
* 1); thus, subset 0 contains the non-dominated solutions, subset 1 contains
* the non-dominated solutions after removing those belonging to subset 0, and
* so on.
*/
public class Getfirstfront {

/**
 * The <code>SolutionSet</code> to rank
 */
//private SolutionSet solutionSet_;

/**
 * An array containing all the fronts found during the search
 */
private SolutionSet ranking_;

/**
 * stores a <code>Comparator</code> for dominance checking
 */
//private static final Comparator dominance_ = new DominanceComparator();

/**
 * stores a <code>Comparator</code> for Overal Constraint Violation
 * Comparator checking
 */
//private static final Comparator constraint_ = new OverallConstraintViolationComparator();

/**
 * Constructor.
 * 
 * @param solutionSet
 *            The <code>SolutionSet</code> to be ranked.
 */
public Getfirstfront(SolutionSet solutionSet) {
	//solutionSet_ = solutionSet;

	// dominateMe[i] contains the number of solutions dominating i
	int[] dominateMe = new int[solutionSet.size()];
	
	Solution a1,a2; 
	
	for(int i=0;i<solutionSet.size();i++){
		dominateMe[i]=0;
	}

	// iDominate[k] contains the list of solutions dominated by k
	for(int i=0;i<solutionSet.size();i++){
		a1=solutionSet.get(i);
		for(int j=i+1;j<solutionSet.size();j++){
			if(dominateMe[i]==0){
				a2=solutionSet.get(j);
				int num1=0,num2=0,num3=0;
				for(int k=0;k<a1.numberOfObjectives();k++){
					if(a1.getObjective(k)>a2.getObjective(k)){
						num1++;
					}
					else if(a1.getObjective(k)==a2.getObjective(k)){
						num2++;
					}
					else
						num3++;
				}
				if(num1==0){
					dominateMe[j]=1;
				}
				else if (num3==0){
					dominateMe[i]=1;
				}
			}else{
				break;
			}
		}
	}
	
	ranking_=new SolutionSet(solutionSet.size());
	for(int i=0;i<solutionSet.size();i++){
		if(dominateMe[i]==0)
			ranking_.add(solutionSet.get(i));
	}
	//List<Integer>[] iDominate = new List[solutionSet_.size()];

	// front[i] contains the list of individuals belonging to the front i
	//List<Integer>[] front = new List[solutionSet_.size() + 1];

	// flagDominate is an auxiliar encodings.variable
	//int flagDominate;

	// Initialize the fronts
	//for (int i = 0; i < front.length; i++)
	//	front[i] = new LinkedList<Integer>();

	/*
	 * //-> Fast non dominated sorting algorithm for (int p = 0; p <
	 * solutionSet_.size(); p++) { // Initialice the list of individuals
	 * that i dominate and the number // of individuals that dominate me
	 * iDominate[p] = new LinkedList<Integer>(); dominateMe[p] = 0; // For
	 * all q individuals , calculate if p dominates q or vice versa for (int
	 * q = 0; q < solutionSet_.size(); q++) { flagDominate
	 * =constraint_.compare(solutionSet.get(p),solutionSet.get(q)); if
	 * (flagDominate == 0) { flagDominate
	 * =dominance_.compare(solutionSet.get(p),solutionSet.get(q)); }
	 * 
	 * if (flagDominate == -1) { iDominate[p].add(new Integer(q)); } else if
	 * (flagDominate == 1) { dominateMe[p]++; } }
	 * 
	 * // If nobody dominates p, p belongs to the first front if
	 * (dominateMe[p] == 0) { front[0].add(new Integer(p));
	 * solutionSet.get(p).setRank(0); } }
	 */

	// -> Fast non dominated sorting algorithm

} // Ranking

/**
 * Returns a <code>SolutionSet</code> containing the solutions of a given
 * rank.
 * 
 * @param rank
 *            The rank
 * @return Object representing the <code>SolutionSet</code>.
 */
public SolutionSet getSubfront() {
	return ranking_;
} // getSubFront

/**
 * Returns the total number of subFronts founds.
 */
} // Ranking
