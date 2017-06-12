// Miriam Ruiz Ferrandez
// Juana L. Redondo
// Pilar M. Ortigosa
// Benjamin Ivorra
package org.uma.jmetal.operator.impl.localsearch;

import org.uma.jmetal.operator.LocalSearchOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.solution.util.RepairDoubleSolution;
import org.uma.jmetal.solution.util.RepairDoubleSolutionAtBounds;
import org.uma.jmetal.util.archive.impl.NonDominatedSolutionListArchive;
import org.uma.jmetal.util.comparator.DominanceComparator;
import org.uma.jmetal.util.pseudorandom.impl.ExtendedPseudoRandomGenerator;
import org.uma.jmetal.util.pseudorandom.impl.JavaRandomGenerator;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * This class implements a SASS local search operator for
 * multi-objective problems.
 *
 * @author Miriam R. Ferrandez <mrferrandez@ual.es>
 */
@SuppressWarnings("serial")
public class SASSLocalSearch<S extends Solution<?>> implements LocalSearchOperator<S> {
	private Problem<S> problem;
	private NonDominatedSolutionListArchive<S> externalList;
	private int improvementRounds ;
	//private Comparator<S> constraintComparator ;
	private Comparator<S> dominanceComparator ;

	//private MutationOperator<S> mutationOperator;
	private int evaluations ;

	private int numberOfImprovements ;
	private int numberOfNonComparableSolutions ;
	
	private RepairDoubleSolution solutionRepair ;
	
	private double radius;

	/**
	 * Constructor.
	 * Creates a new local search object.
	 * @param improvementRounds number of iterations
	 * @param problem problem to resolve
	 */
	public SASSLocalSearch(int improvementRounds, Problem<S> problem){
		this.problem=problem;
		//this.mutationOperator=mutationOperator;
		this.improvementRounds=improvementRounds;
		externalList=new NonDominatedSolutionListArchive<S>();
		dominanceComparator  = new DominanceComparator<S>();
		//constraintComparator = new OverallConstraintViolationComparator<S>();

		numberOfImprovements = 0 ;
		numberOfNonComparableSolutions = 0 ;
		
		solutionRepair = new RepairDoubleSolutionAtBounds();
	}

	/**
	 * Executes the local search.
	 * @param  solution The solution to improve
	 * @return An improved solution
	 */
	@SuppressWarnings("unchecked")
	public S execute(S solution) {
		int ic = 0;
		int best;
		boolean isnotdomexternalList;
		evaluations = 0;
		numberOfNonComparableSolutions = 0 ;
		
		int nvar = problem.getNumberOfVariables();
		int fcnt = 0; // number of failures
		int scnt = 0; // number of successes
		// b : normalized bias term to direct the search
		double[] b = new double[nvar]; 
		JavaRandomGenerator rnd = new JavaRandomGenerator();
		ExtendedPseudoRandomGenerator extpsrnd = new ExtendedPseudoRandomGenerator(rnd);
		
		int Maxfcnt = 5; // Maximum number of consecutive failures
		int Scnt = 5;
		int Fcnt = 3; 
		double ex= 2.0; // Expansion constant 
		double ct = 0.5; // Contraction constant

		int rounds = improvementRounds; // icmax
		
		DoubleSolution yold = (DoubleSolution) solution.copy();
		//problem.evaluate((S)yold);
		DoubleSolution ynew = (DoubleSolution) solution.copy();
		
		// sigmaub : normalized radius of the original individual
		double sigmaub = radius;//1 * Math.sqrt(nvar); //NORMALIZED radius
				// sqrt (sum (yold.getUpperBound(i)-yold.getLowerBound(i))^2)
		// sigma : size
		double sigma = sigmaub;
		// sigmalb : 
		double sigmalb = Double.max(sigmaub/1000.0, 1.0E-5);
		
		double[] zaux = new double[nvar];
		double gaussrnd, lowerBound, upperBound; //, ynorm, ydesnorm;

		while ((ic < rounds)&&(fcnt < Maxfcnt)) {
//			System.out.println("SASS Iteration: "+ic);
			if (scnt > Scnt){
				sigma = ex*sigma;
				////System.out.println(scnt+" mayor que "+Scnt);
				////System.out.println("Sigma= "+sigma);
			}
			if (fcnt > Fcnt){
					sigma = ct*sigma;
					////System.out.println(fcnt+" mayor que "+Fcnt);
					////System.out.println("Sigma= "+sigma);
			}
			if (sigma < sigmalb){
				sigma = sigmaub;
				for (int i = 0; i < nvar; i++){
					b[i] = 0.0;
				}
				////System.out.println("Sigma menor que "+sigmalb);
			}
			if (sigma > sigmaub){
				sigma = sigmaub;
				////System.out.println("Sigma mayor que "+sigmalb);
			}
			
			////System.out.println("Valor f obj: "+ yold.getObjective(0)+" "+yold.getObjective(1));
			// Generate a multivariate Gaussian random vector
			for (int i = 0; i < nvar; i++){
				gaussrnd = extpsrnd.randNormal(b[i], sigma);
				////System.out.println("Gaussian random vector "+gaussrnd);
			//}	
				if (gaussrnd > sigmaub){
					gaussrnd = sigmaub;
				}
				else if(gaussrnd < -sigmaub){
					gaussrnd = -sigmaub;
				}
				lowerBound = yold.getLowerBound(i);
				upperBound = yold.getUpperBound(i);
				zaux[i] = gaussrnd * (upperBound - lowerBound); //RAD NORMALIZADO
				////System.out.println("Gaussian random vector "+zaux[i]);
				//ynorm = (yold.getVariableValue(i) - lowerBound)/(upperBound - lowerBound);
				//ynorm = ynorm + gaussrnd;
				//ydesnorm = lowerBound + (upperBound - lowerBound)*ynorm;
				
				//System.out.println("yold "+yold.getVariableValue(i));
			//for (int i = 0; i < nvar; i++){
				ynew.setVariableValue(i, solutionRepair.repairSolutionVariableValue(
						yold.getVariableValue(i) + zaux[i], lowerBound, upperBound));
				//System.out.println("ynew sum"+ynew.getVariableValue(i));
			}
			problem.evaluate((S)ynew);
			////System.out.println("Nuevo valor f obj: "+ ynew.getObjective(0)+" "+ynew.getObjective(1));
			evaluations ++;
			best = dominanceComparator.compare((S)ynew, (S)yold);
			if (best == -1){ // ynew dominates yold
				scnt ++;
				fcnt = 0;
				////System.out.println("ynew dominates yold");
				yold = (DoubleSolution) ynew.copy();
			}
			else{
				isnotdomexternalList = externalList.add((S)ynew);
				if (isnotdomexternalList == true){ // It is also considered a success 
					scnt ++;  // = 0;
					fcnt = 0; // ++;
					numberOfNonComparableSolutions ++;
					for (int i = 0; i < nvar; i++){
						b[i] = 0.4*zaux[i] + 0.2*b[i];
					}
					////System.out.println("ynew added to nondominated external list");
				}
				else{
					for (int i = 0; i < nvar; i++){
						lowerBound = yold.getLowerBound(i);
						upperBound = yold.getUpperBound(i);
						ynew.setVariableValue(i, solutionRepair.repairSolutionVariableValue(
								yold.getVariableValue(i) - zaux[i], lowerBound, upperBound));
						//System.out.println("ynew resta"+ynew.getVariableValue(i));
					}
					problem.evaluate((S)ynew);
					////System.out.println("Nuevo valor resta f obj: "+ ynew.getObjective(0)+" "+ynew.getObjective(1));
					evaluations ++;
					best = dominanceComparator.compare((S)ynew, (S)yold);
					if (best == -1){ // ynew dominates yold
						scnt ++;
						fcnt = 0;
						////System.out.println("ynew dominates yold resta");
						yold = (DoubleSolution) ynew.copy();
					}
					else{
						isnotdomexternalList = externalList.add((S)ynew);
						if (isnotdomexternalList == true){ // It is also considered a success
							scnt ++;  // = 0;
							fcnt = 0; // ++;
							numberOfNonComparableSolutions ++;
							for (int i = 0; i < nvar; i++){
								b[i] = b[i] - 0.4*zaux[i];
							}
							////System.out.println("ynew added to nondominated external list resta");
						}
						else{
							for (int i = 0; i < nvar; i++){
								b[i] = 0.5*b[i];
							}
							fcnt ++;
							scnt = 0;
							////System.out.println("fcnt ++; fcnt = "+fcnt);
						}
					}
				}
			}
			ic ++;
			//yold = (DoubleSolution) ynew.copy();
//			System.out.println("Salgo de la ite");
		}
			return (S) yold.copy();
		}

		/**
		 * Returns the number of evaluations
		 */
		public int getEvaluations() {
			return evaluations;
		}

		@Override public int getNumberOfImprovements() {
			return numberOfImprovements ;
		}

		@Override public int getNumberOfNonComparableSolutions() {
			return numberOfNonComparableSolutions ;
		}
		
		//@SuppressWarnings("unchecked")
		public NonDominatedSolutionListArchive<S> getExternalList() {
			//List<S> copyExtList = new ArrayList<>();
			//for (int i=0; i<this.externalList.size(); i++){
			//	copyExtList.add((S) this.externalList.get(i).copy());
			//}
			return externalList;
		}
		
		public void clearExternalList() {
			externalList = new NonDominatedSolutionListArchive<S>();
		}
		public boolean addToExternalList(S solution) {
			return externalList.add(solution);
		}
		public void setRadius(double rad) {
			this.radius = rad;
		}
		
	}
