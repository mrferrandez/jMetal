package org.uma.jmetal.algorithm.multiobjective.wasfga.util;

import org.uma.jmetal.algorithm.multiobjective.mombi.util.AbstractUtilityFunctionsSet;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.GenericSolutionAttribute;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


@SuppressWarnings("serial")
public class NewWASFGARanking<S extends Solution<?>> extends GenericSolutionAttribute<S, Integer> 
		implements Ranking<S> {

  private AbstractUtilityFunctionsSet<S> utilityFunctions;
  private List<List<S>> rankedSubpopulations;
  private int numberOfRanks = 0;

  public NewWASFGARanking(AbstractUtilityFunctionsSet<S> utilityFunctions) {
    this.utilityFunctions = utilityFunctions;
  }

  @Override
  public Ranking<S> computeRanking(List<S> population) {

	this.numberOfRanks 		= (population.size() + 1) / this.utilityFunctions.getSize();	
	this.rankedSubpopulations = new ArrayList<>(this.numberOfRanks);
	for (int i = 0; i < this.numberOfRanks; i++) {
		this.rankedSubpopulations.add(new ArrayList<S>());
	}
	List<S> temporalList 	= new LinkedList<>();
	temporalList.addAll(population);
	
	List<Integer> temporalWeights = new LinkedList<>();
	
	Double [][] matrix = new Double[this.utilityFunctions.getSize()][population.size()];
	List <IndexSorter<Double>> is = new LinkedList<>();//IndexSorter<Double> [this.utilityFunctions.getSize()];
	Integer[][] indexes = new Integer [this.utilityFunctions.getSize()][population.size()];
	double [] diff = new double [this.utilityFunctions.getSize()];  
	double maxdiff = 0.0;
	int maxdiffweight = 0;
	
	for (int weight = 0; weight < this.utilityFunctions.getSize(); weight++) {
		temporalWeights.add(weight);
		for (int solutionIdx = 0; solutionIdx < temporalList.size(); solutionIdx++) {
			matrix[weight][solutionIdx] = this.utilityFunctions.evaluate(temporalList.get(solutionIdx), weight);
		}
		//Arrays.sort(matrix[weight]);
		IndexSorter<Double> isweight = new IndexSorter<Double>(matrix[weight]);
		isweight.sort();
		is.add(isweight);
		// Observación: Los elementos de matrix no cambian su orden
		// Note that the array itself is not sorted, but can be accessed in a sorted way 
		// via the values of the getIndexes() method of the IndexSorter class.
//		System.out.print("Unsorted: ");
//		for ( Double d : matrix[weight] ){
//			System.out.print(d);
//			System.out.print("\t");
//		}
//		System.out.println();
//		System.out.print("Sorted");
//		for ( Integer i : isweight.getIndexes() ){
//			System.out.print(matrix[weight][i]);
//			System.out.print("\t");
//		}
//		System.out.println();
		indexes[weight] = isweight.getIndexes();
		// Calculate difference between the best and the second best individuals
		diff[weight] = matrix[weight][indexes[weight][1]] - matrix[weight][indexes[weight][0]];
		if (diff[weight] >= maxdiff){
			maxdiff = diff[weight];
			maxdiffweight = weight;
		}
	}
	int toRemoveIdx = indexes[maxdiffweight][0];
	temporalWeights.remove(maxdiffweight);
	S solutionToInsert = temporalList.remove(toRemoveIdx);//temporalList.get(toRemoveIdx);//
	setAttribute(solutionToInsert, 0); //idx = 0
	this.rankedSubpopulations.get(0).add(solutionToInsert); //idx = 0
	for (int weight = 0; weight < this.utilityFunctions.getSize(); weight++) {
		is.get(weight).remove(toRemoveIdx);
	}
	
	for (int idx = 0; idx < this.numberOfRanks; idx++) {
		int remainingWeights = temporalWeights.size();
		for (int weightRank = 0; weightRank < remainingWeights; weightRank++) {
			maxdiff = 0.0;
			maxdiffweight = 0;
			for (int weightIdx = 0; weightIdx < temporalWeights.size(); weightIdx++) {
				int weight = temporalWeights.get(weightIdx);
				//is.get(weight).remove(toRemoveIdx);
				indexes[weight] = is.get(weight).getIndexes();
				diff[weight] = matrix[weight][indexes[weight][1]] - matrix[weight][indexes[weight][0]];
				if (diff[weight] >= maxdiff){
					maxdiff = diff[weight];
					maxdiffweight = weightIdx; // He cambiado a weightIdx // Antes ponia weight 
				}
			}
			toRemoveIdx = indexes[temporalWeights.get(maxdiffweight)][0];
			temporalWeights.remove(maxdiffweight);
			solutionToInsert = temporalList.remove(toRemoveIdx);//temporalList.get(toRemoveIdx);// 
			setAttribute(solutionToInsert, idx);
			this.rankedSubpopulations.get(idx).add(solutionToInsert);
			for (int weight = 0; weight < this.utilityFunctions.getSize(); weight++) {
				is.get(weight).remove(toRemoveIdx);
			}
		}
		for (int weight = 0; weight < this.utilityFunctions.getSize(); weight++) {
			temporalWeights.add(weight);
		}
	}
	
	return this;
  }

  @Override
  public List<S> getSubfront(int rank) {
    return this.rankedSubpopulations.get(rank);
  }

  @Override
  public int getNumberOfSubfronts() {
    return this.rankedSubpopulations.size();
  }
  
  public AbstractUtilityFunctionsSet<S> getUtilityFunctions() {
    return this.utilityFunctions;
  }
  
  
  
}
