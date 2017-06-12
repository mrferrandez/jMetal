package org.uma.jmetal.algorithm.multiobjective.wasfga.util;

import org.uma.jmetal.algorithm.multiobjective.mombi.util.AbstractUtilityFunctionsSet;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.archive.impl.NonDominatedSolutionListArchive;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.GenericSolutionAttribute;
import org.uma.jmetal.util.solutionattribute.impl.OverallConstraintViolation;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

@SuppressWarnings("serial")
public class WASFGARankingND<S extends Solution<?>> extends GenericSolutionAttribute<S, Integer> 
		implements Ranking<S> {

  private AbstractUtilityFunctionsSet<S> utilityFunctions;
  private List<List<S>> rankedSubpopulations;
  private int numberOfRanks = 0;
  
  private NonDominatedSolutionListArchive<S> nondomList = new NonDominatedSolutionListArchive<S>();
  
  public WASFGARankingND(AbstractUtilityFunctionsSet<S> utilityFunctions) {
    this.utilityFunctions = utilityFunctions;
  }

  @SuppressWarnings("unchecked")
@Override
  public Ranking<S> computeRanking(List<S> population) {

	this.numberOfRanks 		= (population.size() + 1) / this.utilityFunctions.getSize();	
	this.rankedSubpopulations = new ArrayList<>(this.numberOfRanks);
	for (int i = 0; i < this.numberOfRanks; i++) {
		this.rankedSubpopulations.add(new ArrayList<S>());
	}
	List<S> temporalList 	= new LinkedList<>();
	//List<S> constrainedList 	= new LinkedList<>();
	List<S> dominatedList = new ArrayList<>();
	boolean isnotdominated;
	for (int i=0; i < population.size(); i++){
    	isnotdominated = nondomList.add((S)population.get(i).copy());
    	if (isnotdominated){
    		temporalList.add(population.get(i));
		}
    	else{
    		dominatedList.add(population.get(i));
    	}
	}
	
	//rankingNonDominated:
	for (int idx = 0; idx < this.numberOfRanks; idx++) {
		for (int weigth = 0; weigth < this.utilityFunctions.getSize(); weigth++) {
			if (temporalList.size()!=0){
				int toRemoveIdx = 0;
				double minimumValue = this.utilityFunctions.evaluate(temporalList.get(0), weigth);	
				for (int solutionIdx = 1; solutionIdx < temporalList.size(); solutionIdx++) {
					double value = this.utilityFunctions.evaluate(temporalList.get(solutionIdx), weigth);
				
					if (value < minimumValue) {
						minimumValue = value;
						toRemoveIdx = solutionIdx;
					}
				}
				S solutionToInsert = temporalList.remove(toRemoveIdx);
				setAttribute(solutionToInsert, idx);
				this.rankedSubpopulations.get(idx).add(solutionToInsert);
			}
			else{
				int toRemoveIdx = 0;
				double minimumValue = this.utilityFunctions.evaluate(dominatedList.get(0), weigth);	
				for (int solutionIdx = 1; solutionIdx < dominatedList.size(); solutionIdx++) {
					double value = this.utilityFunctions.evaluate(dominatedList.get(solutionIdx), weigth);
				
					if (value < minimumValue) {
						minimumValue = value;
						toRemoveIdx = solutionIdx;
					}
				}
				S solutionToInsert = dominatedList.remove(toRemoveIdx);
				setAttribute(solutionToInsert, idx);
				this.rankedSubpopulations.get(idx).add(solutionToInsert);
			}
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
