package org.uma.jmetal.runner.multiobjective;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.archive.impl.NonDominatedSolutionListArchive;
import org.uma.jmetal.util.comparator.CrowdingDistanceComparator;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.front.Front;
import org.uma.jmetal.util.front.imp.ArrayFront;
import org.uma.jmetal.util.front.util.FrontUtils;
import org.uma.jmetal.util.solutionattribute.impl.CrowdingDistance;

public class CrowdingRunner {
	
	@SuppressWarnings("unchecked")
	public static <S extends Solution<?>> void main(String[] args) throws FileNotFoundException{
		
		// CROWDING DISTANCE
		Front frontP = new ArrayFront(args[0]);
		Front frontEL = new ArrayFront(args[1]);
		//List<S> population = (List<S>) FrontUtils.convertFrontToSolutionList(frontP);
		List<S> population = SolutionListUtils.getNondominatedSolutions((List<S>) FrontUtils.convertFrontToSolutionList(frontP));
		List<S> externalList = (List<S>) FrontUtils.convertFrontToSolutionList(frontEL);
		NonDominatedSolutionListArchive<S> jointNDPop= new NonDominatedSolutionListArchive<S>();
		jointNDPop.initializeRemovedSolutions();
		//List<S> jointPopulation = new ArrayList<>();
		//jointPopulation.addAll(population);
		//jointPopulation.addAll(externalList);
		for (int j=0; j < population.size(); j++){
			jointNDPop.add(population.get(j));
		}
		for (int j=0; j < externalList.size(); j++){
			jointNDPop.add(externalList.get(j));
		}
		
		CrowdingDistance<S> crowdingDistance = new CrowdingDistance<S>() ;
		crowdingDistance.computeDensityEstimator(jointNDPop.getSolutionList());
		
		Collections.sort(jointNDPop.getSolutionList(), new CrowdingDistanceComparator<S>()) ;
		int solutionsToSelect = 30;
		solutionsToSelect = Integer.min(solutionsToSelect, jointNDPop.size());
		List<S> finalPopulation = new ArrayList<>(solutionsToSelect);
	    int i = 0 ;
	    while (finalPopulation.size() < solutionsToSelect) {
	    	finalPopulation.add(jointNDPop.get(i)) ;
	    	i++ ;
	    }
	    new SolutionListOutput(finalPopulation)
	    .setSeparator("\t")
	    .setVarFileOutputContext(new DefaultFileOutputContext("VARcrowdingND.tsv"))
	    .setFunFileOutputContext(new DefaultFileOutputContext("FUNcrowdingND.tsv"))
	    .print();
	}
	
}
