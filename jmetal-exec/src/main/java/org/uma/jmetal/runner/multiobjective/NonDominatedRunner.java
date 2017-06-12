package org.uma.jmetal.runner.multiobjective;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.comparator.CrowdingDistanceComparator;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.front.Front;
import org.uma.jmetal.util.front.imp.ArrayFront;
import org.uma.jmetal.util.front.util.FrontUtils;
import org.uma.jmetal.util.solutionattribute.impl.CrowdingDistance;

public class NonDominatedRunner {
	
	@SuppressWarnings("unchecked")
	public static <S extends Solution<?>> void main(String[] args) throws FileNotFoundException{
		
		//---------------------------------------------------------------------------------
		// CALCULO DEL NUMERO DE SOLUCIONES NO DOMINADAS EN EL FRENTE DADO COMO ENTRADA
		Front front = new ArrayFront(args[0]);
		List<S> solutionList = (List<S>) FrontUtils.convertFrontToSolutionList(front);
		List<S> nondominatedList = SolutionListUtils.getNondominatedSolutions((List<S>) solutionList);
		System.out.println("Number of non dominated solutions: "+nondominatedList.size());
		//---------------------------------------------------------------------------------
	}
	
}
