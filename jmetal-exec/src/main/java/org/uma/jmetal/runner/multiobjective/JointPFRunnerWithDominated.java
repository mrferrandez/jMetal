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

public class JointPFRunnerWithDominated {
	
	@SuppressWarnings("unchecked")
	public static <S extends Solution<?>> void main(String[] args) throws FileNotFoundException{
		
		List<S> jointPFdom= new ArrayList<S>();
		//jointPF.initializeRemovedSolutions();
		
		// READ ALL THE PARETO FRONTS
		// except the last four fronts (NSGAII, MOEAD, SMSEMOA, SPEA2)
		for (int i=0; i<(args.length-4); i++){
			Front front = new ArrayFront(args[i]);
			List<S> population = (List<S>) FrontUtils.convertFrontToSolutionList(front);
			for (int j=0; j < population.size(); j++){
				jointPFdom.add(population.get(j));
			}
		}
		int totalnumobj = jointPFdom.get(0).getNumberOfObjectives();
		double[] max = new double [totalnumobj];
		double[] min = new double [totalnumobj];
		for (int nobj=0; nobj<totalnumobj; nobj++){
			max[nobj] = jointPFdom.get(0).getObjective(nobj);
			min[nobj] = jointPFdom.get(0).getObjective(nobj);
		}
		double aux;
		for (int j=1; j<jointPFdom.size(); j++){
			for (int nobj=0; nobj<totalnumobj; nobj++){
				aux = jointPFdom.get(j).getObjective(nobj);
				if (aux < min[nobj]){
					min[nobj] = aux;
				} else if (aux > max[nobj]){
					max[nobj] = aux;
				}
			}
		}
		for (int i=args.length-4; i<args.length; i++){
			Front front = new ArrayFront(args[i]);
			List<S> population = (List<S>) FrontUtils.convertFrontToSolutionList(front);
			for (int j=0; j < population.size(); j++){
				boolean isIntoROI = true;
				for (int nobj=0; nobj<totalnumobj; nobj++){
					if ((population.get(j).getObjective(nobj)>max[nobj])||(population.get(j).getObjective(nobj)<min[nobj])){
						isIntoROI = false;
					}
				}
				if (isIntoROI){
					jointPFdom.add(population.get(j));
				}
			}
		}
		
	    new SolutionListOutput(jointPFdom)
	    .setSeparator("\t")
	    .setVarFileOutputContext(new DefaultFileOutputContext("VARjointPFdom.tsv"))
	    .setFunFileOutputContext(new DefaultFileOutputContext("FUNjointPFdom.tsv"))
	    .print();
	}
	
}
