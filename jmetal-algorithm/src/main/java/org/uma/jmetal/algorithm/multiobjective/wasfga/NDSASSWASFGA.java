package org.uma.jmetal.algorithm.multiobjective.wasfga;

import org.uma.jmetal.algorithm.multiobjective.mombi.AbstractMOMBI;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.ASFWASFGA;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.AbstractUtilityFunctionsSet;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.Normalizer;
import org.uma.jmetal.algorithm.multiobjective.wasfga.util.WASFGARanking;
import org.uma.jmetal.algorithm.multiobjective.wasfga.util.WeightVector;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.localsearch.SASSLocalSearch;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.archive.impl.NonDominatedSolutionListArchive;
import org.uma.jmetal.util.comparator.EqualSolutionsComparator;
import org.uma.jmetal.util.distance.impl.EuclideanDistanceBetweenSolutionsInObjectiveSpace;
import org.uma.jmetal.util.distance.impl.EuclideanDistanceBetweenSolutionsInSolutionSpace;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.pseudorandom.impl.JavaRandomGenerator;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.DistanceToReferencePoint;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * Implementation of the preference based algorithm named WASF-GA on jMetal5.0
 *
 * @author Juanjo Durillo
 *
 *         This algorithm is described in the paper: A.B. Ruiz, R. Saborido, M.
 *         Luque "A Preference-based Evolutionary Algorithm for Multiobjective
 *         Optimization: The Weighting Achievement Scalarizing Function Genetic
 *         Algorithm". Journal of Global Optimization. May 2015, Volume 62,
 *         Issue 1, pp 101-129
 *         DOI = {10.1007/s10898-014-0214-y}
 *         
 *         Modified version of WASFGA called Non-Dominated WASFGA 
 * 		   Local Search Operator SASS for improving dominated individuals
 * @author Miriam R. Ferrandez 
 *         Juana L. Redondo
 *         
 */
public class NDSASSWASFGA<S extends Solution<?>> extends WASFGA<S> {
	/**
	 *
	 */
//  Quitamos todo lo que hereda de WASFGA
//	private static final long serialVersionUID = 1L;
//	protected int maxEvaluations;
//	protected int evaluations;
//	protected Normalizer normalizer;
	
//	final ASFWASFGA<S> achievementScalarizingFunction;
//	List<Double> referencePoint = null;
	
	//private NonDominatedSolutionListArchive<S> externalList;
	//private int numberOfNonComparableSolutions ;
	private EuclideanDistanceBetweenSolutionsInObjectiveSpace<S> distance ;
	private EuclideanDistanceBetweenSolutionsInSolutionSpace<Solution<Double>> distDomain;
	private JavaRandomGenerator rnd ;
	protected DistanceToReferencePoint<S> dist2refPoint;
	protected double maxDist2refPoint;

	protected SASSLocalSearch<S> localSearch;
	
	/**
	 * Constructor
	 *
	 * @param problem
	 *            Problem to solve
	 */
	public NDSASSWASFGA(Problem<S> problem,
								int populationSize,
								int maxIterations,
								SASSLocalSearch<S> localSearch,
								CrossoverOperator<S> crossoverOperator,
								MutationOperator<S> mutationOperator,
								SelectionOperator<List<S>, S> selectionOperator,
								SolutionListEvaluator<S> evaluator,
								List<Double> referencePoint) {

		super(problem,populationSize,maxIterations,crossoverOperator,mutationOperator,selectionOperator,evaluator,referencePoint);
		//setMaxPopulationSize();
		//this.referencePoint 				= referencePoint;
		//this.achievementScalarizingFunction =  createUtilityFunction();
		//externalList = new NonDominatedSolutionListArchive<S>();
		//numberOfNonComparableSolutions = 0 ;
		distance = new EuclideanDistanceBetweenSolutionsInObjectiveSpace<S>() ;
		distDomain = new EuclideanDistanceBetweenSolutionsInSolutionSpace<Solution<Double>>();
		rnd = new JavaRandomGenerator();
		dist2refPoint = new DistanceToReferencePoint<S>();
		this.localSearch = localSearch;
	}

	public ASFWASFGA<S> createUtilityFunction() {
		double [][] weights;
		if (referencePoint.size()==2)
			weights = WeightVector.initUniformWeights2D(0.005, getMaxPopulationSize());
		else{
			String filePath=System.getProperty("user.dir")+"/weightsWASFGA/W3D_"+Integer.toString(getMaxPopulationSize())+".dat";
			weights = WeightVector.getWeightsFromFile(filePath);
		}
		weights = WeightVector.invertWeights(weights,true);
		ASFWASFGA<S> aux = new ASFWASFGA<>(weights,referencePoint);

		return aux;
	}

	public int getPopulationSize() {
		return getMaxPopulationSize();
	}

	@Override
	public void specificMOEAComputations() {
		initializeBounds(this.getProblem().getNumberOfObjectives());
		updateNadirPoint(this.getPopulation());
		updateReferencePoint(this.getPopulation());
		localSearch.getExternalList().initializeRemovedSolutions();
		this.evaluations += this.getPopulationSize();
//		System.out.println("Number of evaluations: "+ this.evaluations);
//		System.out.println("External List Size: "+ localSearch.getExternalList().size());
//		new SolutionListOutput(this.getPopulation())
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_POP"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_POP"+Integer.toString(iterations)+".tsv"))
//	    .print();
		double value = 0.0;
		double maxval= 0.0;
		for (int i=0; i < this.getMaxPopulationSize(); i++){
			value = distance.getDistance(referencePoint, this.getPopulation().get(i));
			dist2refPoint.setAttribute(this.getPopulation().get(i), value);
			if (value > maxval){
				maxval = value;
			}
		}
		maxDist2refPoint = maxval;
	}
	
//	@Override
//	public void run() {
//		List<S> offspringPopulation;
//		List<S> matingPopulation;
//		
//
//		this.setPopulation(createInitialPopulation());
//		this.evaluatePopulation(this.getPopulation());
//		classify(this.getPopulation());
//		
//		initProgress();
//		//specific GA needed computations
//		this.specificMOEAComputations();
//		while (!isStoppingConditionReached()) {
//			matingPopulation = selection(this.getPopulation());
//			offspringPopulation = reproduction(matingPopulation);
//			offspringPopulation = evaluatePopulation(offspringPopulation);
//			this.setPopulation(replacement(this.getPopulation(), offspringPopulation));
//			updateProgress();
//			// specific GA needed computations
//			this.specificMOEAComputations();
//		}
//	}

	@SuppressWarnings("unchecked")
	@Override
	protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
		List<S> jointPopulation = new ArrayList<>();
		//jointPopulation.addAll(population); //NO AÑADO P porque ya se añaden EXTERNAL LIST e IMPROVED
		jointPopulation.addAll(offspringPopulation);
		
		List<S> improvedPopulation;
		List<S> dominatedPopulation = new ArrayList<>();
		boolean isnotdominated;
		for (int i=0; i < population.size(); i++){
        	isnotdominated = localSearch.addToExternalList((S)population.get(i).copy());
        	if (!isnotdominated){
        		dominatedPopulation.add(population.get(i));
    		}
		}
//		System.out.println("Number of dominated solutions: "+ dominatedPopulation.size());
		dominatedPopulation.addAll(localSearch.getExternalList().getRemovedSolutions());
		
		if (dominatedPopulation.size()>0){
			improvedPopulation = improveWithSASS(dominatedPopulation);
			jointPopulation.addAll(improvedPopulation);
		}
//		System.out.println("External List Size before fitting: "+ localSearch.getExternalList().size());
		fitExternalList();
		
//		new SolutionListOutput( localSearch.getExternalList().getSolutionList())
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_EL"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_EL"+Integer.toString(iterations)+".tsv"))
//	    .print();
		
		jointPopulation.addAll(localSearch.getExternalList().getSolutionList());
//		System.out.println("SASSexternalListPopulationSize After Fitting= "+ localSearch.getExternalList().size());
		
		Ranking<S> ranking = computeRanking(jointPopulation);

		return selectBest(ranking);
	}
	
//	@SuppressWarnings("unchecked")
//	protected void classify(List<S> population) {
//		boolean isnotdominated;
//		for (int i=0; i < population.size(); i++){
//        		isnotdominated = externalList.add((S)population.get(i).copy());
//        		if (!isnotdominated){
////        			numberOfNonComparableSolutions ++;
////        		}
////        		else{
//        			dominatedPopulation.add(population.get(i));
//        		}
//		}
//		System.out.println("Number of non comparable sol: "+ numberOfNonComparableSolutions);
//	}
	
	protected List<S> improveWithSASS(List<S> dominatedPopulation) {
//		new SolutionListOutput(dominatedPopulation)
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_DOM"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_DOM"+Integer.toString(iterations)+".tsv"))
//	    .print();
		List<S> improvedPopulation = new ArrayList<>(dominatedPopulation.size());
		double rad;
		for (int i = 0; i < dominatedPopulation.size(); i++) {
			S individual = dominatedPopulation.get(i);
			double highestDistance = 0.0;
			int index = 0;
			// Calculate radius
			for (int j = 0; j< this.getMaxPopulationSize();j++){
			      double aux = distance.getDistance(individual, this.getPopulation().get(j));
			      if (aux > highestDistance){
			    	  highestDistance = aux;
			    	  index = j;
			      }
			}
			rad=distDomain.getNormalizedDistance((DoubleSolution)individual, (DoubleSolution) this.getPopulation().get(index));
			localSearch.setRadius(rad);
//			System.out.println("RADIUS: " + rad);
			improvedPopulation.add(localSearch.execute(individual));
			this.evaluations += localSearch.getEvaluations();
		}
		//this.setPopulation(SASSPopulation);
		
		//List<S> externalListSASS = localSearch.getExternalList();
		//System.out.println("externalListPopulationSize= "+ externalListPopulation.size());
		//new SolutionListOutput(externalListPopulation)
	    //.setSeparator("\t")
	    //.setVarFileOutputContext(new DefaultFileOutputContext("VAR"+Integer.toString(iterations)+".tsv"))
	    //.setFunFileOutputContext(new DefaultFileOutputContext("FUN"+Integer.toString(iterations)+".tsv"))
	    //.print();
//		System.out.println("Number of improved solutions: "+ improvedPopulation.size());
//
//		new SolutionListOutput(improvedPopulation)
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_IMPR"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_IMPR"+Integer.toString(iterations)+".tsv"))
//	    .print();
		
		return improvedPopulation;
	}
	
	@SuppressWarnings("unchecked")
	protected S improveDominatedSolutions(S dominatedSolution) {
		DoubleSolution improvedSolution = (DoubleSolution) dominatedSolution.copy();
		double bestDistance = Double.MAX_VALUE;
		int index=0;
		double aux;
		for (int i=0; i < localSearch.getExternalList().size(); i++){
			aux = distance.getDistance(dominatedSolution, localSearch.getExternalList().get(i));
			if (aux < bestDistance){
		        bestDistance = aux;
		        index = i;
			}
		}
		double x0,x1,x2, lambda;
		for (int i=0; i < dominatedSolution.getNumberOfVariables(); i++){
			x0 = (double) localSearch.getExternalList().get(index).getVariableValue(i);
			x1 = (double) dominatedSolution.getVariableValue(i);
			lambda = rnd.nextDouble();
			x2 = x1 + lambda * (x0-x1);
			improvedSolution.setVariableValue(i, x2);
		}
		
		return (S) improvedSolution;
	}
	
	protected List<S> improveDominatedSolutions(List<S> dominatedPopulation) {
//		new SolutionListOutput(dominatedPopulation)
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_DOM"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_DOM"+Integer.toString(iterations)+".tsv"))
//	    .print();
//		boolean isnotdominated;
		List<S> improvedPopulation = new ArrayList<>();
		for (int i=0; i < dominatedPopulation.size(); i++){
			improvedPopulation.add(improveDominatedSolutions(dominatedPopulation.get(i)));
		}
		improvedPopulation = evaluatePopulation(improvedPopulation);
		this.evaluations += improvedPopulation.size();
//		for (int i=0; i < improvedPopulation.size(); i++){
//			isnotdominated = externalList.add((S)improvedPopulation.get(i));
//			if (isnotdominated)
//				numberOfNonComparableSolutions ++;
//			else{
//				dominatedPopulation.add(improvedPopulation.get(i));
//			}
//		}
		
//		System.out.println("Number of improved solutions: "+ improvedPopulation.size());
//
//		new SolutionListOutput(improvedPopulation)
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_IMPR"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_IMPR"+Integer.toString(iterations)+".tsv"))
//	    .print();
		return improvedPopulation;
	}
	
	protected void fitExternalList(){
		if (localSearch.getExternalList().size() > this.getPopulationSize()){
			double value;
			//for (int i=0; i < localSearch.getExternalList().size(); i++){
				
			//	dist2refPoint.setAttribute(localSearch.getExternalList().get(i), value);
			//}
			Iterator<S> iterator = localSearch.getExternalList().getSolutionList().iterator();
			while (iterator.hasNext()) {
		        S listIndividual = iterator.next();
		        value = distance.getDistance(referencePoint, listIndividual);
		        //dist2refPoint.getAttribute(listIndividual)
		        if ( value > maxDist2refPoint){
		        	iterator.remove();
		        	System.out.println(value);
		        	System.out.println(maxDist2refPoint);
		        }
			}
		}
	}
	
	protected Ranking<S> computeRanking(List<S> solutionList) {
		this.achievementScalarizingFunction.setNadir(getNadirPoint(solutionList));
		this.achievementScalarizingFunction.setUtopia(getReferencePoint(solutionList));
		Ranking<S> ranking = new WASFGARanking<>(this.achievementScalarizingFunction);
		ranking.computeRanking(solutionList);
		return ranking;
	}
	
	protected void addRankedSolutionsToPopulation(Ranking<S> ranking, int index, List<S> population) {
		population.addAll(ranking.getSubfront(index));
	}
	
	protected void addLastRankedSolutionsToPopulation(Ranking<S> ranking,int index, List<S>population) {
		List<S> front 	= ranking.getSubfront(index);
		int remain 		= this.getPopulationSize() - population.size();
		population.addAll(front.subList(0, remain));
	}
	
	protected List<S> selectBest(Ranking<S> ranking) {
		List<S> population = new ArrayList<>(this.getPopulationSize());
		int rankingIndex = 0;

		while (populationIsNotFull(population)) {
			if (subfrontFillsIntoThePopulation(ranking, rankingIndex, population)) {
				addRankedSolutionsToPopulation(ranking, rankingIndex, population);
				rankingIndex++;
			} else {
				addLastRankedSolutionsToPopulation(ranking, rankingIndex, population);
			}
		}
		return population;
	}

	private boolean subfrontFillsIntoThePopulation(Ranking<S> ranking, int index, List<S> population) {
		return (population.size()+ranking.getSubfront(index).size() < this.getPopulationSize());
	}
	protected AbstractUtilityFunctionsSet<S> getUtilityFunctions() {
		return this.achievementScalarizingFunction;
	}
	
	@Override public List<S> getResult() {
		return getNonDominatedSolutions(getPopulation());
	}
	protected List<S> getNonDominatedSolutions(List<S> solutionList) {
		return SolutionListUtils.getNondominatedSolutions(solutionList);
	}
	
	public List<S> getAllPopulation() {
		return getPopulation();
	}

	@Override public String getName() {
		return "NDSASSWASFGA" ;
	}

	@Override public String getDescription() {
		return "Non Dominated with SASS Weighting Achievement Scalarizing Function Genetic Algorithm" ;
	}
}
