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
 *         Modified version of WASFGA called  WASFGA advanced population
 * @author Miriam R. Ferrandez 
 *         Juana L. Redondo
 *         
 */
public class SASSWASFGA_ad<S extends Solution<?>> extends AbstractMOMBI<S> {
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	protected int maxEvaluations;
	protected int evaluations;
	protected Normalizer normalizer;
	
	final ASFWASFGA<S> achievementScalarizingFunction;
	List<Double> referencePoint = null;
	
	//private NonDominatedSolutionListArchive<S> externalList;
	//private int numberOfNonComparableSolutions ;
	private EuclideanDistanceBetweenSolutionsInObjectiveSpace<S> distance ;
	private EuclideanDistanceBetweenSolutionsInSolutionSpace<Solution<Double>> distDomain;
	private JavaRandomGenerator rnd ;
	protected DistanceToReferencePoint<S> dist2refPoint;
	protected double maxDist2refPoint;
	protected int bestIndex2refPoint;
	
	private List<S> oldPopulation ;
	
	protected SASSLocalSearch<S> localSearch;

	/**
	 * Constructor
	 *
	 * @param problem
	 *            Problem to solve
	 */
	public SASSWASFGA_ad(Problem<S> problem,
								int populationSize,
								int maxIterations,
								SASSLocalSearch<S> localSearch,
								CrossoverOperator<S> crossoverOperator,
								MutationOperator<S> mutationOperator,
								SelectionOperator<List<S>, S> selectionOperator,
								SolutionListEvaluator<S> evaluator,
								List<Double> referencePoint) {

		super(problem,maxIterations,crossoverOperator,mutationOperator,selectionOperator,evaluator);
		setMaxPopulationSize(populationSize);
		this.referencePoint 				= referencePoint;
		this.achievementScalarizingFunction =  createUtilityFunction();
		//externalList = new NonDominatedSolutionListArchive<S>();
		//numberOfNonComparableSolutions = 0 ;
		distance = new EuclideanDistanceBetweenSolutionsInObjectiveSpace<S>() ;
		distDomain = new EuclideanDistanceBetweenSolutionsInSolutionSpace<Solution<Double>>();
		rnd = new JavaRandomGenerator();
		dist2refPoint = new DistanceToReferencePoint<S>();
		this.evaluations = 0;
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
		//this.evaluations += this.getPopulationSize();
		System.out.println("Number of evaluations: "+ this.evaluations);
		//System.out.println("External List Size: "+ this.externalList.size());
		new SolutionListOutput(this.getPopulation())
	    .setSeparator("\t")
	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_POP"+Integer.toString(iterations)+".tsv"))
	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_POP"+Integer.toString(iterations)+".tsv"))
	    .print();
		
		// First of all: select the best individual in P, i. e., the one with the minimum distance
		// to the reference point
		int indexBest = 0;
		double smallestDistance = Double.MAX_VALUE;
		
		double value = 0.0;
		double maxval= 0.0;
		for (int i=0; i < this.getMaxPopulationSize(); i++){
			value = distance.getDistance(referencePoint, this.getPopulation().get(i));
			dist2refPoint.setAttribute(this.getPopulation().get(i), value);
			if (value > maxval){
				maxval = value;
			}
			if (value < smallestDistance){
				smallestDistance = value;
				indexBest = i;
			}
		}
		bestIndex2refPoint = indexBest;
		maxDist2refPoint = maxval;
	}
	
	@Override
	public void run() {
		List<S> offspringPopulation = new ArrayList<>();
		List<S> matingPopulation;

		this.setPopulation(createInitialPopulation());
		this.evaluatePopulation(this.getPopulation());
		this.evaluations += this.getPopulationSize();
		initProgress();
		//specific GA needed computations
		this.specificMOEAComputations();
		while (!isStoppingConditionReached()) {
			if (iterations==1){
				matingPopulation = selection(this.getPopulation());
				offspringPopulation = reproduction(matingPopulation);
				offspringPopulation = evaluatePopulation(offspringPopulation);
				this.evaluations += offspringPopulation.size();
			}
			
			this.setPopulation(replacement(this.getPopulation(), offspringPopulation));
			updateProgress();
			// specific GA needed computations
			this.specificMOEAComputations();
		}
	}

	@SuppressWarnings("unchecked")
	@Override
	protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
		List<S> jointPopulation = new ArrayList<>();
		
		List<S> advancedPopulation;
		if(iterations==1){
			jointPopulation.addAll(population);
			jointPopulation.addAll(offspringPopulation);
		}
		else{
			advancedPopulation = calculateAdvancedSolutions();
			jointPopulation.addAll(advancedPopulation);
		}
		
		oldPopulation = new ArrayList<>();
		List<S> dominatedPopulation = new ArrayList<>();
		boolean isnotdominated;
		for (int i=0; i < population.size(); i++){
			oldPopulation.add((S)population.get(i).copy());
			isnotdominated = localSearch.addToExternalList((S)population.get(i).copy());
        	if (!isnotdominated){
        		dominatedPopulation.add(population.get(i));
    		}
		}
		dominatedPopulation.addAll(localSearch.getExternalList().getRemovedSolutions());
		System.out.println("Number of dominated solutions: "+ dominatedPopulation.size());
		
		if (iterations!=1){
			if (dominatedPopulation.size()>0){
				jointPopulation.addAll(dominatedPopulation);
			}
			S improvedIndividual = improveWithSASS();
			jointPopulation.add(improvedIndividual);
		
			System.out.println("External List Size before fitting: "+ localSearch.getExternalList().size());
			fitExternalList();
		
			new SolutionListOutput( localSearch.getExternalList().getSolutionList())
			.setSeparator("\t")
			.setVarFileOutputContext(new DefaultFileOutputContext("VAR_EL"+Integer.toString(iterations)+".tsv"))
			.setFunFileOutputContext(new DefaultFileOutputContext("FUN_EL"+Integer.toString(iterations)+".tsv"))
			.print();
		
			jointPopulation.addAll(localSearch.getExternalList().getSolutionList());
			System.out.println("SASSexternalListPopulationSize After Fitting= "+ localSearch.getExternalList().size());	
		}
		
		
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
	
//	@SuppressWarnings("unchecked")
//	protected S improveDominatedSolutions(S dominatedSolution) {
//		DoubleSolution improvedSolution = (DoubleSolution) dominatedSolution.copy();
//		double bestDistance = Double.MAX_VALUE;
//		int index=0;
//		double aux;
//		for (int i=0; i < externalList.size(); i++){
//			aux = distance.getDistance(dominatedSolution, externalList.get(i));
//			if (aux < bestDistance){
//		        bestDistance = aux;
//		        index = i;
//			}
//		}
//		double x0,x1,x2, lambda;
//		for (int i=0; i < dominatedSolution.getNumberOfVariables(); i++){
//			x0 = (double) externalList.get(index).getVariableValue(i);
//			x1 = (double) dominatedSolution.getVariableValue(i);
//			lambda = rnd.nextDouble();
//			x2 = x1 + lambda * (x0-x1);
//			improvedSolution.setVariableValue(i, x2);
//		}
//		
//		return (S) improvedSolution;
//	}
	
//	protected List<S> improveDominatedSolutions(List<S> dominatedPopulation) {
//		new SolutionListOutput(dominatedPopulation)
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_DOM"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_DOM"+Integer.toString(iterations)+".tsv"))
//	    .print();
////		boolean isnotdominated;
//		List<S> improvedPopulation = new ArrayList<>();
//		for (int i=0; i < dominatedPopulation.size(); i++){
//			improvedPopulation.add(improveDominatedSolutions(dominatedPopulation.get(i)));
//		}
//		improvedPopulation = evaluatePopulation(improvedPopulation);
//		this.evaluations += improvedPopulation.size();
////		for (int i=0; i < improvedPopulation.size(); i++){
////			isnotdominated = externalList.add((S)improvedPopulation.get(i));
////			if (isnotdominated)
////				numberOfNonComparableSolutions ++;
////			else{
////				dominatedPopulation.add(improvedPopulation.get(i));
////			}
////		}
//		System.out.println("Number of improved solutions: "+ improvedPopulation.size());
//
//		new SolutionListOutput(improvedPopulation)
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_IMPR"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_IMPR"+Integer.toString(iterations)+".tsv"))
//	    .print();
//		return improvedPopulation;
//	}
	
	@SuppressWarnings("unchecked")
	protected S calculateAdvancedSolutions(S oldSolution, S newSolution) {
		DoubleSolution advancedSolution = (DoubleSolution) oldSolution.copy();
		//double bestDistance = Double.MAX_VALUE;
		//int index=0;
		//double aux;
		//for (int i=0; i < externalList.size(); i++){
		//	aux = distance.getDistance(dominatedSolution, externalList.get(i));
		//	if (aux < bestDistance){
		//        bestDistance = aux;
		//        index = i;
		//	}
		//}
		double x0,x1,x2, lambda;
		for (int i=0; i < oldSolution.getNumberOfVariables(); i++){
			x0 = (double) oldSolution.getVariableValue(i);
			x1 = (double) newSolution.getVariableValue(i);
			lambda = rnd.nextDouble();
			x2 = x1 + lambda * (x1-x0);
			x2 = Double.max(x2, ((DoubleSolution)oldSolution).getLowerBound(i));
			x2 = Double.min(x2, ((DoubleSolution)oldSolution).getUpperBound(i));
			advancedSolution.setVariableValue(i, x2);
		}
		
		return (S) advancedSolution;
	}
	
	protected List<S> calculateAdvancedSolutions() {
		//new SolutionListOutput(dominatedPopulation)
	    //.setSeparator("\t")
	    //.setVarFileOutputContext(new DefaultFileOutputContext("VAR_DOM"+Integer.toString(iterations)+".tsv"))
	    //.setFunFileOutputContext(new DefaultFileOutputContext("FUN_DOM"+Integer.toString(iterations)+".tsv"))
	    //.print();
//		boolean isnotdominated;
		List<S> advancedPopulation = new ArrayList<>();
		for (int i=0; i < this.getPopulationSize(); i++){
			advancedPopulation.add(calculateAdvancedSolutions(oldPopulation.get(i),this.getAllPopulation().get(i)));
		}
		advancedPopulation = evaluatePopulation(advancedPopulation);
		this.evaluations += advancedPopulation.size();
//		for (int i=0; i < improvedPopulation.size(); i++){
//			isnotdominated = externalList.add((S)improvedPopulation.get(i));
//			if (isnotdominated)
//				numberOfNonComparableSolutions ++;
//			else{
//				dominatedPopulation.add(improvedPopulation.get(i));
//			}
//		}
		System.out.println("Number of advanced solutions: "+ advancedPopulation.size());

		new SolutionListOutput(advancedPopulation)
	    .setSeparator("\t")
	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_AD"+Integer.toString(iterations)+".tsv"))
	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_AD"+Integer.toString(iterations)+".tsv"))
	    .print();
		return advancedPopulation;
	}
	
	protected void fitExternalList(){
		if (localSearch.getExternalList().size() > this.getPopulationSize()){
			Iterator<S> iterator = localSearch.getExternalList().getSolutionList().iterator();
			while (iterator.hasNext()) {
		        S listIndividual = iterator.next();
		        if (dist2refPoint.getAttribute(listIndividual) > maxDist2refPoint){
		        	iterator.remove();
		        	System.out.println(dist2refPoint.getAttribute(listIndividual));
		        	System.out.println(maxDist2refPoint);
		        }
			}
		}
	}
	
	protected S improveWithSASS() {
		// applied to allPopulation P
		// First of all: select the best individual in P, i. e., the one with the minimum distance
		// to the reference point
//		S bestIndividual ;
//		int indexBest = 0;
//		double smallestDistance = Double.MAX_VALUE;
//		for (int i = 0; i < this.getPopulationSize(); i++) {
//			S individual = this.getAllPopulation().get(i);
//			double aux = (double) dist2refPoint.getAttribute(individual); //CORREGIDO
//			if (aux < smallestDistance){
//				smallestDistance = aux;
//				indexBest = i;
//			}
//		}
		S bestIndividual = this.getAllPopulation().get(bestIndex2refPoint);
		
		S improvedIndividual ;
		double rad;
		double highestDistance = 0.0;
		int indexRad = 0;
		// Calculate radius
		for (int j = 0; j< this.getMaxPopulationSize();j++){
			double aux = distance.getDistance(bestIndividual, this.getPopulation().get(j));
			if (aux > highestDistance){
				highestDistance = aux;
			    indexRad = j;
			}
		}
		rad=distDomain.getNormalizedDistance((DoubleSolution)bestIndividual, (DoubleSolution) this.getPopulation().get(indexRad));
		localSearch.setRadius(rad);
		System.out.println("RADIUS: " + rad);
		improvedIndividual = localSearch.execute(bestIndividual);
		this.evaluations += localSearch.getEvaluations();
		//this.setPopulation(SASSPopulation);
		
		//List<S> externalListSASS = localSearch.getExternalList();
		//System.out.println("externalListPopulationSize= "+ externalListPopulation.size());
		//new SolutionListOutput(externalListPopulation)
	    //.setSeparator("\t")
	    //.setVarFileOutputContext(new DefaultFileOutputContext("VAR"+Integer.toString(iterations)+".tsv"))
	    //.setFunFileOutputContext(new DefaultFileOutputContext("FUN"+Integer.toString(iterations)+".tsv"))
	    //.print();
		////System.out.println("Number of improved solutions: "+ improvedPopulation.size());

		////new SolutionListOutput(improvedPopulation)
	    ////.setSeparator("\t")
	    ////.setVarFileOutputContext(new DefaultFileOutputContext("VAR_IMPR"+Integer.toString(iterations)+".tsv"))
	    ////.setFunFileOutputContext(new DefaultFileOutputContext("FUN_IMPR"+Integer.toString(iterations)+".tsv"))
	    ////.print();
		
		return improvedIndividual;
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
		return "NDWASFGA" ;
	}

	@Override public String getDescription() {
		return "Non Dominated Weighting Achievement Scalarizing Function Genetic Algorithm" ;
	}
}
