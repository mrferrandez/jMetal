package org.uma.jmetal.algorithm.multiobjective.wasfga;

import org.uma.jmetal.algorithm.impl.AbstractEvolutionaryAlgorithm;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.ASFWASFGA;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.AbstractUtilityFunctionsSet;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.Normalizer;
import org.uma.jmetal.algorithm.multiobjective.wasfga.util.WASFGARanking;
import org.uma.jmetal.algorithm.multiobjective.wasfga.util.WeightVector;
import org.uma.jmetal.operator.impl.localsearch.SASSLocalSearch;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.solutionattribute.Ranking;

import java.util.ArrayList;
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
 * Modified version of WASFGA
 * Local Search Operator SASS included instead of Mutation and Crossover Operators
 * @author Miriam R. Ferrandez
 * 
 */
public class SASSWASFGA<S extends Solution<?>, Result> extends AbstractEvolutionaryAlgorithm<S,List<S>> {
	/**
	 *
	 */
	private final int maxIterations;
	private int maxPopulationSize ;

	private int iterations = 0;
	private final SolutionListEvaluator<S> evaluator;
	private final List<Double> idealPoint;
	private final List<Double> nadirPoint;
	
	private static final long serialVersionUID = 1L;
	protected int maxEvaluations;
	protected int evaluations;
	protected Normalizer normalizer;
	
	protected SASSLocalSearch<S> localSearch;
	
	final ASFWASFGA<S> achievementScalarizingFunction;
	List<Double> referencePoint = null;

	/**
	 * Constructor
	 *
	 * @param problem
	 *            Problem to solve
	 */
	public SASSWASFGA(Problem<S> problem,
								int populationSize,
								int maxIterations,
								SASSLocalSearch<S> localSearch,
								SolutionListEvaluator<S> evaluator,
								List<Double> referencePoint) {

	    setProblem(problem);
		this.maxIterations = maxIterations;
		this.maxPopulationSize = populationSize;
		this.referencePoint = referencePoint;
		this.localSearch = localSearch;
		this.achievementScalarizingFunction =  createUtilityFunction();
		
		this.evaluator = evaluator;

		this.nadirPoint     = new ArrayList<Double>(this.getProblem().getNumberOfObjectives());
		this.initializeNadirPoint(this.getProblem().getNumberOfObjectives());
		this.idealPoint = new ArrayList<Double>(this.getProblem().getNumberOfObjectives());
		this.initializeIdealPoint(this.getProblem().getNumberOfObjectives());
	}

	public ASFWASFGA<S> createUtilityFunction() {
		double [][] weights;
		if (referencePoint.size()==2)
			weights = WeightVector.initUniformWeights2D(0.005, this.maxPopulationSize);
		else{
			String filePath=System.getProperty("user.dir")+"/weightsWASFGA/W3D_"+Integer.toString(this.maxPopulationSize)+".dat";
			weights = WeightVector.getWeightsFromFile(filePath);
		}
		weights = WeightVector.invertWeights(weights,true);
		ASFWASFGA<S> aux = new ASFWASFGA<>(weights,referencePoint);

		return aux;
	}
	
	protected boolean populationIsNotFull(List<S> population) {
		return population.size() < this.maxPopulationSize;
	}
	
	@Override
	protected void initProgress() {
		this.iterations = 1;
	}

	@Override
	protected void updateProgress() {
		this.iterations+=1;
	}

	@Override
	protected boolean isStoppingConditionReached() {
		return this.iterations >= this.maxIterations;
	}

	@Override
	protected List<S> evaluatePopulation(List<S> population) {
		population = evaluator.evaluate(population, getProblem());

		return population;
	}

	/**
	 * Copied from AbstractGeneticAlgorithm
	 * This method implements a default scheme create the initial population of genetic algorithm
	 * @return
	 */
	protected List<S> createInitialPopulation() {
		List<S> population = new ArrayList<>(this.maxPopulationSize);
		for (int i = 0; i < this.maxPopulationSize; i++) {
			S newIndividual = getProblem().createSolution();
			population.add(newIndividual);
		}
		return population;
	}
	
	@Override
	protected List<S> selection(List<S> population) {
		//List<S> matingPopulation = new ArrayList<>(population.size());
		//for (int i = 0; i < this.maxPopulationSize; i++) {
		//	S solution = selectionOperator.execute(population);
		//	matingPopulation.add(solution);
		//}

		return population;//matingPopulation;
	}

	@SuppressWarnings("unchecked")
	@Override
	protected List<S> reproduction(List<S> population) {
		List<S> SASSPopulation = new ArrayList<>(this.maxPopulationSize);
		for (int i = 0; i < this.maxPopulationSize; i++) {
			SASSPopulation.add(localSearch.execute(population.get(i)));
			this.evaluations += localSearch.getEvaluations();
		}
		this.setPopulation(SASSPopulation);
		
		List<S> externalListPopulation = (List<S>) localSearch.getExternalList();
		System.out.println("externalListPopulationSize= "+ externalListPopulation.size());
		new SolutionListOutput(externalListPopulation)
	    .setSeparator("\t")
	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR"+Integer.toString(iterations)+".tsv"))
	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN"+Integer.toString(iterations)+".tsv"))
	    .print();
		
		return externalListPopulation;
	}
	
	@Override
	public void run() {
		List<S> offspringPopulation;
		//List<S> matingPopulation;

		this.setPopulation(createInitialPopulation());
		this.evaluatePopulation(this.getPopulation());
		initProgress();
		//specific GA needed computations
		this.specificMOEAComputations();
		while (!isStoppingConditionReached()) {
			System.out.println("NUEVA GENERACION "+ iterations);
			//matingPopulation = selection(this.getPopulation());
			offspringPopulation = reproduction(this.getPopulation());
			//offspringPopulation = evaluatePopulation(offspringPopulation);
			this.setPopulation(replacement(this.getPopulation(), offspringPopulation));
			updateProgress();
			// specific GA needed computations
			this.specificMOEAComputations();
		}
		
		System.out.println("Total number of evaluations: "+ this.evaluations);
	}

	public void specificMOEAComputations() {
		initializeBounds(this.getProblem().getNumberOfObjectives());
		updateNadirPoint(this.getPopulation());
		updateIdealPoint(this.getPopulation());
		//localSearch.clearExternalList(); NO BORRO LA LISTA 
	}

	@Override
	protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
		List<S> jointPopulation = new ArrayList<>();
		jointPopulation.addAll(population);
		jointPopulation.addAll(offspringPopulation);
		Ranking<S> ranking = computeRanking(jointPopulation);
		return selectBest(ranking);
	}
	
	protected Ranking<S> computeRanking(List<S> solutionList) {
		this.achievementScalarizingFunction.setNadir(getNadirPoint(solutionList));
		this.achievementScalarizingFunction.setUtopia(getIdealPoint(solutionList));
		Ranking<S> ranking = new WASFGARanking<>(this.achievementScalarizingFunction);
		ranking.computeRanking(solutionList);
		return ranking;
	}
	
	protected void addRankedSolutionsToPopulation(Ranking<S> ranking, int index, List<S> population) {
		population.addAll(ranking.getSubfront(index));
	}
	
	protected void addLastRankedSolutionsToPopulation(Ranking<S> ranking,int index, List<S>population) {
		List<S> front 	= ranking.getSubfront(index);
		int remain 		= this.maxPopulationSize - population.size();
		population.addAll(front.subList(0, remain));
	}
	
	protected List<S> selectBest(Ranking<S> ranking) {
		List<S> population = new ArrayList<>(this.maxPopulationSize);
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
		return (population.size()+ranking.getSubfront(index).size() < this.maxPopulationSize);
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
		return "WASFGA" ;
	}

	@Override public String getDescription() {
		return "Weighting Achievement Scalarizing Function Genetic Algorithm" ;
	}
	
	private void initializeIdealPoint(int size) {
		for (int i = 0; i < size; i++)
			this.idealPoint.add(Double.POSITIVE_INFINITY);
	}

	private void initializeNadirPoint(int size) {
		for (int i = 0; i < size; i++)
			this.nadirPoint.add(Double.NEGATIVE_INFINITY);
	}
	
	public void initializeBounds(int size) {
		this.initializeNadirPoint(size);
		this.initializeIdealPoint(size);
	}

	protected void updateIdealPoint(S s) {
		for (int i = 0; i < s.getNumberOfObjectives(); i++)
			this.idealPoint.set(i, Math.min(this.idealPoint.get(i),s.getObjective(i)));
	}

	private void updateNadirPoint(S s) {
		for (int i = 0; i < s.getNumberOfObjectives(); i++)
			this.nadirPoint.set(i, Math.max(this.nadirPoint.get(i),s.getObjective(i)));
	}

	public void updateIdealPoint(List<S> population) {
		for (S solution : population)
			this.updateIdealPoint(solution);
	}

	public void updateNadirPoint(List<S> population) {
		for (S solution : population)
			this.updateNadirPoint(solution);
	}
	
	public List<Double> getIdealPoint(List<S> solutionList) {
		this.updateIdealPoint(solutionList);
		return this.idealPoint;
	}

	public List<Double> getNadirPoint(List<S> solutionList) {
		this.updateNadirPoint(solutionList);
		return this.nadirPoint;
	}

}
