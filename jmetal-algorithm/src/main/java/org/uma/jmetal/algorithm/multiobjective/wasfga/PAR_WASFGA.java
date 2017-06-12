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
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
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
 */
public class PAR_WASFGA<S extends Solution<?>> extends AbstractMOMBI<S> {
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	protected int maxEvaluations;
	protected int evaluations;
	protected Normalizer normalizer;
	
	final ASFWASFGA<S> achievementScalarizingFunction;
	List<Double> referencePoint = null;
	
	protected double [][] weights;

	/**
	 * Constructor
	 *
	 * @param problem
	 *            Problem to solve
	 */
	public PAR_WASFGA(Problem<S> problem,
								int populationSize,
								int maxIterations,
								CrossoverOperator<S> crossoverOperator,
								MutationOperator<S> mutationOperator,
								SelectionOperator<List<S>, S> selectionOperator,
								SolutionListEvaluator<S> evaluator,
								List<Double> referencePoint) {

		super(problem,maxIterations,crossoverOperator,mutationOperator,selectionOperator,evaluator);
		setMaxPopulationSize(populationSize);
		this.referencePoint 				= referencePoint;
		this.achievementScalarizingFunction =  createUtilityFunction();
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
		this.weights = weights;
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
	}

	@Override
	protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
		List<S> jointPopulation = new ArrayList<>();
		jointPopulation.addAll(population);
		jointPopulation.addAll(offspringPopulation);
		List<S> ROIPopulation = new ArrayList<>();
		ROIPopulation = computeROI(jointPopulation);
		List<S> resultPopulation = new ArrayList<>();
		if (ROIPopulation.size() > this.getMaxPopulationSize()){
			Ranking<S> ranking = computeRanking(ROIPopulation);
			resultPopulation = selectBest(ranking);
		}
		else
			resultPopulation = ROIPopulation;
		
		return resultPopulation;
	}
	
	protected List<S> computeROI(List<S> population) {
		int m = population.get(0).getNumberOfObjectives();
		double[] max = new double[m];
//		initializeBounds(this.getProblem().getNumberOfObjectives());
//		updateNadirPoint(jointPopulation);
//		updateReferencePoint(jointPopulation);
		this.achievementScalarizingFunction.setNadir(getNadirPoint(population));
		this.achievementScalarizingFunction.setUtopia(getReferencePoint(population));
		double ASFmin = Double.MAX_VALUE;
	    for (int i =0; i < population.size(); i++){
	    	for (int weight=0; weight < achievementScalarizingFunction.getVectorSize(); weight++){
	    		double ASFvalue = achievementScalarizingFunction.evaluate(population.get(i),weight);
	    		if(ASFvalue < ASFmin){
	    			ASFmin = ASFvalue;
	    		}
	    	}
	    }
	    System.out.println("Smin "+ ASFmin);
	     
	    for (int nobj=0; nobj < m; nobj++){
	    	List<Double> refPointAux = new ArrayList<>();
	    	for (int ei=0; ei < m; ei++){
	    		refPointAux.add(referencePoint.get(ei));
	    	}
	    	refPointAux.set(nobj, referencePoint.get(nobj)+(this.getNadirPoint().get(nobj)-this.getReferencePoint().get(nobj))*ASFmin);
	    	ASFWASFGA<S> auxASFunction = new ASFWASFGA<>(this.weights,refPointAux);
	    	auxASFunction.setNadir(getNadirPoint(population));
	    	auxASFunction.setUtopia(getReferencePoint(population));
	    	double auxASFmin = Double.MAX_VALUE;//auxASFunction.evaluate(population.get(0),0);
	        int indexMin = 0;
	        for (int i =0; i < population.size(); i++){
	        	for (int weight=0; weight < achievementScalarizingFunction.getVectorSize(); weight++){
	        		double auxASFvalue = auxASFunction.evaluate(population.get(i),weight);
	        		if(auxASFvalue < auxASFmin){
	        			auxASFmin = auxASFvalue;
	        			indexMin = i;
	        		}
	        	}
	        }
	        max[nobj] = population.get(indexMin).getObjective(nobj);
	        System.out.println("ZauxREF "+ refPointAux.get(0)+" "+ refPointAux.get(1)+" "+ refPointAux.get(2));
	        System.out.println(max[nobj]);
	    }
	    
	    List<S> ROIPopulation = new ArrayList<>();
	    List<S> noROIPopulation = new ArrayList<>();
	    double aux;
	    for (int i =0; i < population.size(); i++){
	    	boolean isIntoROI = true;
	    	for (int nobj=0; nobj < m; nobj++){
	    		aux = population.get(i).getObjective(nobj);
	    		if (aux > max[nobj]){
	    			isIntoROI = false;
	    		}
	    	}
	    	if (isIntoROI){
	    		ROIPopulation.add(population.get(i));
//	    		jointPopulation.remove(jointPopulation.get(i));
	    	}else
	    		noROIPopulation.add(population.get(i));
	    }
	    System.out.println("Número de individuos en ROI");
	    System.out.println(ROIPopulation.size());
	    
	    while (ROIPopulation.size() < this.getMaxPopulationSize()){
	    	ASFmin = Double.MAX_VALUE;
	    	int indexMin = 0;
	    	for (int i =1; i < noROIPopulation.size(); i++){
	    		for (int weight=0; weight < achievementScalarizingFunction.getVectorSize(); weight++){
	    			double ASFvalue = achievementScalarizingFunction.evaluate(noROIPopulation.get(i),weight);
	    			if(ASFvalue < ASFmin){
	    				ASFmin = ASFvalue;
	    				indexMin = i;
	    			}
	    		}
	    	}
	    	ROIPopulation.add(noROIPopulation.get(indexMin));
	    	noROIPopulation.remove(noROIPopulation.get(indexMin));
	    }
	    
		return ROIPopulation;
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
		return "WASFGA" ;
	}

	@Override public String getDescription() {
		return "Weighting Achievement Scalarizing Function Genetic Algorithm" ;
	}
}
