package org.uma.jmetal.algorithm.multiobjective.wasfga;

import org.uma.jmetal.algorithm.multiobjective.mombi.AbstractMOMBI;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.ASFWASFGA;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.AbstractUtilityFunctionsSet;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.Normalizer;
import org.uma.jmetal.algorithm.multiobjective.wasfga.util.PrincipalComponentAnalysis;
import org.uma.jmetal.algorithm.multiobjective.wasfga.util.WASFGARanking;
import org.uma.jmetal.algorithm.multiobjective.wasfga.util.WeightVector;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.distance.impl.EuclideanDistanceBetweenSolutionsInObjectiveSpace;
import org.uma.jmetal.util.distance.impl.HausdorffDistanceBetweenSetsInObjectiveSpace;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
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
public class WASFGA_V<S extends Solution<?>> extends AbstractMOMBI<S> {
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	protected int maxEvaluations;
	protected int evaluations;
	protected Normalizer normalizer;
	
	final ASFWASFGA<S> achievementScalarizingFunction;
	List<Double> referencePoint = null;
	
	protected EuclideanDistanceBetweenSolutionsInObjectiveSpace<S> distance ;
	protected HausdorffDistanceBetweenSetsInObjectiveSpace<S,List<S>> hausdorff;
	protected double[] consecutiveHausdorff = {Double.MAX_VALUE, Double.MAX_VALUE};
	protected double tol = 1.0E-5;
	protected PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
	protected int nreducedvar = 2;
	protected double epsMesh = 1.0;
	private JMetalRandom randomGenerator ;

	/**
	 * Constructor
	 *
	 * @param problem
	 *            Problem to solve
	 */
	public WASFGA_V(Problem<S> problem,
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
		this.maxEvaluations = maxIterations;
		
		distance = new EuclideanDistanceBetweenSolutionsInObjectiveSpace<S>() ;
		hausdorff = new HausdorffDistanceBetweenSetsInObjectiveSpace<S, List<S>>(distance);
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
	protected boolean isStoppingConditionReached() {
		boolean stop = false;
		if (this.iterations >= this.maxEvaluations){
			stop = true;
		} else{
			if ((consecutiveHausdorff[0]<tol)&&(consecutiveHausdorff[1]<tol)){
				stop = true;
			}
		}
		
		return stop;
	}

	@Override
	public void specificMOEAComputations() {
		initializeBounds(this.getProblem().getNumberOfObjectives());
		updateNadirPoint(this.getPopulation());
		updateReferencePoint(this.getPopulation());
//		new SolutionListOutput(this.getPopulation())
//	    .setSeparator("\t")
//	    .setVarFileOutputContext(new DefaultFileOutputContext("VAR_POP"+Integer.toString(iterations)+".tsv"))
//	    .setFunFileOutputContext(new DefaultFileOutputContext("FUN_POP"+Integer.toString(iterations)+".tsv"))
//	    .print();
	}
	
	protected double[][] proyectPCA(List<S> solutionList) {
		int nvar=solutionList.get(0).getNumberOfVariables();
		pca.setup(solutionList.size(), nvar);
		double[][] varsample= new double[solutionList.size()][nvar];
		for (int sample=0; sample<solutionList.size(); sample++){
			for (int index=0; index<nvar; index++){
				varsample[sample][index]=(double) solutionList.get(sample).getVariableValue(index);
			}
			pca.addSample(varsample[sample]);
		}
		//int nreducedvar=2;
		pca.computeBasis(nreducedvar);
		double[][] reducedSpaceList = new double [solutionList.size()][nreducedvar];
		for (int sample=0; sample<solutionList.size(); sample++){
			reducedSpaceList[sample]=pca.sampleToEigenSpace(varsample[sample]);
		}
		return reducedSpaceList;
	}
	
//	protected boolean belongToCell(double[] minCell, double[] maxCell, double[] individual) {
//		boolean belong = true;
//		for (int naxis=0; naxis<individual.length; naxis++){
//			if ((individual[naxis]>maxCell[naxis])&&(individual[naxis]<minCell[naxis])){
//				belong =false;
//				break;
//			}
//		}
//		return belong;
//	}
	
	@Override
	protected List<S> reproduction(List<S> population) {
		List<S> offspringPopulation = new ArrayList<>(this.getMaxPopulationSize());
		double[][] reducedVar = proyectPCA(population);
		double[] max = reducedVar[0];
		double[] min = reducedVar[0];
		for (int index=1; index<population.size(); index++){
			for (int naxis=0; naxis<nreducedvar; naxis++){
				if(reducedVar[index][naxis]>max[naxis]){
					max[0]=reducedVar[index][naxis];
				}
				else if (reducedVar[index][naxis]<min[naxis]){
					min[0]=reducedVar[index][naxis];
				}
			}
		}
		double[] width = new double[nreducedvar];
		for (int naxis=0; naxis<nreducedvar; naxis++){
			width[naxis] = max[naxis]-min[naxis];
		}
		// Desplazamiento en la dirección ortogonal. Puede ser "-" al azar
		double move = width[1];
		if (randomGenerator.nextDouble()<=0.5){
			move = 0.25*width[1];
		}
		else{
			move = -0.25*width[1]; 
		} 
		// Anchura agrandada
		width[0] = 1.25 * width[0];
		// Desplazamiento en la dirección ortogonal aplicado
		width[1] = width[1] + move;
		
		// VORONOI
		int nMesh = (int) epsMesh*population.size();
		double[] dMesh = new double[nreducedvar];
		for (int naxis=0; naxis<nreducedvar; naxis++){ //división de cada eje en intervalos
			dMesh[naxis] = width[naxis]/nMesh;
		}
		double[][] mesh = new double[nreducedvar][nMesh+1];
		//for (int naxis=0; naxis<nreducedvar; naxis++){
		mesh[0][0] = 1.25*min[0];
		mesh[1][0] = min[0] + move;
		//}
		for (int naxis=0; naxis<nreducedvar; naxis++){
			for (int ncell=1; ncell<nMesh+1; ncell++){
				mesh[naxis][ncell] = mesh[naxis][ncell-1]+dMesh[naxis];
			}
		}
		double[][] allCells = new double[nMesh+1][nMesh+1];
		//int[] pos = new int[nreducedvar];
		for (int i=0; i<population.size(); i++){
			for (int naxis=0; naxis<nreducedvar; naxis++){
				for (int pos=(int)(nMesh+1); pos>1; pos=pos/2){
					if(reducedVar[i][naxis]<mesh[pos][naxis]){
						
					}
				}
			}
		}
		
//		for (int i = 0; i < this.getMaxPopulationSize(); i += 2) {
//			List<S> parents = new ArrayList<>(2);
//			int parent1Index = JMetalRandom.getInstance().nextInt(0, this.getMaxPopulationSize()-1);
//			int parent2Index = JMetalRandom.getInstance().nextInt(0, this.getMaxPopulationSize()-1);
//			while (parent1Index==parent2Index)
//				parent2Index = JMetalRandom.getInstance().nextInt(0, this.getMaxPopulationSize()-1);
//			parents.add(population.get(parent1Index));
//			parents.add(population.get(parent2Index));
//
//			//List<S> offspring = crossoverOperator.execute(parents);
//
//			//mutationOperator.execute(offspring.get(0));
//			//mutationOperator.execute(offspring.get(1));
//
//			//offspringPopulation.add(offspring.get(0));
//			//offspringPopulation.add(offspring.get(1));
//		}
		return offspringPopulation;
	}

	@Override
	protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
		List<S> jointPopulation = new ArrayList<>();
		jointPopulation.addAll(population);
		jointPopulation.addAll(offspringPopulation);
		Ranking<S> ranking = computeRanking(jointPopulation);
		List<S> newPopulation = selectBest(ranking);
		consecutiveHausdorff[0] = consecutiveHausdorff[1];
		consecutiveHausdorff[1] = hausdorff.getDistance(population, newPopulation);
		return newPopulation;
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
