//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package org.uma.jmetal.algorithm.multiobjective.smsemoa;

import org.uma.jmetal.algorithm.impl.AbstractGeneticAlgorithm;
import org.uma.jmetal.algorithm.multiobjective.mombi.util.ASFWASFGA;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.qualityindicator.impl.Hypervolume;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.DominanceRanking;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
@SuppressWarnings("serial")
public class PAR_SMSEMOA<S extends Solution<?>> extends AbstractGeneticAlgorithm<S, List<S>> {
  protected final int maxEvaluations;
  protected final double offset ;

  protected int evaluations;

  private Hypervolume<S> hypervolume;
  
  List<Double> referencePoint = null;
  double[][] weights = {{1,1,1}};
  final ASFWASFGA<S> achievementScalarizingFunction;
  
  private final List<Double> utopiaPoint;
  private final List<Double> nadirPoint;

  /**
   * Constructor
   */
  public PAR_SMSEMOA(Problem<S> problem, int maxEvaluations, int populationSize, double offset,
      CrossoverOperator<S> crossoverOperator, MutationOperator<S> mutationOperator,
      SelectionOperator<List<S>, S> selectionOperator, Hypervolume<S> hypervolumeImplementation,List<Double> referencePoint) {
    super(problem) ;
    this.maxEvaluations = maxEvaluations;
    setMaxPopulationSize(populationSize);

    this.offset = offset ;

    this.crossoverOperator = crossoverOperator;
    this.mutationOperator = mutationOperator;
    this.selectionOperator = selectionOperator;

    this.hypervolume = hypervolumeImplementation ;
    this.referencePoint = referencePoint;
    
    
    this.achievementScalarizingFunction  = new ASFWASFGA<>(weights,referencePoint);
    
    this.nadirPoint     = new ArrayList<Double>(this.getProblem().getNumberOfObjectives());
	this.initializeNadirPoint(this.getProblem().getNumberOfObjectives());
	this.utopiaPoint = new ArrayList<Double>(this.getProblem().getNumberOfObjectives());
	this.initializeUtopiaPoint(this.getProblem().getNumberOfObjectives());
  }

  @Override protected void initProgress() {
    evaluations = getMaxPopulationSize() ;
  }

  @Override protected void updateProgress() {
    evaluations++ ;
  }

  @Override protected boolean isStoppingConditionReached() {
    return evaluations >= maxEvaluations ;
  }

  @Override protected List<S> evaluatePopulation(List<S> population) {
    for (S solution : population) {
      getProblem().evaluate(solution);
    }
    return population ;
  }

  @Override protected List<S> selection(List<S> population) {
    List<S> matingPopulation = new ArrayList<>(2);
    for (int i = 0; i < 2; i++) {
      S solution = selectionOperator.execute(population);
      matingPopulation.add(solution);
    }

    return matingPopulation;
  }

  @Override protected List<S> reproduction(List<S> population) {
    List<S> offspringPopulation = new ArrayList<>(1);

    List<S> parents = new ArrayList<>(2);
    parents.add(population.get(0));
    parents.add(population.get(1));

    List<S> offspring = crossoverOperator.execute(parents);

    mutationOperator.execute(offspring.get(0));

    offspringPopulation.add(offspring.get(0));
    return offspringPopulation;
  }

  @Override protected List<S> replacement(List<S> population, List<S> offspringPopulation) {

	
    List<S> jointPopulation = new ArrayList<>();
    jointPopulation.addAll(population);
    jointPopulation.addAll(offspringPopulation);
    List<S> resultPopulation = new ArrayList<>() ;
    
    int m = jointPopulation.get(0).getNumberOfObjectives();
//    double[] min = new double[m];
    double[] max = new double[m];
	initializeBounds(this.getProblem().getNumberOfObjectives());
	updateNadirPoint(jointPopulation);
	updateUtopiaPoint(jointPopulation);
    if (jointPopulation.size() > this.getMaxPopulationSize()){
    this.achievementScalarizingFunction.setNadir(getNadirPoint(jointPopulation));
	this.achievementScalarizingFunction.setUtopia(getUtopiaPoint(jointPopulation));
    double ASFmin = achievementScalarizingFunction.evaluate(jointPopulation.get(0),0);
    //int indexMin = 0;
    for (int i =1; i < jointPopulation.size(); i++){
    	double ASFvalue = achievementScalarizingFunction.evaluate(jointPopulation.get(i),0);
    	if(ASFvalue < ASFmin){
    		ASFmin = ASFvalue;
    	//	indexMin = i;
    	}
    }
    System.out.println("Smin "+ ASFmin);
    for (int nobj=0; nobj < m; nobj++){
    	List<Double> refPointAux = new ArrayList<>();
    	for (int ei=0; ei < m; ei++){
    		refPointAux.add(referencePoint.get(ei));
    	}
    	refPointAux.set(nobj, referencePoint.get(nobj)+(nadirPoint.get(nobj)-utopiaPoint.get(nobj))*ASFmin);
    	ASFWASFGA<S> auxASFunction = new ASFWASFGA<>(weights,refPointAux);
    	auxASFunction.setNadir(getNadirPoint(jointPopulation));
    	auxASFunction.setUtopia(getUtopiaPoint(jointPopulation));
    	double auxASFmin = auxASFunction.evaluate(jointPopulation.get(0),0);
        int indexMin = 0;
        for (int i =1; i < jointPopulation.size(); i++){
        	double auxASFvalue = auxASFunction.evaluate(jointPopulation.get(i),0);
        	if(auxASFvalue < auxASFmin){
        		auxASFmin = auxASFvalue;
        		indexMin = i;
        	}
        }
        max[nobj] = jointPopulation.get(indexMin).getObjective(nobj);
        System.out.println("ZauxREF "+ refPointAux.get(0)+" "+ refPointAux.get(1)+" "+ refPointAux.get(2));
        System.out.println(max[nobj]);
    }
    
    List<S> ROIPopulation = new ArrayList<>();
    List<S> noROIPopulation = new ArrayList<>();
    double aux;
//    double[] min = {0.0, 0.0};//, 10.0
//    double[] max = {0.3, 0.04};//, 55.0
    for (int i =0; i < jointPopulation.size(); i++){
    	boolean isIntoROI = true;
    	for (int nobj=0; nobj < m; nobj++){
    		aux = jointPopulation.get(i).getObjective(nobj);
    		if (aux > max[nobj]){
    			isIntoROI = false;
    		}
    	}
    	if (isIntoROI){
    		ROIPopulation.add(jointPopulation.get(i));
//    		jointPopulation.remove(jointPopulation.get(i));
    	}else
    		noROIPopulation.add(jointPopulation.get(i));
    }
    System.out.println("Número de individuos en ROI");
    System.out.println(ROIPopulation.size());

    

//    while ((ROIPopulation.size() < this.getMaxPopulationSize())&&(ROIPopulation.size() < population.size()+offspringPopulation.size())){
//    	ASFmin = achievementScalarizingFunction.evaluate(noROIPopulation.get(0),0);
//    	int indexMin = 0;
//    	for (int i =1; i < noROIPopulation.size(); i++){
//        	double ASFvalue = achievementScalarizingFunction.evaluate(noROIPopulation.get(i),0);
//        	if(ASFvalue < ASFmin){
//        		ASFmin = ASFvalue;
//        		indexMin = i;
//        	}
//        }
//    	ROIPopulation.add(noROIPopulation.get(indexMin));
//    	noROIPopulation.remove(noROIPopulation.get(indexMin));
//    }
    if (ROIPopulation.size() > this.getMaxPopulationSize()){
    	Ranking<S> ranking = computeRanking(ROIPopulation);
    	List<S> lastSubfront = ranking.getSubfront(ranking.getNumberOfSubfronts()-1) ;
    	
    	//System.out.println("Dentro del IF: Va a eliminar un individuo de la población");
    	//System.out.println(evaluations);

    	lastSubfront = hypervolume.computeHypervolumeContribution(lastSubfront, ROIPopulation) ;
    
    	for (int i = 0; i < ranking.getNumberOfSubfronts()-1; i++) {
    		for (S solution : ranking.getSubfront(i)) {
    			resultPopulation.add(solution);
    		}
    	}

    	for (int i = 0; i < lastSubfront.size()-1; i++) {
    		resultPopulation.add(lastSubfront.get(i)) ;
    	}
    }
    else
    	resultPopulation = ROIPopulation;
    }
    else{
    	resultPopulation = jointPopulation;
    }

    //if (evaluations%100==0){
		new SolutionListOutput(resultPopulation)
		.setSeparator("\t")
		.setVarFileOutputContext(new DefaultFileOutputContext("VAR_POP"+Integer.toString(evaluations)+".tsv"))
		.setFunFileOutputContext(new DefaultFileOutputContext("FUN_POP"+Integer.toString(evaluations)+".tsv"))
		.print();
	//}
    	
    return resultPopulation ;
  }

  @Override public List<S> getResult() {
    return getPopulation();
  }

  protected Ranking<S> computeRanking(List<S> solutionList) {
    Ranking<S> ranking = new DominanceRanking<S>();
    ranking.computeRanking(solutionList);

    return ranking;
  }
  
	public List<Double> getUtopiaPoint() {
		return this.utopiaPoint;
	}

	public List<Double> getNadirPoint() {
		return this.nadirPoint;
	}
	
	public List<Double> getUtopiaPoint(List<S> solutionList) {
		this.updateUtopiaPoint(solutionList);
		return this.utopiaPoint;
	}

	public List<Double> getNadirPoint(List<S> solutionList) {
		this.updateNadirPoint(solutionList);
		return this.nadirPoint;
	}


	private void initializeUtopiaPoint(int size) {
		for (int i = 0; i < size; i++)
			this.getUtopiaPoint().add(Double.POSITIVE_INFINITY);
	}

	private void initializeNadirPoint(int size) {
		for (int i = 0; i < size; i++)
			this.getNadirPoint().add(Double.NEGATIVE_INFINITY);
	}
	
	public void initializeBounds(int size) {
		this.initializeNadirPoint(size);
		this.initializeUtopiaPoint(size);
	}

	protected void updateUtopiaPoint(S s) {
		for (int i = 0; i < s.getNumberOfObjectives(); i++)
			this.getUtopiaPoint().set(i, Math.min(this.getUtopiaPoint().get(i),s.getObjective(i)));
	}

	private void updateNadirPoint(S s) {
		for (int i = 0; i < s.getNumberOfObjectives(); i++)
			this.getNadirPoint().set(i, Math.max(this.getNadirPoint().get(i),s.getObjective(i)));
	}

	public void updateUtopiaPoint(List<S> population) {
		for (S solution : population)
			this.updateUtopiaPoint(solution);
	}

	public void updateNadirPoint(List<S> population) {
		for (S solution : population)
			this.updateNadirPoint(solution);
	}

  @Override public String getName() {
    return "SMSEMOA" ;
  }

  @Override public String getDescription() {
    return "S metric selection EMOA" ;
  }
}
