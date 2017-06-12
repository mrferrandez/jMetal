package org.uma.jmetal.algorithm.multiobjective.wasfga;

import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.localsearch.SASSLocalSearch;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.AlgorithmBuilder;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.util.List;

/**
 * @author Miriam R. Ferrández <mrferrandez@ual.es>
 */
public class WASFGABuilder<S extends Solution<?>> implements AlgorithmBuilder<WASFGA<S>> {
  public enum WASFGAVariant {WASFGA, WASFGA_ad, NDSASSWASFGA, PAR_SASS, PARWASFGA, NDrankWASFGA}

  /**
   * WASFGABuilder class
   */
  private final Problem<S> problem;
  private int maxEvaluations;
  private int populationSize;
  private CrossoverOperator<S>  crossoverOperator;
  private MutationOperator<S> mutationOperator;
  private SelectionOperator<List<S>, S> selectionOperator;
  private SolutionListEvaluator<S> evaluator;

  private WASFGAVariant variant;
  
  private List<Double> referencePoint = null;
  
  private int improvementRounds = 7;
  private SASSLocalSearch<S> localSearch;

  /**
   * WASFGABuilder constructor
   */

public WASFGABuilder(Problem<S> problem, CrossoverOperator<S> crossoverOperator,
      MutationOperator<S> mutationOperator, List<Double> referencePoint) {
    this.problem = problem;
    maxEvaluations = 25000;
    populationSize = 100;
    this.crossoverOperator = crossoverOperator ;
    this.mutationOperator = mutationOperator ;
    this.referencePoint = referencePoint;
    selectionOperator = new BinaryTournamentSelection<S>(new RankingAndCrowdingDistanceComparator<S>()) ;
    evaluator = new SequentialSolutionListEvaluator<S>();

    this.variant = WASFGAVariant.WASFGA ;
    
    this.localSearch = new SASSLocalSearch<S>(improvementRounds, problem);
  }

  public WASFGABuilder<S> setMaxEvaluations(int maxEvaluations) {
    if (maxEvaluations < 0) {
      throw new JMetalException("maxEvaluations is negative: " + maxEvaluations);
    }
    this.maxEvaluations = maxEvaluations;

    return this;
  }

  public WASFGABuilder<S> setPopulationSize(int populationSize) {
    if (populationSize < 0) {
      throw new JMetalException("Population size is negative: " + populationSize);
    }

    this.populationSize = populationSize;

    return this;
  }

  public WASFGABuilder<S> setSelectionOperator(SelectionOperator<List<S>, S> selectionOperator) {
    if (selectionOperator == null) {
      throw new JMetalException("selectionOperator is null");
    }
    this.selectionOperator = selectionOperator;

    return this;
  }

  public WASFGABuilder<S> setSolutionListEvaluator(SolutionListEvaluator<S> evaluator) {
    if (evaluator == null) {
      throw new JMetalException("evaluator is null");
    }
    this.evaluator = evaluator;

    return this;
  }


  public WASFGABuilder<S> setVariant(WASFGAVariant variant) {
    this.variant = variant;

    return this;
  }

  public WASFGA<S> build() {
    WASFGA<S> algorithm = null ;
    if (variant.equals(WASFGAVariant.WASFGA)) {
      algorithm = new WASFGA<S>(problem, populationSize, maxEvaluations,crossoverOperator,
          mutationOperator, selectionOperator, evaluator, referencePoint);
    } else if (variant.equals(WASFGAVariant.WASFGA_ad)) {
      algorithm = new WASFGA_ad<S>(problem,  populationSize, maxEvaluations,crossoverOperator,
          mutationOperator, selectionOperator, evaluator, referencePoint);
    } else if (variant.equals(WASFGAVariant.NDSASSWASFGA)) {
      algorithm = new NDSASSWASFGA<S>(problem, populationSize, maxEvaluations, localSearch, crossoverOperator,
          mutationOperator, selectionOperator, evaluator, referencePoint);
    } else if (variant.equals(WASFGAVariant.PAR_SASS)) {
        algorithm = new PAR_SASS<S>(problem, populationSize, maxEvaluations, localSearch, crossoverOperator,
                mutationOperator, selectionOperator, evaluator, referencePoint);
    } else if (variant.equals(WASFGAVariant.PARWASFGA)) {
        algorithm = new PAR_WASFGA<S>(problem, populationSize, maxEvaluations, crossoverOperator,
                mutationOperator, selectionOperator, evaluator, referencePoint);
    } else if (variant.equals(WASFGAVariant.NDrankWASFGA)) {
        algorithm = new NDrankWASFGA<S>(problem, populationSize, maxEvaluations, crossoverOperator,
                mutationOperator, selectionOperator, evaluator, referencePoint);
    }

    return algorithm ;
  }

  /* Getters */
  public Problem<S> getProblem() {
    return problem;
  }

  public int getMaxIterations() {
    return maxEvaluations;
  }

  public int getPopulationSize() {
    return populationSize;
  }

  public CrossoverOperator<S> getCrossoverOperator() {
    return crossoverOperator;
  }

  public MutationOperator<S> getMutationOperator() {
    return mutationOperator;
  }

  public SelectionOperator<List<S>, S> getSelectionOperator() {
    return selectionOperator;
  }

  public SolutionListEvaluator<S> getSolutionListEvaluator() {
    return evaluator;
  }
}
