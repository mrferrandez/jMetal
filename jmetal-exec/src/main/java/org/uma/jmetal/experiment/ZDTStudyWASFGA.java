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

package org.uma.jmetal.experiment;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.algorithm.multiobjective.smpso.SMPSOBuilder;
import org.uma.jmetal.algorithm.multiobjective.spea2.SPEA2Builder;
import org.uma.jmetal.algorithm.multiobjective.wasfga.WASFGABuilder;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.problem.multiobjective.zdt.*;
import org.uma.jmetal.qualityindicator.impl.*;
import org.uma.jmetal.qualityindicator.impl.hypervolume.PISAHypervolume;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.archive.impl.CrowdingDistanceArchive;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;
import org.uma.jmetal.util.experiment.Experiment;
import org.uma.jmetal.util.experiment.ExperimentBuilder;
import org.uma.jmetal.util.experiment.component.*;
import org.uma.jmetal.util.experiment.util.TaggedAlgorithm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Example of experimental study based on solving the ZDT problems with four algorithms: 
 * WASFGA, PAR-WASFGA, ND-SASS-WASFGA, WASFGA-ad
 *
 * This experiment assumes that the reference Pareto front are known, so the names of files containing
 * them and the directory where they are located must be specified.
 *
 * Six quality indicators are used for performance assessment.
 *
 * The steps to carry out the experiment are:
 * 1. Configure the experiment
 * 2. Execute the algorithms
 * 3. Compute the quality indicators
 * 4. Generate Latex tables reporting means and medians
 * 5. Generate R scripts to produce latex tables with the result of applying the Wilcoxon Rank Sum Test
 * 6. Generate Latex tables with the ranking obtained by applying the Friedman test
 * 7. Generate R scripts to obtain boxplots
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 * modified by Miriam R. Ferrández <mrferrandez@ual.es>
 */

public class ZDTStudyWASFGA {
  private static final int INDEPENDENT_RUNS = 30 ;

  public static void main(String[] args) throws IOException {
    if (args.length != 1) {
      throw new JMetalException("Missing argument: experiment base directory") ;
    }
    String experimentBaseDirectory = args[0] ;

    List<Problem<DoubleSolution>> problemList = Arrays.<Problem<DoubleSolution>>asList(new ZDT1(), new ZDT2(), new ZDT3(), new ZDT4(), new ZDT6());
     //   ;

    List<String> referenceFrontFileNames = Arrays.asList("ZDT1.pf", "ZDT2.pf", "ZDT3.pf", "ZDT4.pf", "ZDT6.pf");//"ZDT2.pf",  "ZDT3.pf", "ZDT4.pf", "ZDT6.pf") ;

    List<TaggedAlgorithm<List<DoubleSolution>>> algorithmList = configureAlgorithmList(problemList, INDEPENDENT_RUNS) ;

    Experiment<DoubleSolution, List<DoubleSolution>> experiment =
        new ExperimentBuilder<DoubleSolution, List<DoubleSolution>>("ZDTStudyWASFGA")
            .setAlgorithmList(algorithmList)
            .setProblemList(problemList)
            .setReferenceFrontDirectory("/pareto_fronts")
            .setReferenceFrontFileNames(referenceFrontFileNames)
            .setExperimentBaseDirectory(experimentBaseDirectory)
            .setOutputParetoFrontFileName("FUN")
            .setOutputParetoSetFileName("VAR")
            .setIndicatorList(Arrays.asList(
                new Epsilon<DoubleSolution>(), new Spread<DoubleSolution>(), new GenerationalDistance<DoubleSolution>(),
                new PISAHypervolume<DoubleSolution>(),
                new InvertedGenerationalDistance<DoubleSolution>(),
                new InvertedGenerationalDistancePlus<DoubleSolution>()))
            .setIndependentRuns(INDEPENDENT_RUNS)
            .setNumberOfCores(8)
//            .setReferencePoint(configureReferencePoint())
            .build();

//    new ExecuteAlgorithms<>(experiment).run();
    new ComputeQualityIndicators<>(experiment).run() ;
    new GenerateLatexTablesWithStatistics(experiment).run() ;
    new GenerateWilcoxonTestTablesWithR<>(experiment).run() ;
    new GenerateFriedmanTestTables<>(experiment).run();
    new GenerateBoxplotsWithR<>(experiment).setRows(3).setColumns(3).setDisplayNotch().run() ;
  }

  /**
   * The algorithm list is composed of pairs {@link Algorithm} + {@link Problem} which form part of a
   * {@link TaggedAlgorithm}, which is a decorator for class {@link Algorithm}.
   *
   * @param problemList
   * @return
   */
  static List<TaggedAlgorithm<List<DoubleSolution>>> configureAlgorithmList(
      List<Problem<DoubleSolution>> problemList,
      int independentRuns) {
    List<TaggedAlgorithm<List<DoubleSolution>>> algorithms = new ArrayList<>() ;
    double[][] refPoint={{0.80, 0.60}, {0.80, 0.80}, {0.30, 0.80}, {0.99, 0.95}, {0.78, 0.61}};
    
//    List<Double> refPoint = configureReferencePoint();
//    referencePoint.add(0.0);
//    referencePoint.add(0.0);
    //referencePoint.add(30.0);
    
    //int[] popSizes={50, 100, 200};
    //int[] gen ={25, 50, 100};

    for (int run = 0; run < INDEPENDENT_RUNS; run++) {//independentRuns
// //     for (int i = 0; i < problemList.size(); i++) {
 //        double mutationProbability = 1.0 / problemList.get(i).getNumberOfVariables();
 //        double mutationDistributionIndex = 20.0;
 //   	 for (int npop=0; npop < popSizes.length; npop++){
 //   		 for (int g=0; g < gen.length; g++){
/*//    		Algorithm<List<DoubleSolution>> algorithm = new WASFGABuilder<>(problemList.get(i), new SBXCrossover(0.9, 20.0),
//                new PolynomialMutation(1.0 / problemList.get(i).getNumberOfVariables(), 20.0),refPoint)
//                .setMaxEvaluations(500)//gen[g]
//                .setPopulationSize(200)//popSizes[npop]
//                .setVariant(WASFGABuilder.WASFGAVariant.WASFGA)
//                .build();
//            algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, "WASFGA", problemList.get(i), run));
*/    		 //}
    	 //}
        
// //       }

      for (int i = 0; i < problemList.size(); i++) {
    	  Algorithm<List<DoubleSolution>> algorithm = new WASFGABuilder<>(problemList.get(i), new SBXCrossover(0.9, 20.0),
                  new PolynomialMutation(1.0 / problemList.get(i).getNumberOfVariables(), 20.0),configureReferencePoint(refPoint[i]))
                  .setMaxEvaluations(500)
                  .setPopulationSize(100)
                  .setVariant(WASFGABuilder.WASFGAVariant.WASFGA)
                  .build();
              algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, "WASFGA", problemList.get(i), run));
            }
      	

  for (int i = 0; i < problemList.size(); i++) {
    	  Algorithm<List<DoubleSolution>> algorithm = new WASFGABuilder<>(problemList.get(i), new SBXCrossover(0.9, 20.0),
                  new PolynomialMutation(1.0 / problemList.get(i).getNumberOfVariables(), 20.0),configureReferencePoint(refPoint[i]))
                  .setMaxEvaluations(500)
                  .setPopulationSize(100)
                  .setVariant(WASFGABuilder.WASFGAVariant.NDSASSWASFGA)
                  .build();
              algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, "NDSASSWASFGA", problemList.get(i), run));
      }
//      for (int i = 0; i < problemList.size(); i++) {
//    	  Algorithm<List<DoubleSolution>> algorithm = new WASFGABuilder<>(problemList.get(i), new SBXCrossover(0.9, 20.0),
//                  new PolynomialMutation(1.0 / problemList.get(i).getNumberOfVariables(), 20.0),configureReferencePoint(refPoint[i]))
//                  .setMaxEvaluations(500)
//                  .setPopulationSize(100)
//                  .setVariant(WASFGABuilder.WASFGAVariant.PAR_SASS)
//                  .build();
//              algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, "PAR_SASS", problemList.get(i), run));
//      }
      for (int i = 0; i < problemList.size(); i++) {
    	  Algorithm<List<DoubleSolution>> algorithm = new WASFGABuilder<>(problemList.get(i), new SBXCrossover(0.9, 20.0),
                  new PolynomialMutation(1.0 / problemList.get(i).getNumberOfVariables(), 20.0),configureReferencePoint(refPoint[i]))
                  .setMaxEvaluations(500)
                  .setPopulationSize(100)
                  .setVariant(WASFGABuilder.WASFGAVariant.WASFGA_ad)
                  .build();
              algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, "WASFGA_ad", problemList.get(i), run));
      }
//      for (int i = 0; i < problemList.size(); i++) {
//    	  Algorithm<List<DoubleSolution>> algorithm = new WASFGABuilder<>(problemList.get(i), new SBXCrossover(0.9, 20.0),
//                  new PolynomialMutation(1.0 / problemList.get(i).getNumberOfVariables(), 20.0),configureReferencePoint(refPoint[i]))
//                  .setMaxEvaluations(500)
//                  .setPopulationSize(100)
//                  .setVariant(WASFGABuilder.WASFGAVariant.NDrankWASFGA)
//                  .build();
//              algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, "NDrankWASFGA", problemList.get(i), run));
//      }
    }
    return algorithms ;
  }
  
  static List<Double> configureReferencePoint(double[] rp){
	  List<Double> referencePoint = new ArrayList<>();
	  for (int i=0; i<rp.length; i++){
		  referencePoint.add(rp[i]);
	  }
	  return referencePoint;
  }
		  
  
}
