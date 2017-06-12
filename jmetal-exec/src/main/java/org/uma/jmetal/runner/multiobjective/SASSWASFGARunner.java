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

package org.uma.jmetal.runner.multiobjective;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.wasfga.SASSWASFGA;
import org.uma.jmetal.operator.impl.localsearch.SASSLocalSearch;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.runner.AbstractAlgorithmRunner;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.ProblemUtils;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

public class SASSWASFGARunner extends AbstractAlgorithmRunner {
  /**
   * @param args Command line arguments.
   * @throws JMetalException
   * @throws FileNotFoundException
   * Invoking command:
  java org.uma.jmetal.runner.multiobjective.NSGAII45Runner problemName [referenceFront]
   */
  public static void main(String[] args) throws JMetalException, FileNotFoundException {
    Problem<DoubleSolution> problem;
    Algorithm<List<DoubleSolution>> algorithm;
    String referenceParetoFront = "" ;
    List<Double> referencePoint = null;

    String problemName ;
//    if (args.length == 1) {
//      problemName = args[0];
//    } else if (args.length == 2) {
//      problemName = args[0] ;
//      referenceParetoFront = args[1] ;
//    } else {
//      problemName = "org.uma.jmetal.problem.multiobjective.zdt.ZDT4";
//      referenceParetoFront = "jmetal-problem/src/test/resources/pareto_fronts/ZDT4.pf" ;
//    }

    problemName ="org.uma.jmetal.problem.multiobjective.BacVitTemp";
    problem = ProblemUtils.<DoubleSolution> loadProblem(problemName);
    
    referencePoint = new ArrayList<>();
    referencePoint.add(0.0);
    referencePoint.add(0.0);
    referencePoint.add(30.0);
    
    int improvementRounds = 7;

    
    SASSLocalSearch<DoubleSolution> localSearch = new SASSLocalSearch<DoubleSolution>(improvementRounds, problem);
    
    algorithm = new SASSWASFGA<DoubleSolution,List<DoubleSolution>>(problem, 100, 3, localSearch,new SequentialSolutionListEvaluator<DoubleSolution>(),referencePoint) ;

    AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
            .execute() ;

    List<DoubleSolution> population = algorithm.getResult() ;
    long computingTime = algorithmRunner.getComputingTime() ;
    @SuppressWarnings("unchecked")
	List<DoubleSolution> allpopulation = ((SASSWASFGA<DoubleSolution,List<DoubleSolution>>) algorithm).getAllPopulation() ;

    JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

    printFinalSolutionSet(population); // only non dominated solutions are printed in FUN.tsv
    // Print all final population
    new SolutionListOutput(allpopulation)
    .setSeparator("\t")
    .setVarFileOutputContext(new DefaultFileOutputContext("VARall.tsv"))
    .setFunFileOutputContext(new DefaultFileOutputContext("FUNall.tsv"))
    .print();
    
    if (!referenceParetoFront.equals("")) {
      printQualityIndicators(population, referenceParetoFront) ;
    }
  }
}
