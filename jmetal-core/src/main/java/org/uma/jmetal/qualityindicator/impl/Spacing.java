//  GeneralizedSpread.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
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

package org.uma.jmetal.qualityindicator.impl;

import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.front.Front;
import org.uma.jmetal.util.front.imp.ArrayFront;
import org.uma.jmetal.util.front.util.FrontUtils;
import org.uma.jmetal.util.point.Point;
import org.uma.jmetal.util.point.impl.ArrayPoint;
import org.uma.jmetal.util.point.util.comparator.LexicographicalPointComparator;
import org.uma.jmetal.util.point.util.comparator.PointDimensionComparator;
import org.uma.jmetal.util.point.util.distance.EuclideanDistance;

import java.io.FileNotFoundException;
import java.util.List;

/**
 * This class implements the spacing indicator.
 * Reference: 
 * 
 * 
 *
 * @author Miriam R. Ferrandez <mrferrandez@ual.es>
 * @author Juana L. Redondo
 */
@SuppressWarnings("serial")
public class Spacing<S extends Solution<?>> extends GenericIndicator<S> {

  /**
   * Default constructor
   */
  public Spacing() {
  }

  /**
   * Constructor
   *
   * @param referenceParetoFrontFile
   * @throws FileNotFoundException
   */
  public Spacing(String referenceParetoFrontFile) throws FileNotFoundException {
    super(referenceParetoFrontFile) ;
  }

  /**
   * Constructor
   *
   * @param referenceParetoFront
   */
  public Spacing(Front referenceParetoFront) {
    super(referenceParetoFront) ;
  }

  /**
   * Evaluate() method
   * @param solutionList
   * @return
   */
  @Override public Double evaluate(List<S> solutionList) {
	  if (solutionList == null) {
	      throw new JMetalException("The pareto front approximation is null") ;
	    }
	  return spacing(new ArrayFront(solutionList), referenceParetoFront);
  }

  /**
   * Returns the spacing indicator value for a given front
   *
   * @param front           The front
   * @param referenceFront The reference pareto front
   */
  public double spacing(Front front, Front referenceFront) {
	  //int numberOfObjectives = front.getPoint(0).getNumberOfDimensions() ;

	  int numberOfPoints = front.getNumberOfPoints();

	  //if (new EuclideanDistance().compute(front.getPoint(0),front.getPoint(front.getNumberOfPoints() - 1)) == 0.0) {
	  //return 1.0;
	  //} else {
      double dmean = 0.0;

      for (int i = 0 ; i < front.getNumberOfPoints(); i++) {
        dmean += FrontUtils.distanceToNearestPoint(front.getPoint(i), front);
      }

      dmean = dmean / (numberOfPoints);

      double mean = 0.0;
      for (int i = 0; i < front.getNumberOfPoints(); i++) {
    	  mean += Math.pow(FrontUtils.distanceToNearestPoint(front.getPoint(i), front) -
            dmean,2);
      }
      mean = mean / Math.pow(dmean, 2);

      return Math.sqrt(mean / (numberOfPoints));
    //}
  }

  @Override public String getName() {
    return "TS" ;
  }

  @Override public String getDescription() {
    return "Spacing quality indicator" ;
  }

  @Override
  public boolean isTheLowerTheIndicatorValueTheBetter() {
    return true ;
  }
}

