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

package org.uma.jmetal.util.distance.impl;

import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.distance.Distance;

import java.util.ArrayList;
import java.util.List;

/**
 * Class for calculating the Euclidean distance between two {@link Solution} objects in objective space.
 *
 * @author <antonio@lcc.uma.es>
 */
public class HausdorffDistanceBetweenSetsInObjectiveSpace<S extends Solution<?>, L extends List<S>>
    implements Distance<L, L> {
  
  protected Distance<S,S> distance;

  public HausdorffDistanceBetweenSetsInObjectiveSpace(Distance<S,S> distance){
	  this.distance = distance;
  }
	
  @Override
  public double getDistance(L solutionList1, L solutionList2) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    double max1;
    double max2;
    double min1;
    double min2;
    double hd;

    for (int sol1 = 0; sol1 < solutionList1.size(); sol1++){
    	min1 = minDistanceToSet(solutionList1.get(sol1), solutionList2);
    	sum1 = sum1 + min1;
    }
    for (int sol2 = 0; sol2 < solutionList2.size(); sol2++){
    	min2 = minDistanceToSet(solutionList2.get(sol2), solutionList1);
    	sum2 = sum2 + min2;
    }
    max1 = maxDistanceInSet(solutionList1);
    max2 = maxDistanceInSet(solutionList2);
    
    hd = sum1/max1 + sum2/max2;
    
    return hd/2.0;
  }
  
  public double minDistanceToSet(S solution, L solutionList) {
	  double min = distance.getDistance(solution, solutionList.get(0));
	  double aux;
	  for (int nsol = 1; nsol < solutionList.size(); nsol++){
		  aux = distance.getDistance(solution, solutionList.get(nsol));
		  if (aux < min){
			  min = aux;
		  }
	  }
	  return min;
  }
  
  public double maxDistanceInSet(L solutionList){
	  double max = 0.0;
	  double aux;
	  for (int i = 0; i < solutionList.size()-1; i++){
		  for (int j = i+1; j < solutionList.size(); j++){
			  aux = distance.getDistance(solutionList.get(i), solutionList.get(j));
			  if (aux > max){
				  max = aux;
			  }
		  }	
	  }	  
	  return max;
  }
	  
}
