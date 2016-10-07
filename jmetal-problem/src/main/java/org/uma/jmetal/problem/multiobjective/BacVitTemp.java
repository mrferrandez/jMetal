package org.uma.jmetal.problem.multiobjective;


import java.io.IOException;
import java.util.Arrays;

import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;
//import com.comsol.model.internal.Model;

import org.uma.jmetal.problem.ConstrainedProblem;
import org.uma.jmetal.problem.impl.AbstractDoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.solutionattribute.impl.NumberOfViolatedConstraints;
import org.uma.jmetal.util.solutionattribute.impl.OverallConstraintViolation;

import java.util.Arrays;
import java.util.List;

/**
* Class representing problem BacVit
*/
@SuppressWarnings("serial")
public class BacVitTemp extends AbstractDoubleProblem {

static String currentFolder = System.getProperty("user.dir");

/**
* Constructor.
* Creates a default instance of the LocationProblem.
* @param solutionType The solution type must "Real" 
*/
public BacVitTemp (){
    setNumberOfVariables(10);
    setNumberOfObjectives(3);
    setNumberOfConstraints(0);
    setName("BacVitTemp");

    List<Double> lowerLimit = Arrays.asList(10.0, 10.0, 0.1, -250.0, -250.0, -250.0, -250.0, -250.0, -250.0, -250.0);
    List<Double> upperLimit = Arrays.asList(50.0, 50.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0);

    setLowerLimit(lowerLimit);
    setUpperLimit(upperLimit);
    
} // BacVitTemp


/** Evaluate() method */
@Override
public void evaluate(DoubleSolution solution) {
	double [] xv = new double[getNumberOfVariables()] ; // 10 decision variables
	double [] f = new double[getNumberOfObjectives()] ; // 3 objectives

	for (int var = 0; var < xv.length; var++){
		xv[var] = solution.getVariableValue(var);
	}
	
	// Calcular el punto en el que en realidad se esta evaluando
	double[] coord=getcoord(xv);
	for (int var = 1; var < xv.length; var++){
		solution.setVariableValue(var, coord[var]);
	}

	ModelUtil.initStandalone(false);
	String[] xvt=tempini(xv);
	String[][] x=tablepres(xv);

	double[] vectsol = run(xvt,x);
	double Bac = vectsol[0];
	double Vit = vectsol[1];
	double Tmax = vectsol[2];

	// First objective
	f[0] = Bac;

	// Second objective
	f[1] = (1-Vit);// + 1.e9*Math.max(Tmax-50.0,0.0); 
	
	// Third objective
	f[2] = Tmax;

	solution.setObjective(0,f[0]);
	solution.setObjective(1,f[1]);
	solution.setObjective(2,f[2]);
	ModelUtil.disconnect();
	return ;

} // evaluate

public static String[] tempini(double[] xv) {
	String[] xvt={Double.toString(xv[0]), Double.toString(xv[1])};
	//System.out.println("Temperatura a tiempo 0 - inicial/refrigeracion:" +xvt[0]+" "+xvt[1]);
	return xvt;
}
public static double[] getcoord(double[] xv) {
	double[] coord=xv;
	double pres=coord[2];
	double pmin=0.1;
	double pmax=900.0;
	for (int ii=3; ii<coord.length; ii++){//ii=4:length(coord)
	    if(pres+coord[ii]<pmin){
	        coord[ii]=pmin-pres;
	    }
	    else{
	    	if (pres+coord[ii]>pmax){
	    		coord[ii]=pmax-pres;
	    	}  
	    }
	    pres=pres+coord[ii];
	}
	return coord;
}
public static String[][] tablepres(double[] xv) {
	
	xv=	Arrays.copyOfRange(xv, 2, xv.length);//xv(3:end);
	double incmax=2.5;
	double pmin=0.1;
	double pmax=900.0;
	double drt=900.0;
	double interv=drt/xv.length;
	double[] xp= new double[xv.length];
	xp[0]=Math.min(Math.max(xv[0],pmin),pmax);
	for (int ii=1; ii<xv.length; ii++){
	    xp[ii]=Math.min(Math.max(xp[ii-1]+xv[ii],pmin),pmax);
	}
		
	String[][] x;
	if ((xp[xv.length-1]-0.1)/incmax>1e-9){
		x= new String[xv.length+2][2];
		x[xv.length+1][0]=Double.toString(drt+(xp[xv.length-1]-0.1)/incmax);
		x[xv.length+1][1]="0.1";
	}
	else{
		x= new String[xv.length+1][2];
	}
	x[0][0]="0";
	x[0][1]="0.1";
	for (int ii=1; ii<xv.length+1; ii++){//ii=1:length(xv)
	    x[ii][0]=Double.toString(interv*ii);
	    x[ii][1]=Double.toString(xp[ii-1]);
	}
	
	
	//System.out.println("Curva de Presion: Tiempo (s) | Presion (MPa)");
	//for (int ii=0; ii<x.length; ii++){
		//System.out.println(x[ii][0]+" "+x[ii][1]);
	//}
		
	return x;
}

public static double[] run(String[] xvt, String[][] x) {//double, 
	
	Model model;
	try {
		model = ModelUtil.loadCopy("root", "Model");
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
		model = ModelUtil.create("Model");
	}
	
	model.param().set("Trefrig", xvt[1]);
	model.param().set("Tini", xvt[0]);
	model.func("int1").set("table", x);
	model.func("int1").set("funcname", "presv");

	model.study().create("std1");
	model.study("std1").create("time", "Transient");

	String timefin=x[x.length-1][0];

	model.study("std1").feature("time").set("tlist", "range(0,10,"+timefin+")");
	model.study("std1").feature("time").set("rtolactive", true);
	model.study("std1").feature("time").set("mesh", new String[]{"geom1", "mesh1"});

	model.sol().create("sol1");
	model.sol("sol1").study("std1");
	model.sol("sol1").feature().create("st1", "StudyStep");
	model.sol("sol1").feature("st1").set("study", "std1");
	model.sol("sol1").feature("st1").set("studystep", "time");
	model.sol("sol1").feature().create("v1", "Variables");
	model.sol("sol1").feature("v1").set("control", "time");

	model.shape("shape1").feature("shfun1");
	
	model.sol("sol1").feature().create("t1", "Time");
	model.sol("sol1").feature("t1").set("tlist", "range(0,10,"+timefin+")");
	model.sol("sol1").feature("t1").set("plot", "off");
	model.sol("sol1").feature("t1").set("plotfreq", "tout");
	model.sol("sol1").feature("t1").set("probesel", "all");
	model.sol("sol1").feature("t1").set("probefreq", "tsteps");
	model.sol("sol1").feature("t1").set("atolglobalmethod", "scaled");
	model.sol("sol1").feature("t1").set("atolglobal", 0.0010);
	model.sol("sol1").feature("t1").set("maxorder", 2);
	model.sol("sol1").feature("t1").set("control", "time");
	model.sol("sol1").feature("t1").feature().create("seDef", "Segregated");
	model.sol("sol1").feature("t1").feature().create("fc1", "FullyCoupled");
	model.sol("sol1").feature("t1").feature().create("d1", "Direct");
	model.sol("sol1").feature("t1").feature("fc1").set("linsolver", "d1");
	model.sol("sol1").feature("t1").feature().remove("fcDef");
	model.sol("sol1").feature("t1").feature().remove("seDef");
	model.sol("sol1").attach("std1");
	
	model.sol("sol1").runAll();

	model.result().numerical().create("av1", "AvSurface");
	model.result().numerical("av1").selection().set(new int[]{2});
	model.result().numerical("av1").set("expr", "Bac");
	model.result().numerical().create("av2", "AvSurface");
	model.result().numerical("av2").selection().set(new int[]{2});
	model.result().numerical("av2").set("expr", "Vit");
	
	model.result().table().create("tbl1", "Table");
	model.result().table("tbl1").comments("Superficie promedio 1 (Bac)");
	model.result().numerical("av1").set("table", "tbl1");
	model.result().numerical("av1").setResult();
	model.result().table().create("tbl2", "Table");
	model.result().table("tbl2").comments("Superficie promedio 2 (Vit)");
	model.result().numerical("av2").set("table", "tbl2");
	model.result().numerical("av2").setResult();
	
	model.result().table("tbl1").label("Bac");
	model.result().table("tbl2").label("Vit");
	
	model.result().export().create("tbl1", "Table");
	model.result().export().create("tbl2", "Table");
	model.result().export("tbl1").set("filename", "rBac");
	//model.result().export("tbl1").set("table", "tbl2");
	model.result().export("tbl2").set("filename", "rVit");
	model.result().create("pg8", "PlotGroup1D");
	model.result("pg8").create("tblp1", "Table");
	model.result("pg8").create("tblp2", "Table");
	model.result("pg8").feature("tblp2").set("table", "tbl2");
	
	model.result().numerical().create("max1", "MaxSurface");
	model.result().numerical("max1").selection().set(new int[]{2});
	model.result().table().create("tbl3", "Table");
	model.result().numerical("max1").set("table", "tbl3");
	model.result().numerical("max1").setResult();
	
	double[][] TBac=model.result().table("tbl1").getReal(); 
	double[][] TVit=model.result().table("tbl2").getReal(); 
	double[][] TTemp=model.result().table("tbl3").getReal(); 
	double Bac=TBac[TBac.length-1][1];
	double Vit=TVit[TVit.length-1][1];
	double Temp=TTemp[0][1];
	for (int ii=1; ii<TTemp.length; ii++){
		if(TTemp[ii][1]>Temp){
			Temp=TTemp[ii][1];
		}
	}
	double[] BacVitTemp={Bac, Vit, Temp};

	return BacVitTemp;
	}

} // 
