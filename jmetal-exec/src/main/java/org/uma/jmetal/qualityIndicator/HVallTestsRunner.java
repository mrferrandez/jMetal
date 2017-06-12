package org.uma.jmetal.qualityIndicator;

public class HVallTestsRunner {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String [] input = new String [4];
		input[0]="HV";
		input[1]=args[0];
		input[3]="TRUE";
		for (int i=1; i<args.length; i++){
			input[2] = args[i];
			try {
				CommandLineIndicatorRunner.main(input);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

}
